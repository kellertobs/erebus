% ApplyBC    EDIFICE: Apply and extract boundary conditions
%
% [L,RHS,BC]  =  ApplyBC(L,RHS,X,BC,FE)
%
%   Function applies boundary conditions to Stokes problem according to
%   options specified in BC, and extracts them to reduce problem size.
%
%   created 20161115 Tobias Keller
%   modified  20170508  Tobias Keller
%   modified  20200227   Tobias Keller


function  [L,RHS,BC]  =  ApplyBC(L,RHS,X,CTX)

FE = CTX.FE;
BC = CTX.BC;


mapv           =  FE.MapV;
mapP           =  FE.MapP;

apply_topvx    =  0;
apply_topvz    =  0;
apply_botvx    =  0;
apply_botvz    =  0;
apply_leftvx   =  0;
apply_leftvz   =  0;
apply_rightvx  =  0;
apply_rightvz  =  0;
apply_Pfix     =  0;

if     strcmp(BC.TopBot,'ns')  % no slip
    apply_topvx  =  1;
    apply_topvz  =  1;
    apply_botvx  =  1;
    apply_botvz  =  1;
    
elseif strcmp(BC.TopBot,'fs')  % free slip
    apply_topvz  =  1;
    apply_botvz  =  1;
    
elseif strcmp(BC.TopBot,'ts')  % top slip only
    apply_topvz  =  1;
    apply_botvx  =  1;
    apply_botvz  =  1;
    
elseif strcmp(BC.TopBot,'bs')  % bottom slip only
    apply_topvx  =  1;
    apply_topvz  =  1;
    apply_botvz  =  1;
    
end


if     strcmp(BC.Sides,'ns')  % no slip
    apply_leftvx  =  1;
    apply_leftvz  =  1;
    apply_rightvx =  1;
    apply_rightvz =  1;
    
elseif strcmp(BC.Sides,'fs')  % free slip
    apply_leftvx  =  1;
    apply_rightvx =  1;
    
elseif strcmp(BC.Sides,'rs')  % right side only free slip
    apply_leftvx  =  1;
    apply_leftvz  =  1;
    apply_rightvx =  1;
    
elseif strcmp(BC.Sides,'ls')  % left side only free slip
    apply_leftvx  =  1;
    apply_rightvx =  1;
    apply_rightvz =  1;
    
end

if strcmp(BC.FreeSurface,'ON')
    apply_topvx    =  0;
    apply_topvz    =  0;
else
    apply_Pfix     =  1;
end

if strcmp(CTX.INIT.MatMode,'lavalake')
    x                  =  (FE.CoordU(FE.MapUn(1,:),1)' - CTX.INIT.MatXLoc) ./ (0.8*CTX.INIT.MatWidth);
    ind                =  x >= -0.5 & x <= 0.5;
    inflow             =  (cos(x(ind).*4*pi)+cos(x(ind).*2*pi))./2 .* - CTX.PHYS.InflowRate;
    inflow             =  inflow - mean(inflow);
    inflow             =  inflow .* (1+cos(CTX.TIME.total.*2*pi/CTX.PHYS.InflowPeriod))./2;
    BC.WTopBot(2,ind)  =  inflow + mean(BC.WTopBot(2,~ind));
    apply_botvz        =  1;
    apply_botvx        =  1;
end

%*****  optimized extracted boundary conditions  **************************

bc_ind   = [];
bc_val   = [];

if apply_topvx    ==  1
    ind      =  mapv(1,2:end-1,1)';
    bc_ind   =  [bc_ind;ind];
    tmp      =  ones(FE.nxU-2,1).*BC.UTopBot(1,2:end-1).'./diag(X(ind,ind));
    bc_val   =  [bc_val;tmp];
end

if apply_topvz    ==  1
    ind      =  mapv(1,:,2)';
    bc_ind   =  [bc_ind;ind];
    tmp      =  ones(FE.nxU,1).*BC.WTopBot(1,:).'./diag(X(ind,ind));
    bc_val   =  [bc_val;tmp];
end

if apply_botvx    ==  1
    ind      =  mapv(end,2:end-1,1)';
    bc_ind   =  [bc_ind;ind];
    tmp      =  ones(FE.nxU-2,1).*BC.UTopBot(2,2:end-1).'./diag(X(ind,ind));
    bc_val   =  [bc_val;tmp];
end

if apply_botvz    ==  1
    ind      =  mapv(end,:,2)';
    bc_ind   =  [bc_ind;ind];
    tmp      =  ones(FE.nxU,1).*BC.WTopBot(2,:).'./diag(X(ind,ind));
    bc_val   =  [bc_val;tmp];
end

if apply_leftvx   ==  1
    ind      =  mapv(:,1,1);
    bc_ind   =  [bc_ind;ind];
    tmp      =  ones(FE.nzU,1).*BC.USides(1,:).'./diag(X(ind,ind));
    bc_val   =  [bc_val;tmp];
end

if apply_leftvz   ==  1
    ind      =  mapv(2:end-1,1,2);
    bc_ind   =  [bc_ind;ind];
    tmp      =  ones(FE.nzU-2,1).*BC.WSides(1,2:end-1).'./diag(X(ind,ind));
    bc_val   =  [bc_val;tmp];
end

if apply_rightvx  ==  1
    ind      =  mapv(:,end,1);
    bc_ind   =  [bc_ind;ind];
    tmp      =  ones(FE.nzU,1).*BC.USides(2,:).'./diag(X(ind,ind));
    bc_val   =  [bc_val;tmp];
end

if apply_rightvz  ==  1
    ind      =  mapv(2:end-1,end,2);
    bc_ind   =  [bc_ind;ind];
    tmp      =  ones(FE.nzU-2,1).*BC.WSides(2,2:end-1).'./diag(X(ind,ind));
    bc_val   =  [bc_val;tmp];
end

% fix P magnitude if no free surface
if apply_Pfix == 1
    P0       =  0;
    ind      =  mapP(round(CTX.INIT.MatZLoc/FE.D*FE.nzP),1+round(CTX.INIT.MatXLoc/FE.W*FE.nxP));
    bc_ind   =  [bc_ind;ind(:)];
    tmp      =  P0.*ones(size(ind(:)))./diag(X(ind,ind));
    bc_val   =  [bc_val;tmp];
end

% don't compute velocity solution in rigid lake walls
if strcmp(CTX.INIT.MatMode,'lavalake')
    Mat = CTX.MP.Mat;
    Eta = CTX.MP.Eta;
    Max = 0.1.*CTX.RHEO.MaxEta;
    if strcmp(FE.ElType,'Q2Q1')
        indu     =  PIPQ2(Eta,FE)>=Max & PIPQ2(Mat,FE)<=1.05;
        indw     =  PIPQ2(Eta,FE)>=Max & PIPQ2(Mat,FE)<=1.05;
        indp     =  PIPQ1(Eta,FE)>=Max & PIPQ1(Mat,FE)<=1.05;
    elseif strcmp(FE.ElType,'Q1Q1')
        indu     =  PIPQ1(Eta,FE)>=Max & PIPQ1(Mat,FE)<=1.05;
        indw     =  PIPQ1(Eta,FE)>=Max & PIPQ1(Mat,FE)<=1.05;
        indp     =  PIPQ1(Eta,FE)>=Max & PIPQ1(Mat,FE)<=1.05;
    elseif strcmp(FE.ElType,'Q1P0')
        indu     =  PIPQ1(Eta,FE)>=Max & PIPQ1(Mat,FE)<=1.05;
        indw     =  PIPQ1(Eta,FE)>=Max & PIPQ1(Mat,FE)<=1.05;
        indp     =  PIPEl(Eta,FE)>=Max & PIPEl(Mat,FE)<=1.05;
    end
    bc_ind   =  [bc_ind;FE.DOFU(indu);FE.DOFW(indw);FE.DOFP(indp)];
    tmpu     =  zeros(size(FE.DOFU(indu)));
    tmpw     =  zeros(size(FE.DOFW(indw)));
    tmpp     =  (CTX.PROP.RhoChi-CTX.SL.RhoRef).*CTX.PHYS.grav.*FE.CoordP(indp,2)./diag(X(FE.DOFP(indp),FE.DOFP(indp)));
    bc_val   =  [bc_val;tmpu;tmpw;tmpp];
end

% assemble and sort all boundary indices and values
[BC.ind,ind]     =  sort(bc_ind);
BC.val           =  bc_val(ind);

% extract boundary conditions to reduce problem size
BC.free          =  1:length(L(1,:));
BC.free(BC.ind)  =  [];
TMP              =  L(:,BC.ind);
RHS              =  RHS - TMP*BC.val;
L                =  L(BC.free,BC.free);

end