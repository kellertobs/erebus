% InitMat    EDIFICE: Set initial distribution of material types
%
% [CTX] = InitMat(CTX)
%
%   Function initializes the distribution of material types with distinct
%   properties throughout domain
%
%   created   20140730  Tobias Keller
%   modified  20170427  Tobias Keller
%   modified  20200227   Tobias Keller


function  [CTX]  =  InitMat(CTX)


%***  prepare structures and variables
FE    =  CTX.FE;
INIT  =  CTX.INIT;
PROP  =  CTX.PROP;
COORD =  CTX.FE.CoordIP;
        
        
%*****  initialize material type on material points  **********************
MP.Mat  =  zeros(FE.NIP,1);


%***  constant composition
if     strcmp(INIT.MatMode(1:5),'const')
    
    MP.Mat(:)  =  INIT.Mat(1);

    
%***  two-layer model (w/ initial perturbation)
elseif strcmp(INIT.MatMode(1:5),'layer')
     
    ind          =  COORD(:,2) >= INIT.MatThick;
    MP.Mat(ind)  =  INIT.Mat(1);
    ind          =  COORD(:,2) <  INIT.MatThick;
    MP.Mat(ind)  =  INIT.Mat(2);
    
%***  three-layer model upper/lower crust over lithosphere   
elseif strcmp(INIT.MatMode(1:5),'crust')
    
    ind          =  COORD(:,2) >= INIT.Moho;
    MP.Mat(ind)  =  INIT.Mat(1);
    ind          =  COORD(:,2) <  INIT.Moho & COORD(:,2) <= INIT.Conrad;
    MP.Mat(ind)  =  INIT.Mat(2);
    ind          =  COORD(:,2) <  INIT.Conrad;
    MP.Mat(ind)  =  INIT.Mat(3);
    
%***  multi-layer stack over detachment layer
elseif strcmp(INIT.MatMode(1:5),'multi')

    h0  =  INIT.MatThick;
    n   =  INIT.MatLayers;
    dh  =  h0/n;
    
    ind          =  COORD(:,2) >= h0;
    MP.Mat(ind)  =  INIT.Mat(1);
    
    for i = 1:2:n
        ind          =  COORD(:,2) < h0-(i-1)*dh & COORD(:,2) >= h0- i   *dh;
        MP.Mat(ind)  =  INIT.Mat(2);
        ind          =  COORD(:,2) < h0- i   *dh & COORD(:,2) >= h0-(i+1)*dh;
        MP.Mat(ind)  =  INIT.Mat(3);
    end
    
%***  idealised lava lake bed
elseif strcmp(INIT.MatMode(1:4),'lava') % 'lavalake'
    
    MP.Mat(:)    =  INIT.Mat(1);
    redge        =  max(INIT.MatXLoc + INIT.MatWidth/2, FE.W - (FE.W-INIT.MatXLoc-INIT.MatWidth/2)./INIT.MatHeight.*COORD(:,2));
    ledge        =  min(INIT.MatXLoc - INIT.MatWidth/2, 0    + (INIT.MatXLoc-INIT.MatWidth/2)./INIT.MatHeight.*COORD(:,2));
    lake         =  COORD(:,1) >= ledge & COORD(:,1) <= redge;
    MP.Mat(lake) =  INIT.Mat(2);
    
end


%***  add block of specified composition
if INIT.AddBlock > 0
    for n = 1:INIT.AddBlock
        ind          =  abs(COORD(:,1)-INIT.BlockXLoc(n)) <= INIT.BlockWidth( n)/2 & ...
                        abs(COORD(:,2)-INIT.BlockZLoc(n)) <= INIT.BlockHeight(n)/2 ;
        MP.Mat(ind)  =  INIT.BlockMat(n);
    end
end
    

%***  add sphere of specified composition
if INIT.AddSphere > 0
    for n = 1:INIT.AddSphere
        ind          =  (COORD(:,1)-INIT.SphereXLoc(n)).^2./INIT.SphereRadiusX(n)^2  ...
                     +  (COORD(:,2)-INIT.SphereZLoc(n)).^2./INIT.SphereRadiusZ(n)^2 <= 1;
        MP.Mat(ind)  =  INIT.SphereMat(n);
    end
end


%*****  initialize material properties and auxiliary fields  **************
MP.Plith      =  PROP.Rho.*CTX.PHYS.grav.*FE.CoordIP(:,2) + CTX.BC.SurfPres;
MP.Eta        =  PROP.Eta(MP.Mat);
MP.EtaVP      =  PROP.Eta(MP.Mat);
MP.EtaVEP     =  PROP.Eta(MP.Mat);
MP.Xi         =  zeros(FE.NIP,1);
MP.Gdt        =  PROP.G(MP.Mat).*CTX.TIME.step;
MP.Edot       =  zeros(FE.NIP,3);
MP.Edot(:,1)  =  ones(FE.NIP,1).*CTX.RHEO.Strainr0;
MP.Edot(:,2)  = -ones(FE.NIP,1).*CTX.RHEO.Strainr0;
MP.Edot(:,3)  =  1e-24;
MP.Tau        =  zeros(FE.NIP,3);
MP.Taur       =  MP.Tau;
MP.EII        =  zeros(FE.NIP,4);
MP.EII(:,1)   =  SecondInv(MP.Edot);
MP.EII(:,2)   =  SecondInv(MP.Edot);
MP.TII        =  SecondInv(MP.Tau);
MP.TIIr       =  SecondInv(MP.Taur);
MP.YieldStr   =  PROP.Coh(MP.Mat);
MP.Dmg        =  zeros(FE.NIP,1);
MP.Rho        =  PROP.Rho.*ones(size(MP.Mat));

if ~isfield(CTX.SL,'RhoRef')
    if strcmp(CTX.BC.FreeSurface,'OFF')
        CTX.SL.RhoRef  =  mean(CTX.MP.Rho);
    else
        CTX.SL.RhoRef  =  0;
    end
end

%***  return structures
CTX.MP  =  MP;


end

