% TDSolver    EREBUS subroutine to solve the time-dependent sub-problem
%
% [CTX]  =  TDSolver(CTX)
%
%   Uses a finite-difference scheme collocated with the Q2 FE mesh to solve
%   advection-diffusion-reaction equations for the time-dependent evolution
%   of temperature, vesicularity and crystallinity.
%
%   created   20140806  Tobias Keller
%   modified  20170427  Tobias Keller
%   modified  20200227  Tobias Keller
%   modified  20200515  Tobias Keller


function  [CTX] = TDSolver(CTX)

% prepare fields
MP       =  CTX.MP;
SL       =  CTX.SL;
SLo      =  CTX.SLo;
FE       =  CTX.FE;
FEo      =  CTX.FEo;
PROP     =  CTX.PROP;
PHYS     =  CTX.PHYS;

% extract constant material properties
Tsol     =  PROP.Tsol;
Tliq     =  PROP.Tliq;

% prepare velocity field corrected for surface-deformed mesh
Vel      = [SL .U,SL .W] - FE .Vel;
Velo     = [SLo.U,SLo.W] - FEo.Vel;

% set boundary conditions for temperature T, bubble fraction phi
CTX.BC.SlugCount  =  floor((CTX.TIME.total+CTX.PHYS.TauIn/2) / (CTX.PHYS.SlugNo*CTX.PHYS.TauIn));
if  CTX.TIME.total >= (CTX.BC.SlugCount*CTX.PHYS.SlugNo-0.5)*(CTX.PHYS.TauIn) && ...
    CTX.TIME.total <  (CTX.BC.SlugCount*CTX.PHYS.SlugNo+0.5)*(CTX.PHYS.TauIn)
    CTX.BC.SlugOn      =  1;
else
    CTX.BC.SlugOn      =  0;
end
bcTmp = {'zf','zf','ct','ct'};
bcPhi = {'zf','zf','ct','ct'};
Surf  =  repmat(FE.SurfQ2.',FE.nzQ2,1);  % surface position
N     =  CTX.SL.RFN;                     % refined time step count

% refined time-stepping for advection-diffusion problem
for n = 1:N
    
    SLo      =  SL;  % store previous solution

    % prepare non-linear material properties
    RhoM     =  PIPQ2(MP.Rho   ,FE);  % magma density

    cL       =  PROP.c + (1-SL.Phi).*PROP.L./(Tliq-Tsol);  % adjusted heat capacity
    ind      =  SL.T < Tsol;
    cL(ind)  =  PROP.c;
    
    kTmp     =  PROP.k./RhoM./cL;  % thermal diffusivity
    
    % compressible term
    Cmpr     = -(1-SL.Phi).*SL.Phi./PQ1Q2(SL.Pt,FE).*RhoM.*CTX.PHYS.grav.*SL.W .* CTX.PROP.Compr;

    % adiabatic cooling term
    Adbt     =  CTX.PROP.alpha.*(SL.T+273.15)./cL.*CTX.PHYS.grav.*SL.W .* CTX.PROP.Compr;
    
    % viscous dissipation term
    Diss     =  max(0,PIPQ2(MP.Diss,FE)./RhoM./cL);
    
    % approximate 3D heat loss
    HeatLoss =  -2.*kTmp.*(SL.T - CTX.INIT.TempExt)./(CTX.INIT.MatWidth/2+PROP.CntAur).^2;
    
    % heat input and extraction source terms
    TauTmp   =  PHYS.TauTmp.*max(1,PIPQ2((MP.EtaVP./CTX.PROP.Eta(2)).^0.5,FE));
    HeatIn   =  -min(0,SL.T-CTX.BC  .Tmp0   )./(FE.hzQ2/PHYS.WIn) .* exp( (FE.CoordQ2(:,2)-FE.D   )./PHYS.delta) .* CTX.BC.SlugOn;
    HeatOut  =  -   (  SL.T-CTX.INIT.TempExt)./TauTmp             .* exp(-(FE.CoordQ2(:,2)-Surf(:))./PHYS.delta);
    
    % heat evolution source term
    SL.GTmp  =  HeatIn + HeatOut + Adbt + Diss + HeatLoss;
    SL.HeatIn=  HeatIn;  SL.HeatOut = HeatOut;
    
    % solve heat evolution equation
    SL.T     =  AdvDiffSolver(SL.T,SLo.T,Vel,Velo,SL.GTmp,SLo.GTmp,kTmp,CTX.TIME.step/N,'Q2',bcTmp,CTX);
  
    % update crystallinity at new temperature
    SL.Chi   =  max(0,min(1,(SL.T-PROP.Tliq)./(PROP.Tsol-PROP.Tliq))).*(1-SL.Phi);
    MP.Chi   =  max(0,min(1,PQ2IP(SL.Chi,CTX.FE)));
    
    % Bubble input and extraction source terms
    TauPhi   =  PHYS.TauPhi.*max(1,PIPQ2((MP.EtaVP./CTX.PROP.Eta(2)).^0.5,FE));
    GasIn    = -min(0,SL.Phi-CTX.BC.Phi0)./(FE.hzQ2/PHYS.WIn) .* exp( (FE.CoordQ2(:,2)-FE.D   )./PHYS.delta) .* CTX.BC.SlugOn;
    GasOut   = -   (  SL.Phi-0          )./TauPhi             .* exp(-(FE.CoordQ2(:,2)-Surf(:))./PHYS.delta);
    
    % vesicularity evolution source term
    SL.GPhi  =  GasIn + GasOut + Cmpr;
    SL.GasIn =  GasIn;  SL.GasOut = GasOut;

    kPhi     =  PROP.kPhi./max(1,PIPQ2((MP.EtaVP./CTX.PROP.Eta(2)).^0.5,FE));  % thermal diffusivity

    % solve vesicularity evolution equation
    SL.Phi   =  AdvDiffSolver(SL.Phi,SLo.Phi,Vel,Velo,SL.GPhi,SLo.GPhi,kPhi,CTX.TIME.step/N,'Q2',bcPhi,CTX);
    SL.Phi   =  max(0,min(1,SL.Phi));
    MP.Phi   =  max(0,min(1,PQ2IP(CTX.SL.Phi,CTX.FE)));
    
    % update density at new T, phi, chi
    MP       =  Density(MP,CTX);

end

CTX.SL  =  SL;
CTX.MP  =  MP;
