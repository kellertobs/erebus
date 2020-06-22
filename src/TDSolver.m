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
LatH     =  PROP.L;
Hcap     =  PROP.c;
kTmp     =  PROP.k;
kPhi     =  PROP.kPhi.*ones(size(SL.Phi));

% prepare velocity field corrected for surface-deformed mesh
Vel      = [SL .U,SL .W] - FE .Vel;
Velo     = [SLo.U,SLo.W] - FEo.Vel;

% set boundary conditions for temperature T, bubble fraction phi
bcTmp = {'zf',CTX.BC.botTmp,CTX.INIT.TempExt,CTX.INIT.TempExt}; 
bcPhi = {'zf',CTX.BC.botPhi,0,0};

Surf  =  repmat(FE.SurfQ2.',FE.nzQ2,1);  % surface position
N     =  CTX.SL.RFN;                     % refined time step count

% refined time-stepping for advection-diffusion problem
for n = 1:N
    
    SLo      =  SL;  % store previous solution

    % prepare non-linear material properties
    RhoM     =  PIPQ2(MP.Rho,FE);  % magma density
    
    cL       =  Hcap + (1-SL.Phi).*LatH./(Tliq-Tsol);  % adjusted heat capacity
    ind      =  SL.T < Tsol;
    cL(ind)  =  Hcap;
    
    kTmp     =  kTmp./RhoM./cL;  % thermal diffusivity
    
    % compressible term
    Cmpr     = -(1-SL.Phi).*SL.Phi.*RhoM.*CTX.PHYS.grav./PQ1Q2(SL.Pt,FE).*SL.W .* CTX.PROP.Compr;
    
    % adiabatic cooling term
    Adbt     =  1e-4.*(SL.T+273.15)./cL.*CTX.PHYS.grav.*SL.W;
    
    % viscous dissipation term
    Diss     =  max(0,PIPQ2(MP.Diss,FE)./RhoM./cL);
    
    % surface layer target temperature and vesicularity
    Tmp0     =  CTX.INIT.TempExt;
    Phi0     =  0;
    
    % heat evolution source term
    SL.GTmp  =  -(SL .T  -Tmp0)./PHYS.TauTmp .* exp(-(FE.CoordQ2(:,2)-Surf(:))./PHYS.delta) + Adbt + Diss;
    
    % solve heat evolution equation
    SL.T     =  AdvDiffSolver(SL.T,SLo.T,Vel,Velo,SL.GTmp,SLo.GTmp,kTmp,CTX.TIME.step/N,'Q2',bcTmp,CTX);
  
    % update crystallinity at new temperature
    SL.Chi   =  max(0,min(1,(SL.T-PROP.Tliq)./(PROP.Tsol-PROP.Tliq))).*(1-SL.Phi);
    MP.Chi   =  max(0,min(1,PQ2IP(SL.Chi,CTX.FE)));
    
    % vesicularity evolution source term
    SL.GPhi  =  -(SL .Phi-Phi0)./PHYS.TauPhi .* exp(-(FE.CoordQ2(:,2)-Surf(:))./PHYS.delta) + Cmpr;

    % solve vesicularity evolution equation
    SL.Phi   =  AdvDiffSolver(SL.Phi,SLo.Phi,Vel,Velo,SL.GPhi,SLo.GPhi,kPhi,CTX.TIME.step/N,'Q2',bcPhi,CTX);
    SL.Phi   =  max(0,min(1,SL.Phi));
    MP.Phi   =  max(0,min(1,PQ2IP(CTX.SL.Phi,CTX.FE)));
    
    % update density at new T, phi, chi
    MP       =  Density(MP,CTX);

end

CTX.SL  =  SL;
CTX.MP  =  MP;
