
function  [CTX] = TCSolver(CTX)

MP       =  CTX.MP;
MPo      =  CTX.MPo;
SL       =  CTX.SL;
SLo      =  CTX.SLo;
FE       =  CTX.FE;
FEo      =  CTX.FEo;
PROP     =  CTX.PROP;
PHYS     =  CTX.PHYS;

Tsol     =  PROP.Tsol;
Tliq     =  PROP.Tliq;
LatH     =  PROP.L;
Hcap     =  PROP.c;
kTmp     =  PROP.k;
kPhi     =  PROP.kPhi.*ones(size(SL.Phi));

Vel      = [SL .U,SL .W] - FE .Vel;
Velo     = [SLo.U,SLo.W] - FEo.Vel;

bcTmp = {[],CTX.BC.botTmp,CTX.INIT.TempExt,CTX.INIT.TempExt};
bcPhi = {[],CTX.BC.botPhi,0,0};

Surf  =  repmat(FE.SurfQ2.',FE.nzQ2,1);
N     =  CTX.SL.RFN;

for n = 1:N
    
    RhoM     =  PIPQ2(MP.Rho,FE);
    
    cL       =  Hcap + (1-SL.Phi).*LatH./(Tliq-Tsol);
    ind      =  SL.T < Tsol;
    cL(ind)  =  Hcap;
    kTmp     =  kTmp./RhoM./cL;
    
    Cmpr     = -(1-SL.Phi).*SL.Phi.*RhoM.*CTX.PHYS.grav./PQ1Q2(SL.Pt,FE).*SL.W;

    Adbt     =  1e-4.*(SL.T+273.15)./cL.*CTX.PHYS.grav.*SL.W;
    Diss     =  max(0,PIPQ2(MP.Diss,FE)./RhoM./cL);
    
    Tmp0     =  CTX.INIT.TempExt;
    Phi0     =  0;
    
    SL.GTmp  =  -(SL .T  -Tmp0)./PHYS.TauTmp .* exp(-(FE.CoordQ2(:,2)-Surf(:))./PHYS.delta) + Adbt + Diss;
    
    SL.GPhi  =  -(SL .Phi-Phi0)./PHYS.TauPhi .* exp(-(FE.CoordQ2(:,2)-Surf(:))./PHYS.delta) + Cmpr;
    
    SL.T     =  ExplicitAdvDiffSolver(SL.T,SLo.T,Vel,Velo,SL.GTmp,SLo.GTmp,kTmp,CTX.TIME.step/N,'Q2',bcTmp,CTX);
    SL.Phi   =  ExplicitAdvDiffSolver(SL.Phi,SLo.Phi,Vel,Velo,SL.GPhi,SLo.GPhi,kPhi,CTX.TIME.step/N,'Q2',bcPhi,CTX);
    SL.Phi   =  max(0,min(1,SL.Phi));
    
    MP.Phi   =  max(0,min(1,PQ2IP(CTX.SL.Phi,CTX.FE)));
    
    SL.Chi   =  max(0,min(1,(SL.T-PROP.Tliq)./(PROP.Tsol-PROP.Tliq))).*(1-SL.Phi);
    MP.Chi   =  max(0,min(1,PQ2IP(SL.Chi,CTX.FE)));

    MP       =  Density(MP,CTX);
    SLo      =  SL;

end

CTX.SL  =  SL;
CTX.MP  =  MP;
