% Rheology    EREBUS subroutine to update material rheology
%
% [MP] = Rheology(MP,CTX)
%
%   Function updates visco-elasto-plastic material rheology according to
%   the latest solution guess. Parameters and options for rheology read
%   from CTX.RHEO struct.
%
%   created   20140730  Tobias Keller
%   modified  20170427  Tobias Keller
%   modified  20170508  Tobias Keller
%   modified  20200227  Tobias Keller
%   modified  20200515  Tobias Keller


function   [MP]  =  Rheology(MP,CTX)

RHEO   =  CTX.RHEO;
PROP   =  CTX.PROP;

Mat    =  CTX.MP.Mat;
EII    =  CTX.MP.EII;


%*****  get shear viscosity  eta  *****************************************

[MP]  =  Viscosity(EII,Mat,PROP,MP,RHEO,CTX);


%*****  apply elasto-plasticity to get effective eta  *********************

[MP]  =  ElastoPlasticity(EII,Mat,PROP,MP,RHEO,CTX);


end  % end function



%*****  Calculate P,T,c,str-dependent Shear Viscosity  ********************

function [MP]  =  Viscosity(EII,Mat,PROP,MP,RHEO,CTX)

TINY  =  1e-32;
thit  =  CTX.SL.theta_it;

EII0  =  RHEO.Strainr0;
n     =  (1-RHEO.PowerLawExp)/RHEO.PowerLawExp;

Eta0  =  PROP.Eta(Mat);
Eta0  =  10.^SmoothField(log10(Eta0),1,10,CTX.FE,'IP');

if strcmp(RHEO.ConstViscosity,'ON')
    MP.Eta  =  Eta0;
else
    
    %*****  compute stiffening/softening by crystals and bubbles  *********
    fact = max(1e-6,1-MP.Chi./RHEO.Chic).^-RHEO.mChi .* max(1e-6,1-MP.Phi./RHEO.Phic).^-RHEO.mPhi;
    Eta0 = Eta0 .* fact;
    
    %*****  compute strainrate-dependence of powerlaw viscosity  **********
    
    EIIVP   =  max(TINY,EII(:,1) - EII(:,3));
    EtaNN   =  Eta0 .* (EIIVP./EII0).^n;
    Eta     =  2./(1./Eta0 + 1./EtaNN);
    Eta     =  Eta.^thit .* MP.Eta.^(1-thit);

end   % end if-condition

MP.Eta  =  max(RHEO.MinEta,min(RHEO.MaxEta,Eta));

end  % end function



%*****  Calculate Elasto-Plastic Weakening of Viscosity  ******************


function [MP]  =  ElastoPlasticity(EII,Mat,PROP,MP,RHEO,CTX)

HUGE  =  1e+32;
TINY  =  1e-32;
lim   =  1e-3.*mean(PROP.Coh);
thit  =  CTX.SL.theta_it;

if strcmp(RHEO.Elasticity,'ON')
    MP.Gdt     =  PROP.G(Mat).*CTX.TIME.step + RHEO.MinGdt;
else 
    MP.Gdt(:)  =  HUGE;
end


%*****  get shear and compaction failure criteria  ************************

if strcmp(RHEO.Plasticity,'ON')
        
    %***  get cohesion/friction coeff.
    Coh    =  PROP.Coh(Mat);
    Coh    =  Coh - MP.pert.*CTX.INIT.PertCoh;
    Frict  =  PROP.Frict(Mat);
    Frict  =  Frict - MP.pert.*CTX.INIT.PertFrict;
    
    %***  get yield stress
    if     strcmp(CTX.FE.ElType(3:4),'Q1')
        Pt  =  PQ1IP(CTX.SL.Pt,CTX.FE) - CTX.BC.SurfPres;
    elseif strcmp(CTX.FE.ElType(3:4),'P0')
        Pt  =  PQ1IP(PElQ1(CTX.SL.Pt,CTX.FE),CTX.FE) - CTX.BC.SurfPres;
    end
    YieldStr     =  (1-MP.Phi) .* (Coh + Frict.*Pt);
    MP.YieldStr  =  max(lim,YieldStr);
    
else
    
    MP.YieldStr   =  HUGE.*ones(size(MP.Plith));
    
end  % end plasticity condition


%*****  compute effective visco-elasto-plastic properties  ****************

EIIVP       =  max(TINY,EII(:,1) - EII(:,3));
EtaP        =  MP.YieldStr./EIIVP + RHEO.MinEta;
EtaVP       =  min(MP.Eta,EtaP);
EtaVP       =  10.^SmoothField(log10(EtaVP),CTX.SL.Smooth,1,CTX.FE,'IP');
MP.EtaVP    =  EtaVP.^thit .* MP.EtaVP.^(1-thit);

MP.EtaVEP   =  (1./MP.EtaVP + 1./MP.Gdt).^-1;
MP.Xi       =  (1 + MP.Gdt ./ MP.EtaVP ).^-1;


end  % end function


