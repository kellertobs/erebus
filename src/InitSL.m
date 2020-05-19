% InitSL    EREBUS subroutine to initialise solution variables
%
% [CTX] = InitSL(CTX)
%
%   Function initializes the solution variables for the Stokes flow  and
%   time-dependent sub-problems comprising velocities U, W, and pressure
%   P, and the temperature T, vesicularity Phi, and crystallinity Chi,
%   respectively.
%
%   created   20140730  Tobias Keller
%   modified  20170427  Tobias Keller
%   modified  20190418  Tobias Keller
%   modified  20200227  Tobias Keller
%   modified  20200515  Tobias Keller


function  [CTX]  =  InitSL(CTX)


%***  prepare structures and variables
SL   =  CTX.SL;
MP   =  CTX.MP;
FE   =  CTX.FE;
BC   =  CTX.BC;
INIT = CTX.INIT;
PROP = CTX.PROP;


%*****  initialize velocity field  ****************************************

if strcmp(BC.Type,'ConstStr')
    BC.USides(1,:)  =  BC.BGStrainr.*FE.W;
    BC.USides(2,:)  =  0;
end

dv               =  BC.USides(2) - BC.USides(1);
U                =  BC.USides(1) + repmat(((1:FE.nxU)-1)*dv/(FE.nxU-1),FE.nzU,1);
W                =  -repmat(((1:FE.nzU)'-1)*dv/FE.aspect/(FE.nzU-1),1,FE.nxU);
     
BC.UTopBot       =  zeros(2,FE.nxU);
BC.WTopBot       =  zeros(2,FE.nxU);
BC.USides        =  zeros(2,FE.nzU);
BC.WSides        =  zeros(2,FE.nzU);

BC.UTopBot(1,:)  =  U(1  ,:);
BC.UTopBot(2,:)  =  U(end,:);
BC.WTopBot(1,:)  =  W(1  ,:);
BC.WTopBot(2,:)  =  W(end,:);

BC.USides(1,:)   =  U(:,1  );
BC.USides(2,:)   =  U(:,end);
BC.WSides(1,:)   =  W(:,1  );
BC.WSides(2,:)   =  W(:,end);

SL.U             =  1e-6.*U(:);
SL.W             =  1e-6.*W(:);
    

%*****  initialize dynamic, lithostatic, total pressure fields  ***********

SL.Pref   =  SL.RhoRef.*CTX.PHYS.grav.*FE.CoordP(:,2) + CTX.BC.SurfPres;
SL.P      =  (PROP.Rho - SL.RhoRef).*CTX.PHYS.grav.*FE.CoordP(:,2);
SL.Pt     =  SL.P + SL.Pref;


%*****  initialise temperature field  *************************************

%***  idealised lava lake bed
if strcmp(INIT.TempMode(1:4),'lava') % 'lavalake'
    
    lake   =  PIPQ2(MP.Mat,FE);
    SL.T   =  INIT.TempExt + (lake>=1.5).*(INIT.TempInt-INIT.TempExt);
    SL.T   =  SmoothField(SL.T,1/8,round(20/CTX.FE.hzQ1^2),FE,'Q2');

end

BC.botTmp  =  SL.T(FE.MapQ2(end,:));

%*****  initialise bubble and crystal fractions  **************************
Phi0    =  [0;0;PROP.Phi0];
bubble  =  SmoothField(PIPQ2(Phi0(MP.Mat),FE),1/8,round(20/CTX.FE.hzQ1^2),FE,'Q2'); 
SL.Phi  =  bubble + PROP.Phi0.*PIPQ2((MP.pert+1)./2,FE).*(SL.T-INIT.TempExt)./(INIT.TempInt-INIT.TempExt).*(1-FE.CoordQ2(:,2)./FE.D);
SL.Phi  =  max(0,min(PROP.Phi0,SL.Phi));
MP.Phi  =  max(0,min(PROP.Phi0,PQ2IP(SL.Phi,CTX.FE)));

SL.GTmp    =  0.*SL.Phi;
SL.GPhi    =  0.*SL.Phi;
BC.botPhi  =  SL.Phi(FE.MapQ2(end,:));

SL.Chi  =  max(0,min(1,(SL.T-PROP.Tliq)./(PROP.Tsol-PROP.Tliq))).*(1-SL.Phi);
MP.Chi  =  max(0,min(1,PQ2IP(SL.Chi,CTX.FE)));


%*** hand back variables and structures
SL.U  =  SL.U(:);
SL.W  =  SL.W(:);
SL.P  =  SL.P(:);

CTX.SL  =  SL;
CTX.MP  =  MP;
CTX.BC  =  BC;


end

