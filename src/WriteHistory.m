% WriteHistory    EREBUS subroutine to update the history of key metrics
%
% [] = WriteHistory(CTX)
%
%   Function calculates and writes out to file a series of key output metrics
%   to track lava lake evolution. The metrics include heat and outgassing
%   fluxes into the base and from the top, and min, mean and max values of
%   solution variables and auxiliary fields on the lake surface and within
%   the interior.
%
%   created   20140730  Tobias Keller
%   modified  20170427  Tobias Keller
%   modified  20200227  Tobias Keller
%   modified  20200515  Tobias Keller


function  []  =  WriteHistory(CTX)

FE     =  CTX.FE;
SL     =  CTX.SL;
MP     =  CTX.MP;
PROP   =  CTX.PROP;

n      =  CTX.TIME.istep;

FileName  =  [CTX.IO.DataDir '/' CTX.IO.RunID '/' CTX.IO.RunID '_hist.mat'];

if exist(FileName,'file')
    load(FileName,'H');
    if n == 1
        clear H;
    end
end

%*****  record time  ******************************************************

H.time(n) =  CTX.TIME.total;
H.step(n) =  CTX.TIME.step;

vol      =  PElQ2(FE.ElVol./4,FE);
map      =  FE.MapQ2;
mapIP    =  FE.MapIP;
RhoM     =  PIPQ2(MP.Rho,FE);
RhoPhi   =  PIPQ2(MP.RhoPhi,FE);
CL       =  CTX.PROP.c + (1-SL.Phi).*PROP.L./(PROP.Tliq-PROP.Tsol);  % adjusted heat capacity
CL(SL.T < PROP.Tsol)  =  CTX.PROP.c;
Tc       =  CTX.RHEO.Chic .* (PROP.Tsol - PROP.Tliq) + PROP.Tliq;
lake     =  SL.T > Tc;
lakeIP   =  logical(PQ2IP(lake,FE));
bot      =  map(end,:);
top      =  map(1  ,:);
T        =  SL.T;
Tbot     =  CTX.BC.botTmp;
T0       =  CTX.INIT.TempExt;
Phi      =  SL.Phi;
Phibot   =  CTX.BC.botPhi;
hbot     =  diff(FE.CoordQ2(bot,1));
htop     =  diff(FE.CoordQ2(top,1));
Abot     =  sum(hbot); % area of bot surface to normalize fluxes
Atop     =  sum(htop); % area of bot surface to normalize fluxes
Surf     =  repmat(FE.SurfQ2.',FE.nzQ2,1);
W        =  -SL.W;
Z        =  FE.CoordQ2(:,2)-Surf(:);
vol      =  (vol(map(1:end-1,1:end-1))+vol(map(2:end,1:end-1))+vol(map(1:end-1,2:end))+vol(map(2:end,2:end)))./4;

%***  record heat fluxes

% heat flux into bottom of domain [W/m2], positive for heat addition
H.bot.HeatIn (n)  =  sum(sum((RhoM(map(1:end-1,1:end-1)).*CL(map(1:end-1,1:end-1)).*SL.HeatIn (map(1:end-1,1:end-1)) ...
                          +   RhoM(map(2:end  ,1:end-1)).*CL(map(2:end  ,1:end-1)).*SL.HeatIn (map(2:end  ,1:end-1)) ...
                          +   RhoM(map(1:end-1,2:end  )).*CL(map(1:end-1,2:end  )).*SL.HeatIn (map(1:end-1,2:end  )) ...
                          +   RhoM(map(2:end  ,2:end  )).*CL(map(2:end  ,2:end  )).*SL.HeatIn (map(2:end  ,2:end  )))./4 .* vol)) / Abot;

% heat flux out of top of domain [W/m2], positive for heat loss
H.top.HeatOut(n)  = -sum(sum((RhoM(map(1:end-1,1:end-1)).*CL(map(1:end-1,1:end-1)).*SL.HeatOut(map(1:end-1,1:end-1)) ...
                          +   RhoM(map(2:end  ,1:end-1)).*CL(map(2:end  ,1:end-1)).*SL.HeatOut(map(2:end  ,1:end-1)) ...
                          +   RhoM(map(1:end-1,2:end  )).*CL(map(1:end-1,2:end  )).*SL.HeatOut(map(1:end-1,2:end  )) ...
                          +   RhoM(map(2:end  ,2:end  )).*CL(map(2:end  ,2:end  )).*SL.HeatOut(map(2:end  ,2:end  )))./4 .* vol)) / Atop;
                     
                     
%***  record gas fluxes

% gas mass flux into bottom of domain [W/m2], positive for gas addition
H.bot.GasIn (n)  =  sum(sum((RhoPhi(map(1:end-1,1:end-1)).*SL.GasIn(map(1:end-1,1:end-1)) ...
                         +   RhoPhi(map(2:end  ,1:end-1)).*SL.GasIn(map(2:end  ,1:end-1)) ...
                         +   RhoPhi(map(1:end-1,2:end  )).*SL.GasIn(map(1:end-1,2:end  )) ...
                         +   RhoPhi(map(2:end  ,2:end  )).*SL.GasIn(map(2:end  ,2:end  )))./4 .* vol)) / Abot;
                   
% gas mass flux out of top of domain [W/m2], positive for gas loss
H.top.GasOut(n)  = -sum(sum((RhoPhi(map(1:end-1,1:end-1)).*SL.GasOut(map(1:end-1,1:end-1)) ...
                         +   RhoPhi(map(2:end  ,1:end-1)).*SL.GasOut(map(2:end  ,1:end-1)) ...
                         +   RhoPhi(map(1:end-1,2:end  )).*SL.GasOut(map(1:end-1,2:end  )) ...
                         +   RhoPhi(map(2:end  ,2:end  )).*SL.GasOut(map(2:end  ,2:end  )))./4 .* vol)) / Atop;


%***  record flow speed diagnostics
Fsp  =  sqrt(SL.U.^2 + SL.W.^2);
H.lake.Fsp(n,1)  =  min (Fsp(lake));
H.lake.Fsp(n,2)  =  mean(Fsp(lake));
H.lake.Fsp(n,3)  =  max (Fsp(lake));

H.top .Fsp(n,1)  =  min (Fsp(map(1,:)));
H.top .Fsp(n,2)  =  mean(Fsp(map(1,:)));
H.top .Fsp(n,3)  =  max (Fsp(map(1,:)));

%***  record T diagnostics
H.lake.T(n,1)  =  min (SL.T(lake));
H.lake.T(n,2)  =  mean(SL.T(lake));
H.lake.T(n,3)  =  max (SL.T(lake));

H.top .T(n,1)  =  min (SL.T(map(1,:)));
H.top .T(n,2)  =  mean(SL.T(map(1,:)));
H.top .T(n,3)  =  max (SL.T(map(1,:)));

%***  record vesicularity diagnostics
H.lake.phi(n,1)  =  min (SL.Phi(lake));
H.lake.phi(n,2)  =  mean(SL.Phi(lake));
H.lake.phi(n,3)  =  max (SL.Phi(lake));

H.top .phi(n,1)  =  min (SL.Phi(map(1,:)));
H.top .phi(n,2)  =  mean(SL.Phi(map(1,:)));
H.top .phi(n,3)  =  max (SL.Phi(map(1,:)));

%***  record crystallinity diagnostics
H.lake.chi(n,1)  =  min (SL.Chi(lake));
H.lake.chi(n,2)  =  mean(SL.Chi(lake));
H.lake.chi(n,3)  =  max (SL.Chi(lake));

H.top .chi(n,1)  =  min (SL.Chi(map(1,:)));
H.top .chi(n,2)  =  mean(SL.Chi(map(1,:)));
H.top .chi(n,3)  =  max (SL.Chi(map(1,:)));

%***  record density diagnostics
H.lake.Rho(n,1)  =  min (MP.Rho(lakeIP));
H.lake.Rho(n,2)  =  mean(MP.Rho(lakeIP));
H.lake.Rho(n,3)  =  max (MP.Rho(lakeIP));

H.top .Rho(n,1)  =  min (MP.Rho(mapIP(1,:)));
H.top .Rho(n,2)  =  mean(MP.Rho(mapIP(1,:)));
H.top .Rho(n,3)  =  max (MP.Rho(mapIP(1,:)));

%***  record viscosity diagnostics
H.lake.Eta(n,1)  =  min    (MP.Eta(lakeIP));
H.lake.Eta(n,2)  =  geomean(MP.Eta(lakeIP));
H.lake.Eta(n,3)  =  max    (MP.Eta(lakeIP));

H.top .Eta(n,1)  =  min    (MP.Eta(mapIP(1,:)));
H.top .Eta(n,2)  =  geomean(MP.Eta(mapIP(1,:)));
H.top .Eta(n,3)  =  max    (MP.Eta(mapIP(1,:)));

%***  record stress diagnostics
H.lake.TII(n,1)  =  min    (MP.TII(lakeIP,1));
H.lake.TII(n,2)  =  geomean(MP.TII(lakeIP,1));
H.lake.TII(n,3)  =  max    (MP.TII(lakeIP,1));

H.top .TII(n,1)  =  min    (MP.TII(mapIP(1,:),1));
H.top .TII(n,2)  =  geomean(MP.TII(mapIP(1,:),1));
H.top .TII(n,3)  =  max    (MP.TII(mapIP(1,:),1));

%***  record strain rate diagnostics
H.lake.EII(n,1)  =  min    (MP.EII(lakeIP,1));
H.lake.EII(n,2)  =  geomean(MP.EII(lakeIP,1));
H.lake.EII(n,3)  =  max    (MP.EII(lakeIP,1));

H.top .EII(n,1)  =  min    (MP.EII(mapIP(1,:),1));
H.top .EII(n,2)  =  geomean(MP.EII(mapIP(1,:),1));
H.top .EII(n,3)  =  max    (MP.EII(mapIP(1,:),1));

%***  record location of maximum spreading
[~,ind]             =  max(MP.EII(mapIP(1,:),1));
H.top.MaxSprLoc(n)  =  FE.CoordIP(mapIP(1,ind),1);

%***  calculate "plateness" after Tackley et al (2000)
EIItotal  =  sum (MP.EII(mapIP(1,:),1));
Edotsort  =  sort(MP.EII(mapIP(1,:),1),'descend');

i   = 0;
EIIsum = 0;
while EIIsum<0.8*EIItotal
    i = i + 1;
    EIIsum = EIIsum + Edotsort(i);
end
f80         =  i/FE.nxIP;     % fraction of surface area in which highest 80% of strainr occurs
H.top.P(n)  =  1 - f80/0.6;   % P=1 for plate-like, P=0 for constant viscosity convection


%*****  save history structure  *******************************************

save(FileName,'H');

end