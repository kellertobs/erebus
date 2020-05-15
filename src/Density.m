
% Density    EDIFICE: Update material density
%
% [CTX] = Density(MP,CTX)
%
%   Function updates T-dependent density according to latest solution guess. 
%
%   created   20140730  Tobias Keller
%   modified  20170427  Tobias Keller
%   modified  20200227   Tobias Keller


function   [MP]  =  Density(MP,CTX)


%*****  get densities from material types  ********************************

Rho0       =  CTX.PROP.Rho;
MP.RhoPhi  =  CTX.PROP.MPhi.*PQ1IP(CTX.SL.Pt,CTX.FE)./CTX.PHYS.RConst./(PQ2IP(CTX.SL.T,CTX.FE)+273.15);
MP.Rho     =  (1-MP.Phi-MP.Chi).*Rho0 + MP.Phi.*MP.RhoPhi + MP.Chi.*CTX.PROP.RhoChi;
MP.Rho     =  MP.Rho .* (1 - CTX.PROP.alpha.*(PQ2IP(CTX.SL.T,CTX.FE)-CTX.PROP.Tsol));

end





