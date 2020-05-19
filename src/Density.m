% Density    EREBUS subroutine to update material density
%
% [CTX] = Density(MP,CTX)
%
%   Function updates density of three-phase magma according to latest solution guess. 
%
%   created   20140730   Tobias Keller
%   modified  20170427   Tobias Keller
%   modified  20200227   Tobias Keller
%   modified  20200515   Tobias Keller


function   [MP]  =  Density(MP,CTX)


%*****  get densities from material types  ********************************

% ideal gas density
MP.RhoPhi  =  CTX.PROP.MPhi.*PQ1IP(CTX.SL.Pt,CTX.FE)./CTX.PHYS.RConst./(PQ2IP(CTX.SL.T,CTX.FE)+273.15);

% T-dependent melt/crystal density
RhoMC      =  ((1-MP.Chi).*CTX.PROP.Rho + MP.Chi.*CTX.PROP.RhoChi) ...
           .* (1 - CTX.PROP.alpha.*(PQ2IP(CTX.SL.T,CTX.FE)-CTX.PROP.Tsol));

% melt/crystal/bubble mixture density
MP.Rho     =  (1-MP.Phi).*RhoMC + MP.Phi.*MP.RhoPhi;

end





