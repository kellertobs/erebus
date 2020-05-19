% InitPert    EREBUS subroutine to initialise random perturbation field
%
% [pert]  =  InitPert(FE,n,sym)
%
%   Initialise random perturbation field to add controlled noise to different
%   fields to introduce heterogeneity or break symmetry. 
%
%   FE  : input finite-element mesh struct FE as set up by InitFE
%   n   : input n number of smoothing steps applied to redden white noise 
%   sym : input switch, 'ON' for perturbation field symmetric around vertical mid line of domain
%
%   created   20170427  Tobias Keller
%   modified  20200227  Tobias Keller
%   modified  20200515  Tobias Keller


function  [pert]  =  InitPert(FE,n,sym)

rng(5,'twister');

%*****  produce random perturbations  *************************************

% initialise white noise
a  =  rand(FE.NIP,1)-0.5;
a  =  SmoothField(a,1,n,FE,'IP');

if strcmp(sym,'ON')
    a(FE.MapIP(:,1:(FE.nxIP-1)/2)) = fliplr(a(FE.MapIP(:,(FE.nxIP-1)/2+2:end)));
end

pert  =  a./max(abs(a));

