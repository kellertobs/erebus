% UpdateMaterialProps    EREBUS subroutine to update material point properties
%
% [CTX] = UpdateMaterialProps(CTX)
%
%   Function updates material properties and auxiliary fields on integration
%   points according to latest solution guess and parameters and options 
%   provided in the CTX structure.
%
%   created   20140730  Tobias Keller
%   modified  20170427  Tobias Keller
%   modified  20200227  Tobias Keller
%   modified  20200515  Tobias Keller


function [CTX]  =  UpdateMaterialPoints(CTX)

MP  =  CTX.MP;


%*****  Update density  ***************************************************

MP  =  Density(MP,CTX);


%*****  Update stress / strain rates  *************************************

MP  =  StressStrainr(MP,CTX);


%*****  Update rheology  **************************************************

MP  =  Rheology(MP,CTX);


CTX.MP  =  MP;

end







