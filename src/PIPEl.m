% PIPEl    EREBUS subroutine to project field from integration points to elements
%
% [Elfield]  =  PIPEl(IPfield,FE)
%
%   Function returns Elfield projected from input IPfield on integration points 
%   space to elements space.
%
%   created   20161115  Tobias Keller
%   modified  20170427  Tobias Keller
%   modified  20200515  Tobias Keller


function  [Elfield]  =  PIPEl(IPfield,FE)

n        =  length(IPfield(1,:));
Elfield  =  zeros(FE.NEl,n);
N        =  FE.NiP0;

for k=1:n
    
    TMP           =  IPfield(:,k);
    Elfield(:,k)  =  (TMP(FE.El2IP)*N') ./ (ones(size(FE.El2IP))*N');
    
end

end