% PIPQ2    EREBUS subroutine to project field from integration points to Q2-nodes
%
% [Q2field]  =  PIPQ2(IPfield,FE)
%
%   Function returns Q2field projected from input IPfield on integration points 
%   space to Q2-nodes space.
%
%   created   20161115  Tobias Keller
%   modified  20170427  Tobias Keller
%   modified  20200515  Tobias Keller


function  [Q2field]  =  PIPQ2(IPfield,FE)

n        =  length(IPfield(1,:));
Q2field  =  zeros(FE.NQ2,n);
N        =  FE.NiQ2;

for k=1:n
    
    TMP           =  IPfield(:,k);
    TMP           =  (TMP(FE.El2IP)*N') ./ (ones(size(FE.El2IP))*N');
    Q2field(:,k)  =  accumarray(FE.El2Q2(:),TMP(:)) ./ accumarray(FE.El2Q2(:),ones(size(TMP(:))));
    
end

end