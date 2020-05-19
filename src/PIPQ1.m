% PIPQ1    EREBUS subroutine to project field from integration points to Q1-nodes
%
% [Q1field]  =  PIPQ1(IPfield,FE)
%
%   Function returns Q1field projected from input IPfield on integration points 
%   space to Q1-nodes space.
%
%   created   20161115  Tobias Keller
%   modified  20170427  Tobias Keller
%   modified  20200515  Tobias Keller


function  [Q1field]  =  PIPQ1(IPfield,FE)

n        =  length(IPfield(1,:));
Q1field  =  zeros(FE.NQ1,n);
N        =  FE.NiQ1;

for k=1:n
    
    TMP           =  IPfield(:,k);
    TMP           =  (TMP(FE.El2IP)*N') ./ (ones(size(FE.El2IP))*N');
    Q1field(:,k)  =  accumarray(FE.El2Q1(:),TMP(:)) ./ accumarray(FE.El2Q1(:),ones(size(TMP(:))));
    
end

end