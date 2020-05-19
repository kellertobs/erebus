% PQ1IP    EREBUS subroutine to project field from Q1-nodes to integration points
%
% [IPfield]  =  PQ1IP(Q1field,FE)
%
%   Function returns IPfield projected from input Q1field on Q1-nodes
%   space to integration points space.
%
%   created   20161115  Tobias Keller
%   modified  20170427  Tobias Keller
%   modified  20200515  Tobias Keller


function   [IPfield]  =  PQ1IP(Q1field,FE)

n        =  length(Q1field(1,:));
IPfield  =  zeros(FE.NIP,n);
N        =  FE.NiQ1;

for k=1:n
    
    TMP                     =  Q1field(:,k);
    TMP                     =  TMP(FE.El2Q1)*N;
    IPfield(FE.El2IP(:),k)  =  TMP(:);
    
end