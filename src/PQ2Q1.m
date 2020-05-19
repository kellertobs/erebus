% PQ2Q1    EREBUS subroutine to project field from Q2-nodes to Q1-nodes
%
% [Q1field]  =  PQ2Q1(Q2field,FE)
%
%   Function returns Q1field projected from input Q2field on Q2-nodes space
%   to Q1-nodes space.
%
%   created   20161115  Tobias Keller
%   modified  20170427  Tobias Keller
%   modified  20200515  Tobias Keller


function   Q1field  =  PQ2Q1(Q2field,FE)

n        =  length(Q2field(1,:));
Q1field  =  zeros(FE.NQ1,n);
[N,~]    =  ShapeFuncts([-1,1,-1,1;-1,-1,1,1],'Q2');


for k=1:n
    
    TMP                     =  Q2field(:,k);
    TMP                     =  TMP(FE.El2Q2)*N;
    Q1field(FE.El2Q1(:),k)  =  TMP(:);
    
end