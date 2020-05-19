% PQ1Q2    EREBUS subroutine to project field from Q1-nodes to Q2-nodes
%
% [Q2field]  =  PQ1Q2(Q1field,FE)
%
%   Function returns Q2field projected from input Q1field on Q1-nodes
%   space to Q2-nodes space.
%
%   created   20161115  Tobias Keller
%   modified  20170427  Tobias Keller
%   modified  20200515  Tobias Keller


function   Q2field  =  PQ1Q2(Q1field,FE)

n        =  length(Q1field(1,:));
Q2field  =  zeros(FE.NQ2,n);
[N,~]    =  ShapeFuncts([-1,0,1,-1,0,1,-1,0,1;-1,-1,-1,0,0,0,1,1,1],'Q1');


for k=1:n
    
    TMP                     =  Q1field(:,k);
    TMP                     =  TMP(FE.El2Q1)*N;
    Q2field(FE.El2Q2(:),k)  =  TMP(:);
    
end