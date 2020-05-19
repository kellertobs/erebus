% PQ1El    EREBUS subroutine to project field from Q1-nodes to elements
%
% [Elfield]  =  PQ1El(Q1field,FE)
%
%   Function returns Q1field projected from input Elfield on elements space
%   to Q1-nodes space.
%
%   created   20161115  Tobias Keller
%   modified  20170427  Tobias Keller
%   modified  20200515  Tobias Keller

function   [Elfield]  =  PQ1El(Q1field,FE)

n        =  length(Q1field(1,:));
Elfield  =  zeros(FE.NEl,n);
[N,~]    =  ShapeFuncts([0;0],'Q1');


for k=1:n
    
    TMP           =  Q1field(:,k);
    TMP           =  TMP(FE.El2Q1)*N;
    Elfield(:,k)  =  TMP(:);
    
end