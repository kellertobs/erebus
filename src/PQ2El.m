% PQ2El    EREBUS subroutine to project field from Q2-nodes to elements
%
% [Elfield]  =  PQ2El(Q2field,FE)
%
%   Function returns Q2field projected from input Elfield on elements
%   space to Q2-nodes space.
%
%   created   20161115  Tobias Keller
%   modified  20170427  Tobias Keller
%   modified  20200515  Tobias Keller


function   Elfield  =  PQ2El(Q2field,FE)

n        =  length(Q2field(1,:));
Elfield  =  zeros(FE.NEl,n);
[N,~]    =  ShapeFuncts([0;0],'Q2');


for k=1:n
    
    TMP           =  Q2field(:,k);
    TMP           =  TMP(FE.El2Q2)*N;
    Elfield(:,k)  =  TMP(:);
    
end