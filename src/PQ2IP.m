% PQ2IP    EREBUS subroutine to project field from Q2-nodes to integration points
%
% [IPfield]  =  PQ2IP(Q2field,FE)
%
%   Function returns IPfield projected from input Q2field on Q2-nodes space
%   to integration points space.
%
%   created   20161115  Tobias Keller
%   modified  20170427  Tobias Keller
%   modified  20200515  Tobias Keller


function   [IPfield]  =  PQ2IP(Q2field,FE)

n        =  length(Q2field(1,:));
IPfield  =  zeros(FE.NIP,n);
N        =  FE.NiQ2;

for k=1:n
    
    TMP                     =  Q2field(:,k);
    TMP                     =  TMP(FE.El2Q2)*N;
    IPfield(FE.El2IP(:),k)  =  TMP(:);
    
end