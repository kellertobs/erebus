
function  Elfield  =  PIPEl(IPfield,FE)

n        =  length(IPfield(1,:));
Elfield  =  zeros(FE.NEl,n);
N        =  FE.NiP0;

for k=1:n
    
    TMP           =  IPfield(:,k);
    Elfield(:,k)  =  (TMP(FE.El2IP)*N') ./ (ones(size(FE.El2IP))*N');
    
end

end