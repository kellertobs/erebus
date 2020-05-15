% Determine the allowed time step for advection, diffusion and free surface
% [Time]  =  CalcTimestep(Time,VS,VF,Temp,Dens,Dim,El,No,MAP,BC)

function  [CTX]   =  CalcTimestep(CTX)

str   =  'dt0';
dt0   =  2*CTX.TIME.stepo;
h     =  CTX.FE.hzQ2;
if CTX.FE.NU == CTX.FE.NQ1
    mapv  =  CTX.FE.MapQ1;
else
    mapv  =  CTX.FE.MapQ2;
end
 
th     =  CTX.SL.theta_dt;

U      =  th.*CTX.SL.U + (1-th).*CTX.SLo.U;
W      =  th.*CTX.SL.W + (1-th).*CTX.SLo.W;


%***  solid velocity advection  *******************************************

advs_dt    =  CTX.SL.CFL .* min( min(h./2./abs(U(:)+1e-32)) , min(h./2./abs(W(:)+1e-32)) );

if advs_dt<=dt0; str = 'advs'; end
dt         =  min(dt0,advs_dt);


%***  surface deformation  ************************************************
    
if strcmp(CTX.BC.FreeSurface,'ON')
    velsurf  =  abs(W(mapv(1,:)));
    velsurf  =  velsurf-mean(velsurf);
    srf_dt   =  CTX.SL.CFL * min(h./10./abs(velsurf(:)));
    
    if srf_dt<=dt; str = 'srf'; end
    dt  =  min(dt,srf_dt);
end


%***  limit timestep  *****************************************************

CTX.TIME.step  =  min(dt0,dt);
CTX.TIME.limit =  str;

end


