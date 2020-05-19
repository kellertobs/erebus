% AdvDiffSolver    EREBUS subroutine to solve advection-diffusion problem
%
% [field] = AdvDiffSolver(field,oldfield,vel,oldvel,src,srco,kappa,dt,type,bval,CTX)
%
%   Function solves time-dependent advection-diffusion problem by explicit
%   time stepping and staggered finite-difference spatial discretisation. 
%   Various advection schemes are available:
%   - UPW2: second order upwind scheme
%   - UPW3: third order upwind scheme
%   - FRM : upwind-biased flux-conservative Fromm scheme
%   - FTV : central finite-volume flux-conservative scheme
%
%   The routine assumes irregular grid spacing but should only be used on
%   minimally distorted meshes close to rectangular geometry.
%
%   field   : input solution field (can have any type and multiple components)
%   fieldo  : input solution field of previous time step
%   vel     : input velocity field for advection
%   velo    : input velocity field of previous time step
%   src     : input source field of same dimension as solution
%   srco    : input source field of previous time step
%   kappa   : input diffusivity field of same dimensions as solution
%   dt      : input discrete time step
%   bval    : input cell array with boundary conditions
%   CTX     : input context structure
%
%   created   20140730   Tobias Keller
%   modified  20170427   Tobias Keller
%   modified  20200227   Tobias Keller
%   modified  20200515   Tobias Keller


function [field] = AdvDiffSolver(field,fieldo,vel,velo,src,srco,kappa,dt,type,bval,CTX)

FE     =  CTX.FE;
FEo    =  CTX.FEo;
ncomp  =  size(field,2);

switch type
    
    case 'El'
        
        map  =  FE.MapEl;
        nx   =  FE.nxEl;
        nz   =  FE.nzEl;
        xc   =  FE.CoordEl(:,1);
        zc   =  FE.CoordEl(:,2);
        xco  =  FEo.CoordEl(:,1);
        zco  =  FEo.CoordEl(:,2);
        if FE.NU == FE.NQ2
            vel     =  PQ2El(vel,FE);
            velo  =  PQ2El(velo,FE);
        else
            vel     =  PQ1El(vel,FE);
            velo  =  PQ1El(velo,FE);
        end

    case 'Q1'
        
        h    =  FE.hzQ1;
        map  =  FE.MapQ1;
        nx   =  FE.nxQ1;
        nz   =  FE.nzQ1;
        xc   =  FE.CoordQ1(:,1);
        zc   =  FE.CoordQ1(:,2);
        xco  =  FEo.CoordQ1(:,1);
        zco  =  FEo.CoordQ1(:,2);
        if FE.NU == FE.NQ2
            vel     =  PQ2Q1(vel,FE);
            velo  =  PQ2Q1(velo,FE);
        end
        
    case 'Q2'
        
        h    =  FE.hzQ2;
        map  =  FE.MapQ2;
        nx   =  FE.nxQ2;
        nz   =  FE.nzQ2;
        xc   =  FE.CoordQ2(:,1);
        zc   =  FE.CoordQ2(:,2);
        xco  =  FEo.CoordQ2(:,1);
        zco  =  FEo.CoordQ2(:,2);
        if FE.NU == FE.NQ1
            vel     =  PQ1Q2(vel,FE);
            velo  =  PQ1Q2(velo,FE);
        end
        
end

ic   =  3:nx+2;
im   =  2:nx+1;
ip   =  4:nx+3;
imm  =  1:nx+0;
ipp  =  5:nx+4;

jc   =  3:nz+2;
jm   =  2:nz+1;
jp   =  4:nz+3;
jmm  =  1:nz+0;
jpp  =  5:nz+4;

ZF   =  {'zf','zf','zf','zf'};  % zero flux boundaries
CT   =  {'ct','ct','ct','ct'};  % continuation boundaries

for k = 1:ncomp
    
    fieldk   =     field(:,k);
    fieldko  =  fieldo(:,k);
    
    xx       =  ghost_field(xc ,nx,nz,CT,map);
    xxo      =  ghost_field(xco,nx,nz,CT,map);

    zz       =  ghost_field(zc ,nx,nz,CT,map);
    zzo      =  ghost_field(zco,nx,nz,CT,map);

    a        =  ghost_field(fieldk ,nx,nz,bval,map);
    ao       =  ghost_field(fieldko,nx,nz,bval,map);

    qx       =  -diff(a(jc,2:end-1),1,2)./diff(xx(jc,2:end-1),1,2);
    qz       =  -diff(a(2:end-1,ic),1,1)./diff(zz(2:end-1,ic),1,1);
    
    lapl     =  -diff(qx,1,2)./diff((xx(jc,2:end-2)+xx(jc,3:end-1))./2,1,2) ...
                -diff(qz,1,1)./diff((zz(2:end-2,ic)+zz(3:end-1,ic))./2,1,1);
            
    qxo      =  -diff(ao(jc,2:end-1),1,2)./diff(xxo(jc,2:end-1),1,2);
    qzo      =  -diff(ao(2:end-1,ic),1,1)./diff(zzo(2:end-1,ic),1,1);
    
    laplo    =  -diff(qxo,1,2)./diff((xxo(jc,2:end-2)+xxo(jc,3:end-1))./2,1,2) ...
                -diff(qzo,1,1)./diff((zzo(2:end-2,ic)+zzo(3:end-1,ic))./2,1,1);
            
    u        =  ghost_field(vel(:,1),nx,nz,ZF,map);
    w        =  ghost_field(vel(:,2),nx,nz,ZF,map);
    
    uo       =  ghost_field(velo(:,1),nx,nz,ZF,map);
    wo       =  ghost_field(velo(:,2),nx,nz,ZF,map);
        
    if strcmp(CTX.SL.Advection(1:3),'UPW')
        
        um       =  min(u(jc,ic),0);
        up       =  max(u(jc,ic),0);
        wm       =  min(w(jc,ic),0);
        wp       =  max(w(jc,ic),0);
        
        umo      =  min(uo(jc,ic),0);
        upo      =  max(uo(jc,ic),0);
        wmo      =  min(wo(jc,ic),0);
        wpo      =  max(wo(jc,ic),0);
        
        if strcmp(CTX.SL.Advection(4),'2')
            
            da_dzm =   3/2 * diff(a(2:end-2,ic),1,1)./diff(zz(2:end-2,ic),1,1) ...
                     - 1/2 * diff(a(1:end-3,ic),1,1)./diff(zz(1:end-3,ic),1,1);
            da_dzp =   3/2 * diff(a(3:end-1,ic),1,1)./diff(zz(3:end-1,ic),1,1) ...
                     - 1/2 * diff(a(4:end-0,ic),1,1)./diff(zz(4:end-0,ic),1,1);
            da_dxm =   3/2 * diff(a(jc,2:end-2),1,2)./diff(xx(jc,2:end-2),1,2) ...
                     - 1/2 * diff(a(jc,1:end-3),1,2)./diff(xx(jc,1:end-3),1,2);
            da_dxp =   3/2 * diff(a(jc,3:end-1),1,2)./diff(xx(jc,3:end-1),1,2) ...
                     - 1/2 * diff(a(jc,4:end-0),1,2)./diff(xx(jc,4:end-0),1,2);
            
            adv     =  um .* da_dxp + up .* da_dxm + wm .* da_dzp + wp .* da_dzm;
            
            da_dzm =   3/2 * diff(ao(2:end-2,ic),1,1)./diff(zzo(2:end-2,ic),1,1) ...
                     - 1/2 * diff(ao(1:end-3,ic),1,1)./diff(zzo(1:end-3,ic),1,1);
            da_dzp =   3/2 * diff(ao(3:end-1,ic),1,1)./diff(zzo(3:end-1,ic),1,1) ...
                     - 1/2 * diff(ao(4:end-0,ic),1,1)./diff(zzo(4:end-0,ic),1,1);
            da_dxm =   3/2 * diff(ao(jc,2:end-2),1,2)./diff(xxo(jc,2:end-2),1,2) ...
                     - 1/2 * diff(ao(jc,1:end-3),1,2)./diff(xxo(jc,1:end-3),1,2);
            da_dxp =   3/2 * diff(ao(jc,3:end-1),1,2)./diff(xxo(jc,3:end-1),1,2) ...
                     - 1/2 * diff(ao(jc,4:end-0),1,2)./diff(xxo(jc,4:end-0),1,2);
            
            advo    =  umo .* da_dxp + upo .* da_dxm + wmo .* da_dzp + wpo .* da_dzm;
            
        elseif strcmp(CTX.SL.Advection(4),'3')
            
            da_dxm  =    2/6 * diff(a(jc,3:end-1),1,2)./diff(xx(jc,3:end-1),1,2) ...
                       + 5/6 * diff(a(jc,2:end-2),1,2)./diff(xx(jc,3:end-1),1,2) ...
                       - 1/6 * diff(a(jc,1:end-3),1,2)./diff(xx(jc,1:end-3),1,2);
            da_dxp  =    2/6 * diff(a(jc,2:end-2),1,2)./diff(xx(jc,2:end-2),1,2) ...
                       + 5/6 * diff(a(jc,3:end-1),1,2)./diff(xx(jc,3:end-1),1,2) ...
                       - 1/6 * diff(a(jc,4:end-0),1,2)./diff(xx(jc,4:end-0),1,2);
            da_dzm  =    2/6 * diff(a(3:end-1,ic),1,1)./diff(zz(3:end-1,ic),1,1) ...
                       + 5/6 * diff(a(2:end-2,ic),1,1)./diff(zz(3:end-1,ic),1,1) ...
                       - 1/6 * diff(a(1:end-3,ic),1,1)./diff(zz(1:end-3,ic),1,1);
            da_dzp  =    2/6 * diff(a(2:end-2,ic),1,1)./diff(zz(2:end-2,ic),1,1) ...
                       + 5/6 * diff(a(3:end-1,ic),1,1)./diff(zz(3:end-1,ic),1,1) ...
                       - 1/6 * diff(a(4:end-0,ic),1,1)./diff(zz(4:end-0,ic),1,1);
            
            adv     =  um .* da_dxp + up .* da_dxm + wm .* da_dzp + wp .* da_dzm;
            
            da_dxm  =    2/6 * diff(ao(jc,3:end-1),1,2)./diff(xxo(jc,3:end-1),1,2) ...
                       + 5/6 * diff(ao(jc,2:end-2),1,2)./diff(xxo(jc,3:end-1),1,2) ...
                       - 1/6 * diff(ao(jc,1:end-3),1,2)./diff(xxo(jc,1:end-3),1,2);
            da_dxp  =    2/6 * diff(ao(jc,2:end-2),1,2)./diff(xxo(jc,2:end-2),1,2) ...
                       + 5/6 * diff(ao(jc,3:end-1),1,2)./diff(xxo(jc,3:end-1),1,2) ...
                       - 1/6 * diff(ao(jc,4:end-0),1,2)./diff(xxo(jc,4:end-0),1,2);
            da_dzm  =    2/6 * diff(ao(3:end-1,ic),1,1)./diff(zzo(3:end-1,ic),1,1) ...
                       + 5/6 * diff(ao(2:end-2,ic),1,1)./diff(zzo(3:end-1,ic),1,1) ...
                       - 1/6 * diff(ao(1:end-3,ic),1,1)./diff(zzo(1:end-3,ic),1,1);
            da_dzp  =    2/6 * diff(ao(2:end-2,ic),1,1)./diff(zzo(2:end-2,ic),1,1) ...
                       + 5/6 * diff(ao(3:end-1,ic),1,1)./diff(zzo(3:end-1,ic),1,1) ...
                       - 1/6 * diff(ao(4:end-0,ic),1,1)./diff(zzo(4:end-0,ic),1,1);
            
            advo    =  umo .* da_dxp + upo .* da_dxm + wmo .* da_dzp + wpo .* da_dzm;
            
        end
        
    elseif strcmp(CTX.SL.Advection(1:3),'FRM')
        
        wp    =  (w(jp,ic)+w(jc,ic))./2;  zp = (zz(jp,ic)+zz(jc,ic))./2;
        wm    =  (w(jm,ic)+w(jc,ic))./2;  zm = (zz(jm,ic)+zz(jc,ic))./2;
        up    =  (u(jc,ip)+u(jc,ic))./2;  xp = (xx(jc,ip)+xx(jc,ic))./2;
        um    =  (u(jc,im)+u(jc,ic))./2;  xm = (xx(jc,im)+xx(jc,ic))./2;
        divv  =  (wp-wm)./(zp-zm) + (up-um)./(xp-xm);
        
        acc   =  a(jc,ic);  xcc = xx(jc,ic);  zcc = zz(jc,ic);
        ajp   =  a(jp,ic); ajpp = a(jpp,ic);  zjp = zz(jp,ic);  zjpp = zz(jpp,ic);
        ajm   =  a(jm,ic); ajmm = a(jmm,ic);  zjm = zz(jm,ic);  zjmm = zz(jmm,ic);
        aip   =  a(jc,ip); aipp = a(jc,ipp);  xip = xx(jc,ip);  xipp = xx(jc,ipp);
        aim   =  a(jc,im); aimm = a(jc,imm);  xim = xx(jc,im);  ximm = xx(jc,imm);

        adv   =     up .*(-(aipp-aip)./(xipp-xip)./8 + (aip + acc)./(xip-xcc)./2 + (acc-aim )./(xcc-xim )./8) ...
              - abs(up).*(-(aipp-aip)./(xipp-xip)./8 + (aip - acc)./(xip-xcc)./4 - (acc-aim )./(xcc-xim )./8) ...
              -     um .*(-(aip -acc)./(xip -xcc)./8 + (acc + aim)./(xcc-xim)./2 + (aim-aimm)./(xim-ximm)./8) ...
              + abs(um).*(-(aip -acc)./(xip -xcc)./8 + (acc - aim)./(xcc-xim)./4 - (aim-aimm)./(xim-ximm)./8) ...
              +     wp .*(-(ajpp-ajp)./(zjpp-zjp)./8 + (ajp + acc)./(zjp-zcc)./2 + (acc-ajm )./(zcc-zjm )./8) ...
              - abs(wp).*(-(ajpp-ajp)./(zjpp-zjp)./8 + (ajp - acc)./(zjp-zcc)./4 - (acc-ajm )./(zcc-zjm )./8) ...
              -     wm .*(-(ajp -acc)./(zjp -zcc)./8 + (acc + ajm)./(zcc-zjm)./2 + (ajm-ajmm)./(zjm-zjmm)./8) ...
              + abs(wm).*(-(ajp -acc)./(zjp -zcc)./8 + (acc - ajm)./(zcc-zjm)./4 - (ajm-ajmm)./(zjm-zjmm)./8);
        adv   =  adv - acc.*divv;
        
        wp    =  (wo(jp,ic)+wo(jc,ic))./2;  zp = (zzo(jp,ic)+zzo(jc,ic))./2;
        wm    =  (wo(jm,ic)+wo(jc,ic))./2;  zm = (zzo(jm,ic)+zzo(jc,ic))./2;
        up    =  (uo(jc,ip)+uo(jc,ic))./2;  xp = (xxo(jc,ip)+xxo(jc,ic))./2;
        um    =  (uo(jc,im)+uo(jc,ic))./2;  xm = (xxo(jc,im)+xxo(jc,ic))./2;
        divv  =  (wp-wm)./(zp-zm) + (up-um)./(xp-xm);
        
        acc   =  ao(jc,ic);  xcc = xxo(jc,ic);  zcc = zzo(jc,ic);
        ajp   =  ao(jp,ic); ajpp = ao(jpp,ic);  zjp = zzo(jp,ic);  zjpp = zzo(jpp,ic);
        ajm   =  ao(jm,ic); ajmm = ao(jmm,ic);  zjm = zzo(jm,ic);  zjmm = zzo(jmm,ic);
        aip   =  ao(jc,ip); aipp = ao(jc,ipp);  xip = xxo(jc,ip);  xipp = xxo(jc,ipp);
        aim   =  ao(jc,im); aimm = ao(jc,imm);  xim = xxo(jc,im);  ximm = xxo(jc,imm);
        
        advo  =     up .*(-(aipp-aip)./(xipp-xip)./8 + (aip + acc)./(xip-xcc)./2 + (acc-aim )./(xcc-xim )./8) ...
              - abs(up).*(-(aipp-aip)./(xipp-xip)./8 + (aip - acc)./(xip-xcc)./4 - (acc-aim )./(xcc-xim )./8) ...
              -     um .*(-(aip -acc)./(xip -xcc)./8 + (acc + aim)./(xcc-xim)./2 + (aim-aimm)./(xim-ximm)./8) ...
              + abs(um).*(-(aip -acc)./(xip -xcc)./8 + (acc - aim)./(xcc-xim)./4 - (aim-aimm)./(xim-ximm)./8) ...
              +     wp .*(-(ajpp-ajp)./(zjpp-zjp)./8 + (ajp + acc)./(zjp-zcc)./2 + (acc-ajm )./(zcc-zjm )./8) ...
              - abs(wp).*(-(ajpp-ajp)./(zjpp-zjp)./8 + (ajp - acc)./(zjp-zcc)./4 - (acc-ajm )./(zcc-zjm )./8) ...
              -     wm .*(-(ajp -acc)./(zjp -zcc)./8 + (acc + ajm)./(zcc-zjm)./2 + (ajm-ajmm)./(zjm-zjmm)./8) ...
              + abs(wm).*(-(ajp -acc)./(zjp -zcc)./8 + (acc - ajm)./(zcc-zjm)./4 - (ajm-ajmm)./(zjm-zjmm)./8);
        advo  =  advo - acc.*divv;
        
    elseif strcmp(CTX.SL.Advection(1:3),'FTV')
        
        u      =  ghost_field(   vel(:,1),nx,nz,ZF,map);
        w      =  ghost_field(   vel(:,2),nx,nz,ZF,map);
        
        uo     =  ghost_field(velo(:,1),nx,nz,ZF,map);
        wo     =  ghost_field(velo(:,2),nx,nz,ZF,map);
        
        wp     =  (w(jp,ic)+w(jc,ic))./2;  zp = (zz(jp,ic)+zz(jc,ic))./2;
        wm     =  (w(jm,ic)+w(jc,ic))./2;  zm = (zz(jm,ic)+zz(jc,ic))./2;
        up     =  (u(jc,ip)+u(jc,ic))./2;  xp = (xx(jc,ip)+xx(jc,ic))./2;
        um     =  (u(jc,im)+u(jc,ic))./2;  xm = (xx(jc,im)+xx(jc,ic))./2;
        divv   =  (wp-wm)./(zp-zm) + (up-um)./(xp-xm);
        
        acc  =  a(jc,ic);
        ajp  =  a(jp,ic);
        ajm  =  a(jm,ic);
        aip  =  a(jc,ip);
        aim  =  a(jc,im);
        
        adv  =  ((acc+aip)./2.*up - (acc+aim)./2.*um)./(xp-xm) ...
              + ((acc+ajp)./2.*wp - (acc+ajm)./2.*wm)./(zp-zm);
        adv  =  adv - acc.*divv;
             
        wp     =  (wo(jp,ic)+wo(jc,ic))./2;  zp = (zzo(jp,ic)+zzo(jc,ic))./2;
        wm     =  (wo(jm,ic)+wo(jc,ic))./2;  zm = (zzo(jm,ic)+zzo(jc,ic))./2;
        up     =  (uo(jc,ip)+uo(jc,ic))./2;  xp = (xxo(jc,ip)+xxo(jc,ic))./2;
        um     =  (uo(jc,im)+uo(jc,ic))./2;  xm = (xxo(jc,im)+xxo(jc,ic))./2;
        divv   =  (wp-wm)./(zp-zm) + (up-um)./(xp-xm);
        
        acc = ao(jc,ic);
        ajp = ao(jp,ic);
        ajm = ao(jm,ic);
        aip = ao(jc,ip);
        aim = ao(jc,im);
        
        advo =  ((acc+aip)./2.*up - (acc+aim)./2.*um)./(xp-xm) ...
              + ((acc+ajp)./2.*wp - (acc+ajm)./2.*wm)./(zp-zm);
        advo =  advo - acc.*divv;

    else
        adv  = 0;
        advo = 0;
    end

    fieldk(map) =  fieldko(map) - ((adv+advo)./2 - kappa(map).*(lapl+laplo)./2) .* dt;
    fieldk      =  fieldk       +  (src+srco)./2 .* dt;
    field(:,k)  =  fieldk;
    
end
   
end


function  [fgh]  =  ghost_field(f,nx,nz,bval,map)

ic  =  3:nx+2;
jc  =  3:nz+2;

fgh           =  zeros(nz+4,nx+4);
fgh(jc,ic)    =  f(map);

if strcmp(bval{1},'zf')
    fgh(2,ic)  =  fgh(3,ic);
    fgh(1,ic)  =  fgh(3,ic);
elseif strcmp(bval{1},'ct')
    fgh(2,ic)  =  2.*fgh(3,ic) - fgh(4,ic);
    fgh(1,ic)  =  2.*fgh(2,ic) - fgh(3,ic) ;   
else
    fgh(2,ic)  =  bval{1};
    fgh(1,ic)  =  bval{1};
end

if strcmp(bval{2},'zf')
    fgh(nz+3,ic)  =  fgh(nz+2,ic);
    fgh(nz+4,ic)  =  fgh(nz+2,ic);
elseif strcmp(bval{2},'ct')
    fgh(nz+3,ic)  =  2.*fgh(nz+2,ic) - fgh(nz+1,ic);
    fgh(nz+4,ic)  =  2.*fgh(nz+3,ic) - fgh(nz+2,ic);
else
    fgh(nz+3,ic)  =  bval{2};
    fgh(nz+4,ic)  =  bval{2};
end

if strcmp(bval{3},'zf')
    fgh(jc,2)  =  fgh(jc,3);
    fgh(jc,1)  =  fgh(jc,3);
elseif strcmp(bval{3},'ct')
    fgh(jc,2)  =  2.*fgh(jc,3) - fgh(jc,4);
    fgh(jc,1)  =  2.*fgh(jc,2) - fgh(jc,3);   
else
    fgh(jc,2)  =  bval{3};
    fgh(jc,1)  =  bval{3};
end

if strcmp(bval{4},'zf')
    fgh(jc,nx+3)  =  fgh(jc,nz+2);
    fgh(jc,nx+4)  =  fgh(jc,nz+2);
elseif strcmp(bval{4},'ct')
    fgh(jc,nx+3)  =  2.*fgh(jc,nz+2) - fgh(jc,nz+1);
    fgh(jc,nx+4)  =  2.*fgh(jc,nz+3) - fgh(jc,nz+2);
else
    fgh(jc,nx+3)  =  bval{4};
    fgh(jc,nx+4)  =  bval{4};
end

end





