% Smoothfield    EREBUS subroutine to regularise fields by smoothing
%
% [field]  =  SmoothField(field,kappa,nsteps,FE,type)
%
%   Uses an explicit finite-difference lagrangian operation to smooth input field. 
%
%   field   : input field (can have any type, multiple components)
%   kappa   : input relative diffusivity for smoothing (0 <= kappa <= 1)
%   nsteps  : input number of smoothing steps to be applied
%   FE      : input structure with finite-element mesh information
%   type    : input field type ('El', 'Q1', 'Q2', 'IP')
%
%   field   : output regularised field
%
%   created   20140806  Tobias Keller
%   modified  20170427  Tobias Keller
%   modified  20200515  Tobias Keller

function    [field]  =  SmoothField(field,kappa,nsteps,FE,type)

ncomp = size(field,2);

switch type
    
    case 'El'
        
        map  =  FE.MapEl;
        nx   =  FE.nxEl;
        nz   =  FE.nzEl;
        xc   =  FE.CoordEl(:,1);
        zc   =  FE.CoordEl(:,2);

    case 'Q1'
        
        map  =  FE.MapQ1;
        nx   =  FE.nxQ1;
        nz   =  FE.nzQ1;
        xc   =  FE.CoordQ1(:,1);
        zc   =  FE.CoordQ1(:,2);
        
    case 'Q2'
        
        map  =  FE.MapQ2;
        nx   =  FE.nxQ2;
        nz   =  FE.nzQ2;
        xc   =  FE.CoordQ2(:,1);
        zc   =  FE.CoordQ2(:,2);
        
    case 'IP'
        
        map  =  FE.MapIP;
        nx   =  FE.nxIP;
        nz   =  FE.nzIP;
        xc   =  FE.CoordIP(:,1);
        zc   =  FE.CoordIP(:,2);        
end

for n = 1:nsteps
    
    for k = 1:ncomp
        
        % extract k-th component
        fk                   =  field(:,k);
        
        % ghost field with zero-gradient boundaries
        ff                   =  zeros(nz+2,nx+2);
        ff(2:end-1,2:end-1)  =  fk(map);
        ff(1  ,2:end-1)      =  fk(map(1  ,:));
        ff(end,2:end-1)      =  fk(map(end,:));
        ff(2:end-1,  1)      =  fk(map(:,  1));
        ff(2:end-1,end)      =  fk(map(:,end));
        
        % extract coordinates with ghosted boundary extensions
        xx                   =  zeros(nz+2,nx+2);
        xx(2:end-1,2:end-1)  =  xc(map);
        xx(1  ,2:end-1)      =  2.*xc(map(1  ,:)) - xc(map(2    ,:));
        xx(end,2:end-1)      =  2.*xc(map(end,:)) - xc(map(end-1,:));
        xx(2:end-1,  1)      =  2.*xc(map(:,  1)) - xc(map(:,    2));
        xx(2:end-1,end)      =  2.*xc(map(:,end)) - xc(map(:,end-1));
 
        zz                   =  zeros(nz+2,nx+2);
        zz(2:end-1,2:end-1)  =  zc(map);
        zz(1  ,2:end-1)      =  2.*zc(map(1  ,:)) - zc(map(2    ,:));
        zz(end,2:end-1)      =  2.*zc(map(end,:)) - zc(map(end-1,:));
        zz(2:end-1,  1)      =  2.*zc(map(:,  1)) - zc(map(:,    2));
        zz(2:end-1,end)      =  2.*zc(map(:,end)) - zc(map(:,end-1));
        
        h = mean([mean(mean(diff(xx(2:end-1,:),1,2))),mean(mean(diff(zz(:,2:end-1),1,1)))]);
        
        qz          =  -diff(ff(:,2:end-1),1,1)./diff(zz(:,2:end-1),1,1);
        qx          =  -diff(ff(2:end-1,:),1,2)./diff(xx(2:end-1,:),1,2);
        lagr        =  -diff(qz,1,1)./diff((zz(1:end-1,2:end-1)+zz(2:end,2:end-1))./2,1,1) ...
                       -diff(qx,1,2)./diff((xx(2:end-1,1:end-1)+xx(2:end-1,2:end))./2,1,2);
        fk(map)     =  fk(map) + kappa/2*(h/2)^2 .* lagr;
        
        field(:,k)  =  fk;
    end
    
end

end





