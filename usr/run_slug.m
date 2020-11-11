clear variables;

% read in default parameter choices and runtime options
par_EREBUS

% set model run identifier
CTX.IO.RunID           = 'slug4s';
CTX.IO.TryContinue     =  1;
CTX.IO.ContFrame       =  -1;

CTX.IO.LivePlot        = 'ON';
CTX.IO.nwrite          =  20;

CTX.FE.nx              =  50;
CTX.FE.nz              =  250;
CTX.FE.ElType          =  'Q2Q1';        

CTX.FE.W               =  20;    
CTX.FE.D               =  100;   
CTX.FE.LagrMesh        =  'SRF';

CTX.INIT.PertSmooth    =  round(15/(CTX.FE.D/CTX.FE.nz)^2); 

CTX.PHYS.grav          =  9.81;
CTX.PHYS.RConst        =  8.314;
CTX.PHYS.TauPhi        =  1*60;
CTX.PHYS.TauTmp        =  4*3600;
CTX.PHYS.delta         =  1;

CTX.PHYS.SlugNo        =  1; 
CTX.PROP.Phi0          =  0.10.*sqrt(CTX.PHYS.SlugNo);
CTX.PHYS.WIn           =  0.02.*sqrt(CTX.PHYS.SlugNo);
CTX.PHYS.SlugLength    =  CTX.FE.D./(CTX.PHYS.SlugNo).^2;                  
CTX.PHYS.TauIn         =  CTX.PHYS.SlugLength/CTX.PHYS.WIn;

CTX.INIT.MatHeight     =  0;
CTX.INIT.MatWidth      =  5;
CTX.INIT.MatXLoc       =  CTX.FE.W/2;
CTX.INIT.MatZLoc       =  CTX.INIT.MatHeight;
CTX.PROP.CntAur        =  10;                           
CTX.PROP.Compr         =  1.0;

CTX.INIT.AddBlock      =  CTX.PHYS.SlugNo+1;                    
CTX.INIT.BlockWidth    =  CTX.INIT.MatWidth/2*0.7.*ones(1,CTX.PHYS.SlugNo+1);
CTX.INIT.BlockHeight   =  CTX.PHYS.SlugLength.*ones(1,CTX.PHYS.SlugNo+1);
CTX.INIT.BlockZLoc     =  CTX.FE.D-((0:CTX.PHYS.SlugNo:CTX.PHYS.SlugNo^2))*CTX.PHYS.SlugLength;
CTX.INIT.BlockXLoc     =  CTX.INIT.MatXLoc.*ones(1,CTX.PHYS.SlugNo+1);
CTX.INIT.BlockMat      =  3.*ones(1,CTX.PHYS.SlugNo+1);

CTX.INIT.AddSphere     =  0;                    
CTX.INIT.SphereRadiusX =  CTX.INIT.MatWidth/2*0.75/2;
CTX.INIT.SphereRadiusZ =  CTX.INIT.BlockWidth/2;
CTX.INIT.SphereZLoc    =  CTX.FE.D;
CTX.INIT.SphereXLoc    =  CTX.INIT.MatXLoc;
CTX.INIT.SphereMat     =  3;

CTX.RHEO.Plasticity      =  'ON';                % switch for plastic failure mode
CTX.RHEO.Elasticity      =  'OFF';                % switch for elastic mode

CTX.RHEO.Strainr0        =  1.e-2;               % reference strain rate [1/s]
CTX.RHEO.PowerLawExp     =  2;                   % non-Newtonian powerlaw exponent

% run simulation code
run('../src/EREBUS.m')
