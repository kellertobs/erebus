clear variables;

% read in default parameter choices and runtime options
par_EREBUS

% set model run identifier
CTX.IO.RunID           = 'slug_1';
CTX.IO.TryContinue     =  0;

CTX.IO.LivePlot        = 'ON';
CTX.IO.nwrite          =  20;

CTX.FE.nx              =  20;
CTX.FE.nz              =  200;
CTX.FE.ElType          =  'Q2Q1';        

CTX.FE.W               =  20;    
CTX.FE.D               =  200;   
CTX.FE.LagrMesh        =  'SRF';

CTX.INIT.PertSmooth    =  round(15/(CTX.FE.D/CTX.FE.nz)^2); 

CTX.PHYS.TauPhi        =  1*60;    
CTX.PHYS.TauTmp        =  4*3600;
CTX.PHYS.delta         =  1;    
CTX.PROP.Phi0          =  0.30;                
CTX.PHYS.InflowRate    =  0.10;
CTX.PHYS.InflowPeriod  =  18*60;

CTX.INIT.TopoMode      =  'lake';
CTX.INIT.TopoHeight    =  0;
CTX.INIT.TopoWidth     =  4;
CTX.INIT.TopoXLoc      =  CTX.FE.W/2;

CTX.INIT.MatHeight     =  15;
CTX.INIT.MatWidth      =  5;
CTX.INIT.MatXLoc       =  CTX.FE.W/2;
CTX.INIT.MatZLoc       =  CTX.INIT.MatHeight;
CTX.PROP.CntAur        =  10;                           

CTX.INIT.AddBlock      =  5;                    
CTX.INIT.BlockWidth    =  CTX.INIT.MatWidth/2*0.75.*ones(1,5);
CTX.INIT.BlockHeight   =  CTX.INIT.MatHeight.*ones(1,5);
CTX.INIT.BlockZLoc     =  CTX.FE.D-([0,3,6,9,12]+0.5)*CTX.INIT.MatHeight;
CTX.INIT.BlockXLoc     =  CTX.INIT.MatXLoc.*ones(1,5);
CTX.INIT.BlockMat      =  3.*ones(1,5);

CTX.INIT.AddSphere     =  0;                    
CTX.INIT.SphereRadiusX =  CTX.INIT.BlockWidth/2;
CTX.INIT.SphereRadiusZ =  CTX.INIT.BlockWidth/2;
CTX.INIT.SphereZLoc    =  CTX.INIT.MatZLoc+0*CTX.INIT.MatHeight;
CTX.INIT.SphereXLoc    =  CTX.INIT.MatXLoc;
CTX.INIT.SphereMat     =  3;

% run simulation code
run('../src/EREBUS.m')
