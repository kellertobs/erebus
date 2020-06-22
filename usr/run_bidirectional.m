clear variables;

% read in default parameter choices and runtime options
par_EREBUS

% set model run identifier
CTX.IO.RunID           = 'conduit';
CTX.IO.TryContinue     =  1;

CTX.IO.LivePlot        = 'ON';
CTX.IO.nwrite          =  20;

CTX.FE.nx              =  30/2;
CTX.FE.nz              =  300/2;
CTX.FE.ElType          =  'Q2Q1';        

CTX.FE.W               =  20;    
CTX.FE.D               =  200;   
CTX.FE.LagrMesh        =  'SRF';

CTX.INIT.PertSmooth    =  round(15/(CTX.FE.D/CTX.FE.nz)^2); 

CTX.PHYS.TauPhi        =  0.5*60;    
CTX.PHYS.TauTmp        =  6*3600;
CTX.PHYS.delta         =  1;    
CTX.PROP.Phi0          =  0.10;                
CTX.PHYS.InflowRate    =  0.0125;
CTX.PHYS.InflowPeriod  =  1e32;

CTX.INIT.MatHeight     =  15;
CTX.INIT.MatWidth      =  5;
CTX.INIT.MatXLoc       =  CTX.FE.W/2;
CTX.INIT.MatZLoc       =  CTX.INIT.MatHeight;
CTX.PROP.CntAur        =  10;                           

CTX.INIT.AddBlock      =  1;                    
CTX.INIT.BlockWidth    =  CTX.INIT.MatWidth/2*0.75;
CTX.INIT.BlockHeight   =  CTX.FE.D-0*CTX.INIT.MatHeight;
CTX.INIT.BlockZLoc     =  CTX.FE.D/2+0*CTX.INIT.MatHeight/2;
CTX.INIT.BlockXLoc     =  CTX.INIT.MatXLoc;
CTX.INIT.BlockMat      =  3;

CTX.INIT.AddSphere     =  0;                    
CTX.INIT.SphereRadiusX =  CTX.INIT.BlockWidth/2;
CTX.INIT.SphereRadiusZ =  CTX.INIT.BlockWidth/2;
CTX.INIT.SphereZLoc    =  CTX.INIT.MatZLoc+0*CTX.INIT.MatHeight/2;
CTX.INIT.SphereXLoc    =  CTX.INIT.MatXLoc;
CTX.INIT.SphereMat     =  3;

% run simulation code
run('../src/EREBUS.m')
