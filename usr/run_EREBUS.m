clear variables;

% read in default parameter choices and runtime options
par_EREBUS

% set model run identifier
CTX.IO.RunID  = 'lavalake_rand';

% choose to plot results live during simulation,
% set interval of time steps to plot and store output
CTX.IO.LivePlot   = 'ON';
CTX.IO.nwrite     =  10;

% set numerical mesh size
CTX.FE.nx         =  100/1;
CTX.FE.nz         =  100/1;
CTX.FE.ElType     =  'Q2Q1';                   % finite-element type (Q1P0, Q1Q1, Q2Q1)

CTX.FE.W            =  40;                    % Width of model domain [m]
CTX.FE.D            =  40;                   % Depth of model domain [m]

CTX.PHYS.TauPhi      =  1*60;    
CTX.PHYS.TauTmp      =  2*3600;
CTX.PHYS.delta       =  1;    
CTX.PHYS.InflowRate  =  0.05;
CTX.PHYS.InflowPeriod = 1e32;

CTX.INIT.MatHeight  =  30;
CTX.INIT.MatWidth   =  8;
CTX.INIT.MatXLoc    =  CTX.FE.W/2;
CTX.INIT.MatZLoc    =  CTX.INIT.MatHeight;

CTX.INIT.AddBlock      =  1;                     % add block of material
CTX.INIT.BlockWidth    =  [CTX.INIT.MatWidth/2*0.6];
CTX.INIT.BlockHeight   =  [CTX.FE.D-CTX.INIT.MatHeight];
CTX.INIT.BlockZLoc     =  [CTX.FE.D/2+CTX.INIT.MatHeight/2];
CTX.INIT.BlockXLoc     =  [CTX.INIT.MatXLoc];
CTX.INIT.BlockMat      =  [3];

CTX.INIT.AddSphere     =  1;                     % add sphere/ellipse of material
CTX.INIT.SphereRadiusX =  [CTX.INIT.BlockWidth/2];
CTX.INIT.SphereRadiusZ =  [CTX.INIT.BlockWidth/2];
CTX.INIT.SphereZLoc    =  [CTX.INIT.MatZLoc];
CTX.INIT.SphereXLoc    =  [CTX.INIT.MatXLoc];
CTX.INIT.SphereMat     =  [3];

% run simulation code
run('../src/EREBUS.m')
