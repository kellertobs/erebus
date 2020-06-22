
%**************************************************************************
%**********  Set Parameters for EREBUS   **********************************
%**************************************************************************

%*****  Run options  ******************************************************

CTX.IO.SrcDir        =  '../src';                % source directory
CTX.IO.DataDir       =  '../out';                % output directory
CTX.IO.TryContinue   =  0;                       % try to continue prev. run

%*****  Solver options  ***************************************************

CTX.SL.maxits        =  1;                       % max non-linear iterations
CTX.SL.atol          =  1.e-9;                   % absolute convergence tolerance (precond. residual norm)
CTX.SL.rtol          =  1.e-2;                   % relative convergence tolerance
CTX.SL.Advection     =  'FRM';                  % advection method
CTX.SL.CFL           =  0.25;                     % Courant number for coarse time step
CTX.SL.RFN           =  10;                      % Factor for refined adv-diff time step

CTX.SL.theta_it      =  1.0;                     % iterative relaxation parameter
CTX.SL.theta_dt      =  0.5;                     % time-step weighting parameter (0.5 = Cranck-Nicolson; 0 = Forward Euler; 1 = Backward Euler)

CTX.SL.Smooth        =  0.25;                    % Smoothing factor to regularise viscosity
CTX.SL.StabFact      =  1.e-24;                  % Stabilisation factor for P-coefficient matrix
CTX.SL.RhoRef        =  0;                       % reference density (= 0 for full P, != 0 for dynamic P)


%*****  finite-element mesh and Lagrangian particles options  *************

CTX.FE.ElType        = 'Q2Q1';                   % finite-element type (Q1P0, Q1Q1, Q2Q1)
CTX.FE.nblock        =  12800;                   % block size for vectorised matrix assembly 

CTX.FE.nx            =  200;                     % number of elements in x-direction on FE mesh
CTX.FE.nz            =  300;                     % number of elements in z-direction on FE mesh

CTX.FE.D             =  100;                     % Depth of model domain [m]
CTX.FE.W             =   10;                     % Width of model domain [m]

CTX.FE.LagrMesh      =  'SRF';                   % switch for Lagrangian mesh (ON/OFF), or Lagrangian free surface (SRF)


%*****  Time stepping options  ********************************************

CTX.TIME.spyr        =  3600*24*365.25;          % seconds per year

CTX.TIME.step        =  0.01;                    % initial time step size [s]
CTX.TIME.end         =  2*3600;                  % stopping time for simulation run [s]


%*****  I/O and live plotting options  ************************************

CTX.IO.nwrite        =  10;                      % time steps between output frames
CTX.IO.LivePlot      = 'ON';                     % switch ON for live output plotting
CTX.IO.PlotStyle     = 'srf';                    % 'img' = image(), 'srf' = surface(), 
                                                 % '..._lr' = reflected left,
                                                 % '..._rr' = reflected right

%*****  Physical parameters  **********************************************

CTX.PHYS.grav        =  9.81;                    % gravity [m/s2]
CTX.PHYS.RConst      =  8.314;                   % universal gas constant [J/K/mol]
CTX.PHYS.TauPhi      =  1.5*60;                  % surface outgassing time scale [s]
CTX.PHYS.TauTmp      =  3*3600;                  % surface cooling time scale [s]
CTX.PHYS.delta       =  1;                       % surface outgassing/cooling layer depth [m]
CTX.PHYS.InflowRate  =  0.025;                   % bubbly magma inflow speed [m/s]
CTX.PHYS.InflowPeriod = 1e32;                    % bubbly magma inflow period [s] (set very high for constant inflow)


%*****  set initial condition for temperature field  **********************
CTX.INIT.TempMode    = 'lavalake';               % temperature initiation mode
CTX.INIT.TempInt     = 1000;                     % lake interior temperature
CTX.INIT.TempExt     = 100;                      % lake exterior temperature


%*****  initial conditions for material types  ****************************

CTX.INIT.MatMode       =  'lavalake';            % materials initiation mode
CTX.INIT.Mat           =  [1,2];                 % material types for lake geometry
CTX.INIT.MatHeight     =  CTX.FE.D/2;            % lake depth
CTX.INIT.MatWidth      =  8;                     % conduit width
CTX.INIT.MatXLoc       =  CTX.FE.W/2;            % horizontal position of conduit mouth
CTX.INIT.MatZLoc       =  CTX.FE.D/2;            % vertical position of base of lake

CTX.INIT.AddBlock      =  0;                     % add block of material
CTX.INIT.AddSphere     =  0;                     % add sphere/ellipse of material


%*****  Options for initial random perturbations  *************************

CTX.INIT.PertSmooth     =  round(25/(CTX.FE.nz/CTX.FE.D)^2);  % smoothing of random noise
CTX.INIT.PertSymmetric  =  'OFF';                % switch to symmetric noise distribution
CTX.INIT.PertFrict      =  0;                    % amplitude of friction angle perturbation
CTX.INIT.PertCoh        =  0;                    % amplitude of cohesion perturbation


%*****  Set material properties of material types  ************************

CTX.PROP.n      =  [      1;       2;       3];  % material ID numbers
CTX.PROP.Eta    =  [  1.e+4;   1.e+4;   1.e+4];  % material viscosity [Pas]
CTX.PROP.G      =  [  1.e+9;   1.e+9;   1.e+9];  % material shear modulus [Pa]
CTX.PROP.Coh    =  [  1.e+4;   1.e+4;   1.e+4];  % material cohesion [Pa]
CTX.PROP.Frict  =  [    0.5;     0.5;     0.5];  % material friction coefficient [1]

CTX.PROP.k      =  5;                            % thermal conductivity [W/m/K]
CTX.PROP.c      =  1200;                         % heat capacity [J/kg/K]
CTX.PROP.alpha  =  1e-4;                         % thermal expansivity [J/kg/K]
CTX.PROP.L      =  400e3;                        % latent heat of melting [J/kg]
CTX.PROP.Tsol   =  800;                          % solidus temperature [deg C]
CTX.PROP.Tliq   =  1100;                         % liquidus temperature [deg C]

CTX.PROP.Phi0   =  0.10;                         % initial vesicularity [vol]
CTX.PROP.kPhi   =  1e-5;                         % bubble diffusivity [m2/s]
CTX.PROP.Rho    =  2400;                         % melt density [kg/m3]
CTX.PROP.RhoChi =  2600;                         % crystal density
CTX.PROP.MPhi   =  0.025;                        % gas molar mass [kg/mol]
CTX.PROP.Compr  =  1;                            % 1 for compressible, 0 for incompressible gas


%*****  Options for visco-elastic/brittle-plastic rheology  ***************

CTX.RHEO.ConstViscosity  =  'OFF';               % switch for constant viscosity mode
CTX.RHEO.Plasticity      =  'ON';                % switch for plastic failure mode
CTX.RHEO.Elasticity      =  'OFF';               % switch for elastic mode

CTX.RHEO.Strainr0        =  5.e-3;               % reference strain rate [1/s]
CTX.RHEO.PowerLawExp     =  3;                   % non-Newtonian powerlaw exponent

CTX.RHEO.MaxEta          =  1.e+9;               % maximum cutoff viscosity [Pas]
CTX.RHEO.MinEta          =  1.e+2;               % minimum cutoff viscosity [Pas]
CTX.RHEO.MinGdt          =  1.e+2;               % minimum cutoff elastic strength [Pas]

CTX.RHEO.mChi            = 2;                    % crystal stiffening exponent
CTX.RHEO.Chic            = 0.50;                 % crystal close packing fraction
CTX.RHEO.mPhi            = -2;                   % bubble stiffening/softening exponent
CTX.RHEO.Phic            = 0.50;                 % bubble close packing fraction


%***** Boundary conditions for velocity, pressure  ************************

CTX.BC.FreeSurface     =  'ON';                 % switch for free surface top boundary

CTX.BC.Type            =  'ConstStr';            % specify constant strainrate or constant velocity boundaries
CTX.BC.BGStrainr       =   0e-15;                % specify boundary background strainrate [1/s]

CTX.BC.TopBot          =  'fs';                  % top/bot boundaries: free slip 'fs'; no slip 'ns'; free '00', bottom flux 'bt'
CTX.BC.Sides           =  'fs';                  % side boundaries   : free slip 'fs'; no slip 'ns'; simple shear 'ss'; free '00'

CTX.BC.SurfPres        =  1e+5;                  % surface pressure [Pa] (1 atm = 1e+5 Pa)






 
