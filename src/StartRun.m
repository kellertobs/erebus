% StartRun  EREBUS subroutnie to start new simulation run or continue from previous one
%
% [CTX] = StartRun(CTX)
%
%   Prepares data structures and sets initial conditions from input
%   options. Requires as input an application context 'CTX' produced with a
%   EREBUS parameter script
%
%   created   20161115  Tobias Keller
%   modified  20170427  Tobias Keller
%   modified  20190418  Tobias Keller
%   modified  20200227  Tobias Keller
%   modified  20200515  Tobias Keller


function  [CTX] = StartRun(CTX)

if CTX.IO.TryContinue == 1
    
    %*****  check input file  *********************************************
    
    if ~isfield(CTX.IO,'ContID')
        CTX.IO.ContID  =  CTX.IO.RunID;
    end
    
    if isfield(CTX.IO,'ContFrame') && CTX.IO.ContFrame > 0
        CTX.IO.ContName = [CTX.IO.DataDir '/' CTX.IO.ContID '/' CTX.IO.ContID '_' num2str(CTX.IO.ContFrame) '.mat'];
    else
        CTX.IO.ContName = [CTX.IO.DataDir '/' CTX.IO.ContID '/' CTX.IO.ContID '_cont.mat'];
    end
    
    id = exist(CTX.IO.ContName,'file');
    
    if id ~= 2
        CTX.IO.TryContinue = 0;
        fprintf(1,'!!!  Try continue failed: start run from scratch  !!!\n\n');
    else
        
        %*****  continue from output file  ********************************
        
        fprintf(1,'-------       continue from file      -------\n\n');
        load(CTX.IO.ContName);
        
        CTX.TIME.istep  =  CTX.TIME.istep + 1;
        CTX.IO.frame    =  CTX.IO.frame + 1;
        
        CTX.SL.S              =  zeros(2*CTX.FE.NU+CTX.FE.NP,1);
        CTX.SL.S(CTX.FE.DOFU) = CTX.SL.U;
        CTX.SL.S(CTX.FE.DOFW) = CTX.SL.W;
        CTX.SL.S(CTX.FE.DOFP) = CTX.SL.P;
        
    end
    
else
    
    %*****  create output directory if not present yet
    
    if ~exist([CTX.IO.DataDir '/' CTX.IO.RunID],'dir'); mkdir([CTX.IO.DataDir '/' CTX.IO.RunID]); end

    
    %*****  initialise timing / counting parameters  **********************
    
    CTX.TIME.total  =  0;
    CTX.TIME.istep  =  1;
    
    CTX.IO.frame    =  0;
    
    CTX.SL.it       =  0;
    CTX.SL.fnorm    =  1e6;
    
    
    %*****  initialise finite element grid  *******************************
    
    CTX.FE  =  InitFE(CTX.FE);
    CTX     =  InitTopo(CTX);
    CTX.FEo =  CTX.FE;
    
    
    %*****  initialise material distribution and properties  **************
    
    CTX      =  InitMat(CTX);
    CTX.MPo  =  CTX.MP;
    
    
    %*****  initialise random perturbation field  *************************
    
    CTX.MP.pert  =  InitPert(CTX.FE,CTX.INIT.PertSmooth,CTX.INIT.PertSymmetric);
    
    
    %*****  initialise solution variables  ********************************
    
    CTX      =  InitSL(CTX);
    CTX.SLo  =  CTX.SL;
    
    
    %*****  update material properties and aux. fields  *******************
    
    CTX  =  UpdateMaterialPoints(CTX);
    
    
    %*****  Initialise solution vector and scaling matrix  ****************
    
    CTX.SL.S  =  zeros(2*CTX.FE.NU+CTX.FE.NP,1);
    
    
    %*****  plot initial condition and save to file  **********************
    
    if strcmp(CTX.IO.LivePlot,'ON'); LivePlotting(CTX); end
    SaveToFile(CTX);
    CTX.IO.frame  =  1;
    
    
end

end

