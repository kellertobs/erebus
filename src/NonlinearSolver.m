% NonlinearSolver    EREBUS subroutine to solve coupled non-linear Stokes
%                    and time-dependent problems
%
% [CTX] = NonlinearSolver(CTX)
%
%   Function finds solution to non-linear governing equations. It operates
%   according to options set in CTX.SL; it iteratively calls LinearSolver()
%   to update the Stokes velocity-pressure solution, then AdvectFE() to advect
%   the FE mesh to conform to free surface deformation, and TDSolver() to solve
%   the time-dependent sub-problem, before non-linear coefficients and auxiliary 
%   fields are updated according to the new solution guess in UpdateMaterialPoints().
%
%   created   20161115  Tobias Keller
%   modified  20170427  Tobias Keller
%   modified  20170508  Tobias Keller
%   modified  20200227  Tobias Keller
%   modified  20200515  Tobias Keller


function  [CTX] = NonlinearSolver(CTX)
    
% prepare structures
CTX.SLo = CTX.SL;
CTX.MPo = CTX.MP;
CTX.FEo = CTX.FE;
CTX.TIME.stepo = CTX.TIME.step;

% initialise variables
CTX.SL.it     =  0;
CTX.SL.fnorm  =  1;
fnorm0        =  1;

% main non-linear iteration loop
while (  CTX.SL.fnorm        > CTX.SL.atol   ... % convergence criteria
      && CTX.SL.fnorm/fnorm0 > CTX.SL.rtol   ...
      && CTX.SL.it           < CTX.SL.maxits )
    
    CTX.SL.it = CTX.SL.it+1;  % update iteration count
    
    % report iteration count
    if CTX.SL.it < 10
        fprintf(1,'    %i',CTX.SL.it)
    else
        fprintf(1,'   %i',CTX.SL.it)
    end
    
    % solve Stokes velocity-pressure problem
    CTX     =  StokesSolver(CTX);
    
    % update velocity-dependent time step
    CTX     =  CalcTimestep(CTX);
    
    % advect finite-element mesh to trace free surface
    if ~strcmp(CTX.FE.LagrMesh,'OFF')
        CTX =  AdvectFE(CTX);
    end
    
    % solve time-dependent problem
    CTX     =  TDSolver(CTX);
    
    % update solution-dependent fields
    CTX     =  UpdateMaterialPoints(CTX);
    
    % update baseline residual norm
    if CTX.SL.it == 1 || CTX.SL.fnorm > fnorm0
        fnorm0 = CTX.SL.fnorm;
    end
    
end

% update record of surface evolution (positive for height above surface)
CTX.FE.SurfEvo =  [CTX.FE.SurfEvo,-CTX.FE.SurfQ2];


