% StokesSolver    EREBUS subroutine to calculate solution to Stokes equations
%
% [CTX] = StokesSolver(MP,CTX)
%
%   Function finds solution to linear governing equations. It operates
%   according to options set in CTX.SL; it calls AssembleOperator() and
%   ApplyBC() to prepare the coefficient matrix and boundary conditions,
%   then uses "\" operator to solve linear system of equations.
%
%   created 20161115 Tobias Keller
%   modified  20170427  Tobias Keller
%   modified  20170508  Tobias Keller
%   modified  20200227  Tobias Keller
%   modified  20200515  Tobias Keller


function  [CTX] = StokesSolver(CTX)

SL  =  CTX.SL;
MP  =  CTX.MP;
FE  =  CTX.FE;
S   =  SL.S;


%*****  assemble operator and right hand side vector  *********************

[L,RHS]     =  AssembleOperator(MP,SL,CTX);

X           =  max(1e-32,min(1e+32,sqrt(abs(diag(L)))));
X           =  diag(sparse(1./X));

L           =  X*L*X;
RHS         =  X*RHS;

[L,RHS,BC]  =  ApplyBC(L,RHS,X,CTX);


%*****  get non-linear residual  ******************************************

F           =  zeros(size(S));
F(BC.free)  =  L*(S(BC.free)./diag(X(BC.free,BC.free))) - RHS(BC.free);
F(BC.ind)   =  0.;

SL.F.U      =  F(FE.DOFU);
SL.F.W      =  F(FE.DOFW);
SL.F.P      =  F(FE.DOFP);

SL.fnorm    =  norm(F(BC.free),2)/(norm(RHS(BC.free),2)+1.e-16);

if SL.it == 1 || SL.fnorm > SL.fnorm0
    SL.fnorm0  =  SL.fnorm;
end

fprintf(1,'    abs =')
fprintf(1,' %4.4e',SL.fnorm)
fprintf(1,'  rel =')
fprintf(1,' %4.4e\n',SL.fnorm./SL.fnorm0)

if SL.fnorm < SL.atol && SL.it > 1
    CTX.SL  =  SL;
    CTX.BC  =  BC;
    return
end

if SL.fnorm > 1.e+20 && SL.it > 1
    SaveToFile(CTX);
    error('!!!  Error: Diverging solution, stop and try again  !!!')
end


%*****  solve linear system of equations  *********************************

S           =  zeros(size(S));
S(BC.free)  =  L \ RHS(BC.free);
S(BC.ind)   =  BC.val;

SL.S        =  X*S;

SL.U        =  SL.S(FE.DOFU);
SL.W        =  SL.S(FE.DOFW);
SL.P        =  SL.S(FE.DOFP);

% SL.P(SL.P~=0) =  SL.P(SL.P~=0) - SL.P(FE.MapPn(1,1));
SL.Pref     =  CTX.SL.RhoRef.*CTX.PHYS.grav.*FE.CoordP(:,2) + CTX.BC.SurfPres;
SL.Pt       =  SL.P  + SL.Pref;


CTX.SL  =  SL;
CTX.BC  =  BC;

end

