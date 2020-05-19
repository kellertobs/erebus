% PrintDiagnostics    EREBUS subroutine to report model diagnostics
%
% []  =  PrintDiagnostics(CTX)
%
%   Function reports model diagnostics by printing to standard output
%
%   modified  20170427  Tobias Keller
%   modified  20200227  Tobias Keller
%   modified  20200515  Tobias Keller


function  []  =  ReportDiagnostics(CTX)

fprintf('    max z-speed    = ')
fprintf(1,'%1.3e'  ,max(abs(CTX.SL.W))*CTX.TIME.spyr);
fprintf('    max x-speed    = ')
fprintf(1,'%1.3e\n',max(abs(CTX.SL.U))*CTX.TIME.spyr);

fprintf('    max bubble f.  = ')
fprintf(1,'%1.3e'  ,max(CTX.SL.Phi));
fprintf('    min bubble f.  = ')
fprintf(1,'%1.3e\n',min(CTX.SL.Phi));

fprintf('    max crystal f. = ')
fprintf(1,'%1.3e'  ,max(CTX.SL.Chi));
fprintf('    min crystal f. = ')
fprintf(1,'%1.3e\n',min(CTX.SL.Chi));

fprintf('    max temp.      = ')
fprintf(1,'%1.3e'  ,max(CTX.SL.T));
fprintf('    min temp.      = ')
fprintf(1,'%1.3e\n\n',min(CTX.SL.T));

fprintf('    max eta      = ')
fprintf(1,'%1.3e',max(CTX.MP.EtaVEP));
fprintf('    min eta      = ')
fprintf(1,'%1.3e\n',min(CTX.MP.EtaVEP));

fprintf('    max stress   = ')
fprintf(1,'%1.3e',max(CTX.MP.TII(:,1)));
fprintf('    min stress   = ')
fprintf(1,'%1.3e\n',min(CTX.MP.TII(:,1)));

fprintf('    max strainr  = ')
fprintf(1,'%1.3e',max(CTX.MP.EII(:,1)));
fprintf('    min strainr  = ')
fprintf(1,'%1.3e\n\n',min(CTX.MP.EII(:,1)));

fprintf('    max El def.  = ')
fprintf(1,'%1.3f',CTX.FE.MaxDef );
fprintf('    min El def.  = ')
fprintf(1,'%1.3f\n',CTX.FE.MinDef );
fprintf('    max El vol.  = ')
fprintf(1,'%1.3f',max(CTX.FE.ElVol)/mean(CTX.FE.ElVol) );
fprintf('    min El vol.  = ')
fprintf(1,'%1.3f\n\n',min(CTX.FE.ElVol)/mean(CTX.FE.ElVol) );

end