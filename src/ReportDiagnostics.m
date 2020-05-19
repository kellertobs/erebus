% ReportDiagnostics    EREBUS subroutine to report model diagnostics at
%                      each time step
%
% []  =  ReportDiagnostics(CTX)
%
%   Function reports model diagnostics by printing to standard output.
%
%   modified  20170427  Tobias Keller
%   modified  20200227  Tobias Keller
%   modified  20200515  Tobias Keller


function  []  =  ReportDiagnostics(CTX)

Tc       =  CTX.RHEO.Chic .* (CTX.PROP.Tsol - CTX.PROP.Tliq) + CTX.PROP.Tliq;
lake     =  CTX.SL.T > Tc;
lakeIP   =  logical(PQ2IP(lake,CTX.FE));

fprintf('    max z-speed  = ')
fprintf(1,'%1.4f'  ,max(abs(CTX.SL.W(lake))));
fprintf('    max x-speed  = ')
fprintf(1,'%1.4f\n',max(abs(CTX.SL.U(lake))));

fprintf('    max bubbles  = ')
fprintf(1,'%1.4f'  ,max(CTX.SL.Phi(lake)));
fprintf('   mean bubbles  = ')
fprintf(1,'%1.4f'  ,mean(CTX.SL.Phi(lake)));
fprintf('    min bubbles  = ')
fprintf(1,'%1.4f\n',min(CTX.SL.Phi(lake)));

fprintf('    max crystals = ')
fprintf(1,'%1.4f'  ,max(CTX.SL.Chi(lake)));
fprintf('   mean crystals = ')
fprintf(1,'%1.4f'  ,mean(CTX.SL.Chi(lake)));
fprintf('    min crystals = ')
fprintf(1,'%1.4f\n',min(CTX.SL.Chi(lake)));

fprintf('    max temp.    = ')
fprintf(1,'%4.1f'  ,max(CTX.SL.T(lake)));
fprintf('   mean temp.    = ')
fprintf(1,'%4.2f'  ,mean(CTX.SL.T(lake)));
fprintf('    min temp.    = ')
fprintf(1,'%4.2f\n\n',min(CTX.SL.T(lake)));

fprintf('    max eta      = ')
fprintf(1,'%1.4e',max(CTX.MP.EtaVEP(lakeIP)));
fprintf('    min eta      = ')
fprintf(1,'%1.4e\n',min(CTX.MP.EtaVEP(lakeIP)));

fprintf('    max stress   = ')
fprintf(1,'%1.4e',max(CTX.MP.TII(lakeIP,1)));
fprintf('    min stress   = ')
fprintf(1,'%1.4e\n',min(CTX.MP.TII(lakeIP,1)));

fprintf('    max strainr  = ')
fprintf(1,'%1.4e',max(CTX.MP.EII(lakeIP,1)));
fprintf('    min strainr  = ')
fprintf(1,'%1.4e\n\n',min(CTX.MP.EII(lakeIP,1)));

fprintf('    max El def.  = ')
fprintf(1,'%1.4f',CTX.FE.MaxDef );
fprintf('    min El def.  = ')
fprintf(1,'%1.4f\n',CTX.FE.MinDef );
fprintf('    max El vol.  = ')
fprintf(1,'%1.4f',max(CTX.FE.ElVol)/mean(CTX.FE.ElVol) );
fprintf('    min El vol.  = ')
fprintf(1,'%1.4f\n\n',min(CTX.FE.ElVol)/mean(CTX.FE.ElVol) );

end