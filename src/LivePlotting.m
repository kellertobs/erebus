% LivePlotting    EREBUS subroutine to plot model results live during run
%
% []  =  LivePlotting(CTX)
%
%   Function creates figures and plots model solution and material point 
%   properties live during model run. 
%
%   modified  20170427  Tobias Keller
%   modified  20200227  Tobias Keller
%   modified  20200515  Tobias Keller


function  []  =  LivePlotting(CTX)

if strcmp(CTX.IO.LivePlot,'ON')
    
    if CTX.FE.nxEl >= CTX.FE.nzEl
        sp21 = 211;
        sp22 = 212;
        sp31 = 311;
        sp32 = 312;
        sp33 = 313;
    else
        sp21 = 121;
        sp22 = 122;
        sp31 = 131;
        sp32 = 132;
        sp33 = 133;
    end
    
    n=1;
    
    FileName  =  [CTX.IO.DataDir '/' CTX.IO.RunID '/' CTX.IO.RunID '_hist.mat'];
    load(FileName,'H');
    
    figure(n); n=n+1; clf;
    subplot(3,1,1);
    time  = H.time/60;
    area  = pi*(CTX.FE.W/2)^2;
    gasin = H.bot.GasIn*area; gasout = H.top.GasOut*area;
    plot(time,gasin,'b',time,gasout,'r'); axis tight;
    axis([min(time),max(time)+1e-12,min(min(gasin),min(gasout))./10,max(max(gasin),max(gasout)).*1.1+1e-12])
    title('Gas flux [kg/s]')
    subplot(3,1,2);
    heatin = H.bot.HeatIn*area/1e6; heatout = H.top.HeatOut*area/1e6;
    plot(time,heatin,'b',time,heatout,'r'); axis tight;
    axis([min(time),max(time)+1e-12,min(min(heatin),min(heatout))./10,max(max(heatin),max(heatout)).*1.1+1e-12])
    title('Heat flux [MW]')
    subplot(3,1,3);
    plot(time,H.top.Fsp(:,2),'b',time,H.top.Fsp(:,3),'r'); axis tight;
    axis([min(time),max(time)+1e-12,0,max(H.top.Fsp(:,3)).*1.1+1e-12])
    xlabel('Time [min]');
    title('Surface flow speed [m/s]')
    
    figure(n); n=n+1; clf;
    subplot(sp21);
    PlotField(CTX.MP.Mat,CTX.FE,CTX.IO.PlotStyle);
    title('Material types')
    subplot(sp22);
    PlotField(CTX.MP.Rho,CTX.FE,CTX.IO.PlotStyle);
    title('Density [kg/m3]')
    drawnow
    
    figure(n); n=n+1; clf;
    subplot(sp21);
    PlotField( CTX.SL.U,CTX.FE,[CTX.IO.PlotStyle,'U']);
    title('x-Velocity [m/s]')
    subplot(sp22);
    PlotField(-CTX.SL.W,CTX.FE,CTX.IO.PlotStyle);
    title('z-Velocity [m/s]')
    drawnow
    
    figure(n); n=n+1; clf;
    subplot(sp21);
    PlotField(CTX.SL.Phi.*100,CTX.FE,CTX.IO.PlotStyle);
    title('Vesicularity [vol %]')
    subplot(sp22);
    PlotField(min(CTX.RHEO.Chic,CTX.SL.Chi).*100,CTX.FE,CTX.IO.PlotStyle);
    title('Crystallinity [vol %]')
    drawnow
    
    figure(n); n=n+1; clf;
    subplot(sp21);
    PlotField(CTX.SL.GPhi,CTX.FE,[CTX.IO.PlotStyle,'U']);
    title('Degassing Rate [1/s]')
    subplot(sp22);
    PlotField(CTX.SL.GTmp,CTX.FE,CTX.IO.PlotStyle);
    title('Cooling Rate [deg C/s]')
    drawnow
    
    figure(n); n=n+1; clf;
    subplot(sp31);
    PlotField(CTX.MP.DivV,CTX.FE,CTX.IO.PlotStyle);
    title('Volume strain rate [log10 1/s]')
    subplot(sp32);
    PlotField(max(-6,log10(CTX.MP.EII(:,1))),CTX.FE,CTX.IO.PlotStyle);
    title('Shear strain rate [log10 1/s]')
    subplot(sp33);
    PlotField(max(-2,log10(CTX.MP.TII)),CTX.FE,CTX.IO.PlotStyle);
    title('Shear stress [log10 Pa]')
    drawnow
    
    figure(n); n=n+1; clf;
    subplot(sp21);
    PlotField(CTX.SL.P./1e6,CTX.FE,CTX.IO.PlotStyle);
    title('Dynamic pressure [MPa]')
    subplot(sp22);
    PlotField(CTX.SL.Pt./1e6,CTX.FE,CTX.IO.PlotStyle);
    title('Total pressure [MPa]')
    drawnow
    
    figure(n); n=n+1; clf;
    subplot(2,2,1);
    PlotField(log10(CTX.MP.Eta),CTX.FE,CTX.IO.PlotStyle);
    title('Eta V [Pas]')
    subplot(2,2,2);
    PlotField(log10(CTX.MP.EtaVP),CTX.FE,CTX.IO.PlotStyle);
    title('Eta VP [Pas]')
    subplot(2,2,3);
    PlotField((CTX.MP.Xi),CTX.FE,CTX.IO.PlotStyle);
    title('Chi VEP [1]')
    subplot(2,2,4);
    PlotField(log10(CTX.MP.EtaVEP),CTX.FE,CTX.IO.PlotStyle);
    title('Eta VEP [Pas]')
    drawnow
    
    figure(n); n=n+1; clf;
    subplot(sp21);
    PlotField(log10(CTX.MP.YieldStr),CTX.FE,CTX.IO.PlotStyle);
    title('Yield Stress [log10 Pa]')
    subplot(sp22);
    PlotField(CTX.MP.TII - CTX.MP.YieldStr,CTX.FE,CTX.IO.PlotStyle);
    title('Failure function (<=0)')
    drawnow
    
    figure(n); n=n+1; clf;
    subplot(2,2,1);
    PlotField(max(-6,log10(CTX.MP.EII(:,1))),CTX.FE,CTX.IO.PlotStyle);
    title('Dev. strainr total [1/s]')
    subplot(2,2,2);
    PlotField(max(-6,log10(CTX.MP.EII(:,2))),CTX.FE,CTX.IO.PlotStyle);
    title('Dev. strainr viscous [1/s]')
    subplot(2,2,3);
    PlotField(max(-6,log10(CTX.MP.EII(:,3))),CTX.FE,CTX.IO.PlotStyle);
    title('Dev. strainr elastic [1/s]')
    subplot(2,2,4);
    PlotField(max(-6,log10(CTX.MP.EII(:,4))),CTX.FE,CTX.IO.PlotStyle);
    title('Dev. strainr plastic [1/s]')
    drawnow
    
end

end