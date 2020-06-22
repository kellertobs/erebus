% ReportTiming    EREBUS subroutine to report model timing at each time step
%
% []  =  ReportTiming(CTX)
%
%   modified  20170427  Tobias Keller
%   modified  20200227   Tobias Keller


function  []  =  ReportTiming(CTX)

time = CTX.TIME.total;
step = CTX.TIME.step;
spmn = 60;
sphr = 60*60;
spdy = 60*60*24;
spyr = 60*60*24*365.25;

if time < spmn
    tunit = ' sec';
elseif time >=  spmn && time < sphr
    tunit = ' min';
    time  = time/spmn;
elseif time >=  sphr && time < spdy
    tunit = ' hr';
    time  = time/sphr;
elseif time >=  spdy && time < spyr
    tunit = ' day';
    time  = time/spdy;
elseif time >=  spyr && time < 1e3*spyr
    tunit = ' yr';
    time  = time/spyr;
elseif time >=  1000*spyr && time < 1e6*spyr
    tunit = ' kyr';
    time  = time/1e3/spyr;
else
    tunit = ' Myr';
    time  = time/1e6/spyr;
end

if step < spmn
    sunit = ' sec';
elseif step >=  spmn && step < sphr
    sunit = ' min';
    step  = step/spmn;
elseif step >=  sphr && step < spdy
    sunit = ' hr';
    step  = step/sphr;
elseif step >=  spdy && step < spyr
    sunit = ' day';
    step  = step/spdy;
elseif step >=  spyr && step < 1e3*spyr
    sunit = ' yr';
    step  = step/spyr;
elseif step >=  1000*spyr && step < 1e6*spyr
    sunit = ' kyr';
    step  = step/1e3/spyr;
else
    sunit = ' Myr';
    step  = step/1e6/spyr;
end

fprintf('***** ')
fprintf(1,' %i',CTX.TIME.istep)
fprintf(';   time = ')
fprintf(1,'%e %s',time,tunit)
fprintf(';   timestep = ')
fprintf(1,'%e %s.\n\n',step,sunit)

end