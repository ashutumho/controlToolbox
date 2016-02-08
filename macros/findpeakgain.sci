//find the peak gain and frequency
//given frequency must be in "rad/sec"
function [mag,frq]=findPeakGain(sys,fmin,fmax)
    [lhs,rhs]=argn(0)
    if and(typeof(sys)<>[ "rational" "state-space" ]) then
        error(msprintf(gettext("%s: Wrong type for input argument #%d: Linear dynamical system expected.\n"),"findPeakGain",1))
    end
    eps=1.e-7;
    
    select typeof(sys)
    case "rational" then
    case "state-space" then
        sys=ss2tf(sys)
    else
        error(97,1),
    end;
//store the data given from transfer function
    [freqc,repc] = repfreq(sys)
    [fase,ddb] = phasemag(repc,'c') // c is swap by m
    if rhs == 3 then
        fmin = fmin/(2*%pi)
        fmax = fmax/(2*%pi)
        j = find(freqc>= fmin & freqc<= fmax)
        ddb = ddb(j)
        freqc = freqc(j)
    end        
    magindb = max(ddb)
    i = find(ddb == magindb)
    i = min(i)
 if freqc(i)< eps then
     frqval = 0;
 else
     frqval = freqc(i)
 end
    mag = 10^(magindb/20)
    frq = 2*%pi* frqval
endfunction
