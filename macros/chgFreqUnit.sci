function output =chgFreqUnit(varargin)
//Convert frequency unit of the frequency domain data
//
//Calling Sequence
//sysfreq = chgFreqUnit(sys,changedfrequnit)
//
//Parameters
//sys : a matrix of double, represents the frequency response points of systems
//changedfrequnit: a string , represents the changed frequency unit . The new frequency units are rad/Timeunit,cycles/TimeUnit,rad/s,kHz,MHz,GHz,rpm,. The default frequency unit is rad/TimeUnit.
//sysfreq : a matrix of double, represents the update frequency response points of systems
//
//Description
//
// sysfreq = chgFreqUnit(sys,changedfrequnit) updated frequency unit of the sys. The response at corresponding point will be same.
//
//Examples
//sys = frd(1:10,1:10)
//sysfreq = chgFreqUnit(sys,'Khz')
//sysfreq1 = chgFreqUnit(sys,'rpm')
//
//Authors
//Ashutosh Kumar Bhargava

    [lhs,rhs]=argn(0)
    // check the input elements
//    if size(varargin(1),"r") == 0 | size(varargin(2),"r") == 0 then
//        error(msprintf(gettext("Incompatible input argument.")))
    if rhs == 0 then
        error(msprintf(gettext("Function has no input argument..")))
    elseif rhs > 2 | rhs == 1 then
        error(msprintf(gettext("Incorrect number of input arguments.")))
    end
    sysData = varargin(1)
    freqData = varargin(2)
    select typeof(sysData)
    case "st" then
    else
        error(msprintf(gettext("Incompatible input argument.")))
    end
    if typeof(freqData)~='string' then
         error(msprintf(gettext("Incompatible input argument.")))
    elseif strcmp(freqData,'rad/TimeUnit') == 0 | strcmp(freqData,'rad/TimeUnits') == 0 then
        newFreq = 1
    elseif strcmp(freqData,'cycles/TimeUnit') == 0 | strcmp(freqData,'cycles/TimeUnits') == 0 then
        newFreq = 2*%pi
    elseif strcmp(freqData,'rad/s') == 0 then
        newFreq = 1
    elseif strcmp(freqData,'Hz') == 0 | strcmp(freqData,'Hzs') == 0 then
        newFreq = 2*%pi
    elseif strcmp(freqData,'kHz') == 0 | strcmp(freqData,'kHzs') == 0 then
        newFreq = 2*%pi*10^3
    elseif strcmp(freqData,'MHz') == 0 | strcmp(freqData,'MHzs') == 0 then
        newFreq = 2*%pi*10^6
    elseif strcmp(freqData,'GHz') == 0 | strcmp(freqData,'MHzs') == 0 then
        newFreq = 2*%pi*10^9
    elseif strcmp(freqData,'rpm') == 0 then
        newFreq = 2*%pi/60
    else
        error(msprintf(gettext("specified frequency units is rad/TimeUnit,cycles/TimeUnits,rad/s,Hz,kHz,MHz,GHz,rpm.")))
    end
    sysData.freq = sysData.freq/newFreq
    output( "freq" ) = sysData.freq
    output( "resp" ) = sysData.resp
    printf('\tfreq\t resp\n')
    printf('\t------ \t ------\n')
    output = [sysData.freq' sysData.resp']
endfunction
