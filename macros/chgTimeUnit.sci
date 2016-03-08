
function output = chgTimeUnit(varargin)
//Convert time unit of the system
//
//Calling Sequence
//systime = chgFreqUnit(sys,changedtimeunit)
//
//Parameters
//sys : a matrix of double, represents the frequency response points of systems
//changedtimeunit: a string , represents the changed time unit . The new time units are nanoseconds, microseconds, milliseconds, minutes, hours, days, weeks, months, years. The default time unit is second.
//systime : a matrix of double, represents the update system 
//
//Description
//
// systime = chgatimeUnit(sys,changedtimeunit) update time  unit of the sys. The response of the system will be same.
//
//Examples
//s = poly(0,'s'); z = poly(0,'z');
//sys = syslin('c',(s+2)/(s^2+3*s+1))
//systime = chgTimeUnit(sys,'hour')
//systime1 = chgTimeUnit(sys,'nanosecond')
//
//Authors
//Ashutosh Kumar Bhargava
 

        [lhs,rhs]=argn(0)
    // check the input elements
    s = poly(0,'s')
    if size(varargin(1),"r") == 0 | size(varargin(2),"r") == 0 then
        error(msprintf(gettext("Incompatible input argument.")))
    elseif rhs == 0 then
        error(msprintf(gettext("Function has no input argument..")))
    elseif rhs > 2 | rhs == 1 then
        error(msprintf(gettext("Incorrect number of input arguments.")))
    end
    sysData = varargin(1)
    timeData = varargin(2)
    select typeof(sysData)
    case "rational" then
    case "state-space" then
    else
        error(msprintf(gettext("Incompatible input argument.")))
    end
    
    if strcmp(timeData,'nanosecond') == 0 | strcmp(timeData,'nanoseconds') == 0 then
        newTime = 10^-09
    elseif strcmp(timeData,'microsecond') == 0 | strcmp(timeData,'microseconds') == 0 then
        newTime = 10^-6
    elseif strcmp(timeData,'millisecond') == 0 | strcmp(timeData,'milliseconds') == 0 then
        newTime = 10^-3
    elseif strcmp(timeData,'second') == 0 | strcmp(timeData,'seconds') == 0 then
        newTime = 1
    elseif strcmp(timeData,'minute') == 0 | strcmp(timeData,'minutes') == 0 then
        newTime = 60
    elseif strcmp(timeData,'hour') == 0 | strcmp(timeData,'hours') == 0 then
        newTime = 3600
    elseif strcmp(timeData,'day') == 0 | strcmp(timeData,'days') == 0 then
        newTime = 86400
    elseif strcmp(timeData,'week') == 0 | strcmp(timeData,'weeks') == 0 then
        newTime = 604800
    elseif strcmp(timeData,'month') == 0 | strcmp(timeData,'months') == 0 then
        newTime = 2592000
    elseif strcmp(timeData,'year') == 0 | strcmp(timeData,'years') == 0 then
        newTime = 31536000
    else
        error(msprintf(gettext("specified time units is nanoseconds, microseconds, milliseconds, seconds, minutes, hours, days, weeks, months, years .")))
    end
    dimData = size(sysData)
    lengthData = length(dimData)
    numbOfElements = 1
    for ii = 1:lengthData
        numbOfElements = numbOfElements*dimData(ii)
    end
    if typeof(sysData) == 'rational' then
        if sysData.dt == 'c' then
            numSys = sysData.num
            denSys = sysData.den
            for jj = 1:numbOfElements
                tempNum = numSys(jj)
                tempDen = denSys(jj)
                tempNumCoeff = coeff(tempNum)
                tempDenCoeff = coeff(tempDen)
                tempNumLength = length(tempNumCoeff)
                tempDenLength = length(tempDenCoeff)
                indexNumb = max(tempNumLength,tempDenLength)
                if tempNumLength == 1 & tempDenLength == 1 then // constant 
                    tempSysData(jj,1) = tempNum/tempDen
                elseif tempNumLength == tempDenLength then // degree of num = degree of den 
                    for ii = 1:indexNumb-1
                        tempNumCoeff(ii) = tempNumCoeff(ii)*newTime^(tempNumLength-ii)
                        tempDenCoeff(ii) = tempDenCoeff(ii)*newTime^(tempDenLength-ii)
                    end
                elseif tempNumLength > tempDenLength then //degree of num> degree of den
                    for ii = 1:indexNumb-1
                        tempNumCoeff(ii) = tempNumCoeff(ii)*newTime^(tempNumLength-ii)
                    end
                    for ii = 1:tempDenLength
                        tempDenCoeff(ii) = tempDenCoeff(ii)*newTime^(indexNumb-ii)
                    end
                elseif tempNumLength < tempDenLength then // degree of num < degree of den
                    for ii = 1:indexNumb-1
                        tempDenCoeff(ii) = tempDenCoeff(ii)*newTime^(tempDenLength-ii)
                    end
                    for ii = 1:tempNumLength
                        tempNumCoeff(ii) = tempNumCoeff(ii)*newTime^(indexNumb-ii)
                    end                    
                end
                tempNumCoeff = poly(tempNumCoeff,'s',"coeff")
                tempDenCoeff = poly(tempDenCoeff,'s',"coeff")
                tempSysData(jj,1) = tempNumCoeff/tempDenCoeff
            end
            //disp(tempSysData)
            if lengthData <= 2 then
                tempStor = hypermat([dimData],tempSysData(:,1))
            elseif lengthData > 2 then
                tempStor = hypermat([dimData])
                tempStor(:,:,:,:) = tempSysData
            end
            tempStor.dt = 'c'
            output = tempStor
        else
            if sysData.dt == 'd' then
                sysData.dt = 1
            end
            sysData.dt = sysData.dt/newTime
            output = sysData
        end
    elseif typeof(sysData) == 'state-space' then
        if sysData.dt == 'c' then
            [Aa,Bb,Cc,Dd] = abcd(sysData)
            Aa = newTime*Aa
            Bb = newTime*Bb
            output = syslin('c',Aa,Bb,Cc,Dd)
        else
            if sysData.dt == 'd' then
                sysData.dt = 1
            end
            sysData.dt = sysData.dt/newTime
            output = sysData
        end
    end
    //disp(tempSysData)
    //output = sysData
endfunction
