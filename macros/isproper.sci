function output = isproper(varargin)
    
    [lhs,rhs]=argn(0)
    if rhs == 0 then
        error(msprintf(gettext("Function has no input argument..")))
    elseif rhs >=3 then
        error(msprintf(gettext("Incorrect number of input arguments.")))
    end
    if rhs == 2 then
        if strcmp(varargin(2),"elem") == 1 then
            error(msprintf(gettext("Wrong type for argument #%d: String expected.")))
        end
    end
    
    sysData = varargin(1)
    elemData = varargin(2)
    
    select typeof(sysData)
    case "st" then
    case "rational" then
    case "state-space" then
    else
        error(msprintf(gettext("Incompatible input argument.")))
    end
    if typeof(sysData) == 'rational' then
        sizeData = size(sysData)
        lengthData = length(sizeData)
        numbOfElements = 1
        for ii = 1:lengthData
            numbOfElements = numbOfElements*sizeData(ii)
        end
        numData = sysData.num
        denData = sysData.den
        for ii = 1:numbOfElements
            tempNumCoeff = coeff(numData(ii))
            tempDenCoeff = coeff(denData(ii))
            numLength = length(tempNumCoeff)
            denLength = length(tempDenCoeff)
            if numLength > denLength then
                tempOutput(ii,1) = 0
            elseif numLength <= denLength then
                tempOutput(ii,1) = 1 
            end
        end
        if rhs == 1 then
            tempData = find(tempOutput == 0)
            if size(tempData,"r") == 0 then
                output = 0
             else
                 output = 1
             end
        elseif rhs == 2 then
            if lengthData <= 2 then
                output = hypermat([sizeData],tempOutput(:,1))
            elseif lengthData > 2 then
                tempData = hypermat([sizeData])
                tempData(:,:,:,:) = tempOutput(:,1)
                output = tempData
            end
        end
    elseif typeof(sysData) == 'state-space' then
        output = 1
    elseif typeof(sysData) == 'st' then
        output = 1
    end
endfunction
