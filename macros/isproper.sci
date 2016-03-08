function output = isproper(varargin)
//finds that the system is proper
//
//Calling Sequence
//bool = isproper(sys)
//bool  = isproper(sys,'elem')
//
//Parameters
//sys : represents matrix of transfer functions,SISO,MIMO systems and state-space functions
//bool : a boolean, represents the output
//elem: a string
//
//Description 
//
// bool = isproper(sys) the output will be true when the sys will be proper. A system is proper if the number of zero is less than or equal to its number of poles. For MIMO systems individual transfer function must be proper then the output will be true. 
//
//bool = isproper(sys,'elem') if the sys is array type then it will return  logical array of equal dimension of sys.
//
//Examples
//s = poly(0,'s'); sys = syslin('c',1/(s+1));bool = isproper(sys)
//sys = pid(rand(2,2),0,3); 
//bool = isproper(sys)
//bool = isproper(sys,'elem')
//sys =  ssrand(2,3,8,list('st',2,5,5)); bool = isproper(sys)
//
//Authors
//Ashutosh Kumar Bhargava     
//
//
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
    if rhs == 2 then
        elemData = varargin(2)
    end
    
    
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
