function output =isfinite(varargin)
//Checks whether the system has finite coefficients
//
//Calling Sequence
//bool = isfinite(sys)
//bool = isfinite(sys,'elem')
//bool = isfinite(mat)
//
//Parameters
//sys : represents the valid form of state space model or transfer function
//mat : a matrix of double, represents the array, matrix or hyper-matrix
//bool : a boolean, represents the output
//elem : a string, represents the output data presentation    
//
//Description
//bool = isfinite(sys), the output will be 1 when all the coefficient  of the transfer function or state-space is finite.
//bool = isfinite(sys,'elem'), the dimension of the output will be same as the dimension of the sys and each array elements represent the logical response of the respective transfer function or state-space
//bool = isfinite(mat), it checks all the elements of the matrix and produces the output with the same dimension as the matrix. 
//
//Examples
// mat = [1 2 %inf 0 -%inf]; bool = isfinite(mat)
//s = poly(0,'s') ; sys = syslin('c',1/(s+1)); bool = isfinite(sys)
//sys = =ssrand(2,3,8,list('st',2,5,5)); bool = isfinite(sys)
//bool = isfinite(sys,'elem') 
//
//Authors
// Ashutosh Kumar Bhargava

    [lhs,rhs]=argn(0)
    if rhs == 0 then
        error(msprintf(gettext("Function has no input argument..")))
    elseif rhs >=3 then
        error(msprintf(gettext("Incorrect number of input arguments.")))
    end
    if rhs == 2 then
        if strcmp(varargin(2),"elem") == 1 then//varargin(2) == "elem" then
            error(msprintf(gettext("Wrong type for argument #%d: String expected.")))
        end
    end
    sysData = varargin(1)
    select typeof(sysData)
    case "constant" then
        if rhs == 2 & strcmp(varargin(2),"elem") == 1 then
            error(msprintf(gettext("Incorrect number of input arguments.")))
        end
    case "hypermat" then
        if rhs == 2 & strcmp(varargin(2),"elem") == 1 then
            error(msprintf(gettext("Incorrect number of input arguments.")))
        end
    case "rational" then
    case "state-space" then
    else
        error(msprintf(gettext("Incompatible input argument.")))
    end
    sizeData = size(sysData)
    lengthData = length(sizeData)
    numbOfElements = 1
    for ii =1:lengthData
        numbOfElements = numbOfElements*sizeData(ii)
    end
    // test the transfer function or state-space function
    if typeof(sysData)=='rational' | typeof(sysData)=='state-space' then
        if typeof(sysData) == 'state-space' then
            sysData = ss2tf(sysData)
        end
        numData = sysData.num
        denData = sysData.den
        for ii = 1:numbOfElements
            tempNumData = coeff(numData(ii))
            tempDenData = coeff(denData(ii))
            numInfData = find(tempNumData==%inf)
            denInfData = find(tempDenData==%inf)
            if size(numInfData,"r")~=0 | size(denInfData,"r") ~= 0 then
                tempOutput(ii,1) = %F
            else
                tempOutput(ii,1) = %T
            end
        end
    // test the constant or hyper matrix type data     
    elseif typeof(sysData) == 'constant' | typeof(sysData) == 'hypermat' then
        if size(sysData,"r") == 0 then
            output = []
            tempOutput = []
            break
        else
            tempData(:,1) = sysData
            for ii = 1:numbOfElements
                if tempData(ii,1) == %inf then
                    tempOutput(ii,1) = 0
                else
                    tempOutput(ii,1) = 1
                end
            end
        end
    end
    // generate the output
    if rhs == 1 then
        tempData = find(tempOutput == %F)
        // when data type is Rational or State-space type
        if typeof(sysData)=='rational' | typeof(sysData)=='state-space' then
            if size(tempData,"r") == 0 then 
                output = 1
            else
                output = 0
            end
        // when data type is Constant or Rational type    
        elseif typeof(sysData) == 'constant' | typeof(sysData) == 'hypermat' then
            if lengthData <= 2 then
                output = hypermat([sizeData],tempOutput(:,1))
            elseif lengthData > 2 then 
                tempData = hypermat([sizeData])
                tempData(:,:,:,:) = tempOutput(:,1)
                output = tempData
            end
        end
    // when data type is Rational or State-spcae type and user want to check each coeff of array    
    elseif rhs == 2 then
        if lengthData <= 2 then
            output = hypermat([sizeData],tempOutput(:,1))
        elseif lengthData > 2 then 
            tempData = hypermat([sizeData])
            tempData(:,:,:,:) = tempOutput(:,1)
            output = tempData
        end
    end

endfunction
