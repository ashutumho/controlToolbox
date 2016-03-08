
function output = realp(varargin)
    [lhs,rhs]=argn(0)
    if rhs == 0 then
        error(msprintf(gettext("Function has no input argument..")))
    elseif rhs >=3 | rhs == 1 then
        error(msprintf(gettext("Incorrect number of input arguments.")))
    end
    varbName = varargin(1)
    varbData = varargin(2)
    if typeof(varbName) <> 'string' then
        error((gettext("Wrong type for argument: String expected.")))//msprintf
     elseif typeof(varbData) <> 'constant' then
         error((gettext("Wrong type for argument: Real matrix (2D) expected.")))//msprintf
    end
    varbPhase = phasemag(varbData)
    findImag = find(varbPhase <> 0 & varbPhase <> 180)
    if size(findImag,"r") ~= 0 then
        error((gettext("Wrong type for argument: Real matrix (2D) expected.")))
    end
    sizeData = size(varbData)
    numbOfElements = 1 
    for ii = 1 : 2
        numbOfElements = numbOfElements*sizeData(ii)
    end
    for ii = 1: numbOfElements
        minData(ii,1) = -%inf
        freeData(ii,1) = 1
    end
    maxData = -minData
    minData = hypermat([sizeData],minData)
    maxData = hypermat([sizeData],maxData)
    freeData = hypermat([sizeData],freeData)
    output("Name") = varbName
    output("Value") = varbData
    output("Minimum") = minData
    output("Maximum") = maxData
    output("Free") = freeData
    //output = varbData
    
    
endfunction
