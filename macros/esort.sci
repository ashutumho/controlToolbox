
function [poles , index] =esort(varargin)
//sort the continuous domain pole by real part
//
//Calling Sequence
//polesort = esort(poles)
//[ploesort indexnumb] = esort(poles)
//
//Parameters
//poles :  a matrix of double, represents the  poles of a continuous domain
//polesort : a matrix of double, represents the  sorted poles of given continuous domain
//indexnumb :  a matrix of double, represents the indices of given continuous domain poles
//Description
//polesort = esort(ploes) sort the continuous domain poles in decreasing order  by its real part.
// [ploesort indexnumb] = esort(poles) sort the continuous domain poles and also gives the indices of the given pole
//Examples
//poles = [-4 ; 1-%i ; 2 ; 1+%i;-3;-5;-2-%i;-2+%i ]
//polesort = esort(poles)
//[ploesort indexnumb] = esort(poles)
//Authors
// Ashutosh Kumar Bhargava

    [lhs,rhs]=argn(0)
    // check the input elements
    if rhs == 0 then
        error(msprintf(gettext("Function has no input argument..")))
    elseif rhs >=2 then
        error(msprintf(gettext("Incorrect number of input arguments.")))
    end
    // storing the data type
    sysData = varargin($)
    select typeof(sysData)
    case "constant" then
    case "hypermat" then
    else
        error(msprintf(gettext("Incompatible input argument.")))
    end
    // storing the dimension
    sizeData = size(sysData)
    lengthData = length(sizeData)
    numbOfElements = 1
    for ii =1:lengthData
        numbOfElements = numbOfElements*sizeData(ii)
    end
    //storing the number in row vector
    for ii = 1:numbOfElements
        tempSysData(ii,1) = sysData(ii)
        //tempSysData(ii,2) = ii
        if phasemag(sysData(ii)) ~= 0 & phasemag(sysData(ii)) ~= 180 & phasemag(sysData(ii)) ~= 90 & phasemag(sysData(ii)) ~= -90 then
            poleIndex = find(sysData == sysData(ii))
            conjpoleIndex = find(sysData == conj(sysData(ii)))
            if length(poleIndex) ~= length(conjpoleIndex) then
                error(msprintf(gettext("Poles must be appear in complex conjugate form.")))
            end
        end
    end
    if size(sysData,"r") == 0 then
        realSort = []
        realIndex = 1
    else
        [realSort realIndex]= gsort(real(tempSysData(:,1)),'g','d')
        for ii = 1:numbOfElements
            tempComp = imag(tempSysData(realIndex(ii),1))*%i
            realSort(ii) = realSort(ii)+ tempComp
        end
    end
    phs = phasemag(realSort,'c')
    if numbOfElements > 1 then
        for ii = 1:numbOfElements-1
            if (phs(ii) == -phs(ii+1)) & (phs(ii) < phs(ii+1)) then
                temprealSort = realSort(ii)
                realSort(ii) = realSort(ii+1)
                realSort(ii+1) = temprealSort
                temprealIndex = realIndex(ii)
                realIndex(ii) = realIndex(ii+1)
                realIndex(ii+1) = temprealIndex
            elseif phs(ii) == -phs(ii+1) then
                ii = ii+1
            end
        end
        
    end
    poles = realSort
    index = realIndex
endfunction
    
