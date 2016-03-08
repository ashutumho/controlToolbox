
function [poles , index] =dsort(varargin)
//sort the discrete domain pole by magnitude
//
//Calling Sequence
//poldsort = dsort(poles)
//[plodsort indexnumb] = dsort(poles)
//
//Parameters
//poles :  a matrix of double, represents the  poles of a discrete domain
//poldsort : a matrix of double, represents the  sorted poles of given discrete domain
//indexnumb :  a matrix of double, represents the indices of given discrete domain poles
//
//Description
//
//poldsort = dsort(ploes) sort the discrete domain poles in decreasing order  by its magnitude.
//
// [plodsort indexnumb] = dsort(poles) sort the discrete domain poles and also gives the indices of the given pole
//
//Examples
//poles = [-4 ; 1-%i ; 2 ; 1+%i;-3;-5;-2-%i;-2+%i ]
//poldsort = dsort(poles)
//[plodsort indexnumb] = dsort(poles)
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
        magSort = []
        magIndex = 1
    else
        [magSort magIndex]= gsort(abs(tempSysData(:,1)),'g','d')
        for ii = 1:numbOfElements
            magSort(ii) = tempSysData(magIndex(ii),1)
        end
    end
    //disp(magSort)
    poles = magSort
    index = magIndex
endfunction
