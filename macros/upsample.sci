function h = upsample(varargin)
//upsample the input signal data, discrete domain transfer functions and state-space model
//
//Calling Sequence
//output = upsample(sys,N)
//output = upSample(sys,N,Phase)
//
//Parameters
//sys : a matrix of double or rational, represents the real valued sample data, discrete time transfer function and state-space model
//N : a positive integer , represents the upsample of the sys
//Phase: a positive integer , range will be 0 to N-1
//output: a matrix of double or rational, represents the upsampled sys
//
//Description
//
// output = upsample(sys,N) in the output N-1 zero will be interpolated in sys if it is a signal data. If sys is a discrete domain transfer function then the order of output transfer function will increase by N times and sampling time will be reduced by 1/N times. And if the sys is state space type then it will increase the number of states by N times as well as other matrix.
// 
//output = upSample(sys,N,Phase) it applied on the signal data. It adds the phase offset in the upsampled signal.
//
//Examples
//sys = [1 2;3 4]; output = upsample(sys,3)
//output = upsample(sys,3,2)
//sys =rand(2,2,2); output = upsample(sys,3)
//output = upsample(sys,3,2)
//sys =rand(1,2,2); output = upsample(sys,3)
//output = upsample(sys,3,2)
//z = poly(0,'z'); sys = syslin('d',(z+1)/(z^2+z+1));sys.dt = 0.1
//output = upsample(sys,3)
//output.dt 
//A = rand(3,3);B = rand(3,2);C = rand(2,3);D = zeros(2,2);
//sys = syslin('d',A,B,C,D);sys.dt = 0.1;
//output = upsample(sys,3)
//output.dt
//
//Authors
//Ashutosh Kumar Bhargava
//
    [lhs,rhs] = argn(0)
    if rhs == 0 then
        error(msprintf(gettext("Function has no input argument..")))
    elseif rhs == 1 | rhs > 3 then
        error(msprintf(gettext("Incorrect number of input arguments.")))
    end
    if typeof(varargin(1))<>["constant","rational","hypermat","state-space"] then
        error(msprintf(gettext("Incompatible input argument.")))
    elseif typeof(varargin(2)) <> "constant" then 
        error(msprintf(gettext("Resampling rate must be real positive integer ")))
    elseif (rhs == 3 & typeof(varargin(1)) == "state-space")|(rhs == 3 & typeof(varargin(1)) == "rational") then 
        error(msprintf(gettext("Too many input arguments.")))
    elseif rhs == 3 then
        if typeof(varargin(3)) <> "constant" | varargin(3) < 0 | varargin(3) == %inf then
            error(msprintf(gettext("Phase offset must be positive integer from 0 to N-1.(N = Upsample ) ")))
        end
    end
    upSamplingData = varargin(1)
    upSamplingStep = varargin(2)
    if upSamplingStep == 0 | upSamplingStep == %inf | upSamplingStep < 0 then
        error(msprintf(gettext("Resampling rate must be positive integer ")))
    end
    upSamplingStep = upSamplingStep -1
    dimData = size(upSamplingData)                                            // storing the dimension of the matrix
    numbOfElements = 1
    for ii = 1:length(dimData)
        numbOfElements = numbOfElements*dimData(ii)                           // storing the total number of elements
    end
    // upsampling the constant type of data
    if typeof(upSamplingData) == "constant" | typeof(upSamplingData) == "hypermat" then
        if rhs == 3 then
            upSamplingPhaseOffset = varargin(3)
        else
            upSamplingPhaseOffset = 0
        end
        //disp(numbOfElements)
//        disp(upSamplingData)
//        disp(upSamplingStep)
        //disp(upSamplingPhaseOffset)
        //for ii = 1 : numbOfElements
        if numbOfElements == 1 then
            if upSamplingPhaseOffset == 0 then
                tempData = zeros(upSamplingStep+1,1)
                tempData(1,1) = upSamplingData
            else
                tempData = zeros(upSamplingStep+1,1)
                tempData(upSamplingPhaseOffset+1,1) = upSamplingData
            end
        elseif length(dimData) == 2 & dimData(1) == 1 then
//            printf('\n ola \n')
//            disp(max(dimData))
            tempData = zeros(1,max(dimData)*(upSamplingStep+1))
            if upSamplingPhaseOffset == 0 then
                upSamplingPhaseOffset = 1
                tempPhaseOffset = 0
            else
                tempPhaseOffset = 1
            end
            tempIndex = 1
            for ii = upSamplingPhaseOffset:upSamplingStep+1:length(tempData)
                if tempPhaseOffset == 0 then
                    tempData(1,ii) = upSamplingData(1,tempIndex)
                else
                    tempData(1,ii+1) = upSamplingData(1,tempIndex)
                end
                tempIndex = tempIndex+1
            end
            //disp(tempData)
        else
//            numbOfElements = numbOfElements*(upSamplingStep+1)
//            disp(numbOfElements)
            tempData = zeros(numbOfElements*(upSamplingStep+1),1)
            tempNumber = 0
            tempIndex = 1
            for ii = 1:numbOfElements
                if upSamplingPhaseOffset == 0 then
                    tempData((ii-1)*(upSamplingStep+1)+1,1) = upSamplingData(ii) 
                elseif upSamplingPhaseOffset ~= 0 then
                    tempData((ii-1)*(upSamplingStep+1)+upSamplingPhaseOffset+1,1) = upSamplingData(ii) 
                end
            end
            if dimData(1) == 1 then
                dimData(2) = dimData(2)*(upSamplingStep+1)
            else
                dimData(1) = dimData(1)*(upSamplingStep+1)
            end
            tempData = hypermat([dimData],tempData)
            //disp(tempData)
        end
    elseif typeof(varargin(1)) == "rational" then
        z = poly(0,'z')
        samplingData = upSamplingData.dt
        //disp(samplingData)
        if samplingData == 'c' then
            error(msprintf(gettext("Upsample command works on discrete domain systems or matrix of sample data")))
        elseif samplingData == 'd' | samplingData == 0 then
            error(msprintf(gettext("Discrete time system must have proper sampling timing.")))
        elseif samplingData < 0 then
            error(msprintf(gettext("Discrete time system must have positive real value.")))
        end
        numData = upSamplingData.num
        denData = upSamplingData.den
        for ii = 1 : numbOfElements
            tempNumData = numData(ii)
            tempDenData = denData(ii)
            tempNumCoeff = coeff(tempNumData)
            tempDenCoeff = coeff(tempDenData)
            numNewCoeff = zeros(1,length(tempNumCoeff)*(upSamplingStep+1))
            denNewCoeff = zeros(1,length(tempDenCoeff)*(upSamplingStep+1))
            for jj = 1 : length(tempNumCoeff)
                numNewCoeff(1,(jj-1)*(upSamplingStep+1)+1) = tempNumCoeff(jj)
            end
            for jj = 1 : length(tempDenCoeff)
                denNewCoeff(1,(jj-1)*(upSamplingStep+1)+1) = tempDenCoeff(jj)
            end      
            tempNumData = poly(numNewCoeff,'z','coeff')
            tempDenData = poly(denNewCoeff,'z','coeff')
            tempData(ii,1) = tempNumData/tempDenData
        end
        samplingData = samplingData/(upSamplingStep+1)
        if length(dimData) <=2 then
            tempData = hypermat([dimData],tempData(:,1))
        elseif length(dimData) > 2 then
            tempData1 = hypermat([dimData])
            tempData1(:,:,:,:) = tempData(:,1)
            tempData = tempData1 
        end
        tempData.dt = samplingData
        disp(tempData)
    elseif typeof(varargin(1)) == "state-space" then
        samplingData = upSamplingData.dt
        //disp(samplingData)
        if samplingData == 'c' then
            error(msprintf(gettext("Upsample command works on discrete domain systems or matrix of sample data")))
        elseif samplingData == 'd' | samplingData == 0 then
            error(msprintf(gettext("Discrete time system must have proper sampling timing.")))
        elseif samplingData < 0 then
            error(msprintf(gettext("Discrete time system must have positive real value.")))
        end
        [Aa Bb Cc Dd] = abcd(upSamplingData)
        dimAa = size(Aa)
        //disp(dimAa(1,1))
        newDim = dimAa(1,1)*(upSamplingStep+1)
        tempAa = zeros(newDim,newDim)
        for ii = 1 : dimAa(1)
            tempAa((newDim-dimData(1,1)):newDim,ii) = Aa(:,ii)
        end
        tempIndex = 1
        for ii = (dimAa(1)+1):newDim
            tempAa(tempIndex,ii) = 1
            tempIndex = tempIndex+1
        end
        dimBb = size(Bb)
        newDim = dimBb(1,1)*(upSamplingStep+1)
        tempBb = zeros(newDim,dimBb(1,2))
        //disp((newDim-dimBb(1,1)+1))
        for ii = 1 : dimBb(1,2)
            tempBb((newDim-dimBb(1,1))+1:newDim,ii) = Bb(:,ii)
        end
        //disp(tempBb)
        dimCc = size(Cc)
        newDim = dimCc(1,2)*(upSamplingStep+1)
        tempCc = zeros(dimCc(1,1),newDim)
        for ii = 1 : dimCc(1,1)
            tempCc(ii,1:dimCc(1,2)) = Cc(ii,:)
        end
        //disp(tempCc)
        tempData = syslin('d',tempAa,tempBb,tempCc,Dd)
        tempData.dt = samplingData/(upSamplingStep+1)
        //disp(tempData)
                
    end
    h = tempData
endfunction
