
function [output, impdata ] =  pid(varargin)
//Generate parallel form of PID controller, verify the parallel form of PID controller
//
//Calling Seqence
//output = pid(Kp)
//output = pid(Kp,Ki)
//output = pid(Kp,Ki,Kd)
//output = pid(Kp,Ki,Kd,Tf)
//output = pid(Kp,Ki,Kd,Tf,Ts)
//output = pid(tf)
//[output impdata] = pid(Kp)
//[output impdata] = pid(.....,'Notes','...')
//Parameters
//Kp : a matrix of double, represents the  Proportional gain of PID controller
//Ki : a matrix of double, represents the  Integral gain of PID controller 
//Kd- a matrix of double, represents the  Derivative gain of PID controller
//Tf : a matrix of double, represents the Time constant of first order filter of PID controller 
//Ts: a double, represents the Sampling timing
//sys: represents the valid form of PID controller   
//Name : a string represent  name of  Control system
//UserData : a string represent UserData of control system
//Notes : a string represent Notes about the control system
//Description
//output = pid(Kp,Ki,Kd,Tf) generates the parallel form of PID in continuous domain. The mathematical equation is given below
//
//   <latex>
	//    \begin{eqnarray}
    //      \[
    //          output=K_{p}+K_{i}\frac{1}{s}+\frac{K_{d}s}{T_{f}s+1}
    //      \]
	//    \end{eqnarray}
	//   </latex>
//
//All the PID parameters (Kp,Ki,Kd,& Tf) must be real number and output is the transfer function obtained by the given PID parameters data. 
//output = pid(Kp,Ki,Kd,Tf,Ts)  generates the parallel form of PID in discrete domain and Ts is the sampling timing of the controller. The mathematical equation is given below
//
//   <latex>
	//    \begin{eqnarray}
    //      \[
    //          output=K_{p}+K_{i}Iformula+\frac{K_{d}}{T_{f}+Dformula}
    //      \]
	//    \end{eqnarray}
	//   </latex>
//
//There are numbers of way to select the sampling formula for discrete time integration and derivation (Iformula,dformula respectively). By default value of the Iformula and Dformula is forward Euler. Other sampling formula is Backward Euler and Trapezoidal. If Dformula is forward Euler and first order derivative filter time constant Tf~=0 then it must satisfy the Tf >Ts/2 for proper derivative filter pole. The following way to select the sample method
//
//F - Forward Euler
//
//B - Backward Euler
//
//T - Trapezoidal   
//
//output = pid(Kp) generates proportional controller in continuous domain
//
//output = pid(0,Ki) generates integral controller in continuous domain
//
//output = pid(0,0,Kd) generates derivative controller without derivative filter
//
//output = pid(0,0,Kd,Tf) generates derivative filer with derivative filter of Tf time constant
//
//output = pid(Kp,Ki) generates proportional+integral controller in continuous domain
//
//output = pid(Kp,0,Kd) generates proportional + derivative controller without first order filter in continuous domain
//
//output = pid(Kp,0,Kd,Tf) generates proportional + derivative controller with first order filter in continuous domain
//
//output = pid(Kp,Ki,Kd) generates proportional +integral+ derivative controller without first order filter in continuous domain
//
//output = pid(Kp,Ki,Kd,Tf) generates proportional +integral+ derivative controller with first order filter in continuous domain
//
//output = pid(Kp,0,0,0,Ts) generates proportional controller in discrete domain
//
//output = pid(0,Ki,0,0,Ts) generates integarl controller in discrete domain
//
//output = pid(0,0,Kd,0,Ts) generates derivative without first order filter controller in discrete domain
//
//output = pid(0,0,Kd,Tf,Ts) generates derivative controller with first order filter in discrete domain
//
//output = pid(Kp,Ki,0,0,Ts) generates proportiona+Integral  controller in discrete domain
//
//output = pid(Kp,0,Kd,0,Ts) generates proportional+derivative controller without first order filter in discrete domain
//
//output = pid(Kp,0,Kd,Tf,Ts) generates proportional+derivative controller with first order filter in discrete domain
//
//output = pid(Kp,Ki,Kd,Tf,Ts) generates proportional+integral+derivative controller without first order filter in discrete domain
//
//output = pid(Kp,Ki,Kd,Tf,Ts) generates proportional+integral+derivative controller with first order filter in discrete domain
//
//output = pid(.......,Ts,'Iformula','__','Dformula','__') generates discrete domain controller with selection of sampling formula for Iformula and Dformula  
//
//output  = pid(sys) checks the given transfer function and verify it
//
//[output impdata] = pid(.....,'Notes','...','UserData','.....','Notes','....') it adds extra descriptive information about the controller in impdata 
//
// Examples
// output = pid(2)
// output = pid(2,3)
// output = pid(2,3,4)
// output = pid(2,3,4,5)
// output = pid(2,0,4)
// output = pid(2,0,4,5)
// output = pid(2,0,0,0,0.1)
// output = pid(2,3,0,0,0.1,'Iformula','B')
// output = pid(2,3,4,0,0.1,'Iformula','B','Dformula','T')
// output = pid(2,3,4,5,0.1,'Iformula','T','Dformula','B')
// s = poly(0,'s'); sys = syslin('c',3*(s+1)*(s+2)/s);
// [output impdata]= pid(sys)
//z = poly(0,'z'); sys = syslin(0.1,(19.146667 - 38.793333*z + 19.66*z^2)/(  1.6 - 3.4*z + 1.8*z^2))
// [output impdata] = pid(sys,'Iformula','T','Dformula','B')
// Authors
//
// Ashutosh Kumar Bhargava
//
//
//
        [lhs,rhs] = argn(0)
        s = poly(0,'s')
        z = poly(0,'z')
        funcprot(0)
//____________________________________________________________________________________________________________________________________________________________________________
        // find the type of data that are transferred
        if typeof(varargin(1))=='constant' | typeof(varargin(1))=='hypermat'  then
            datatype = 1
        elseif typeof(varargin(1))=='rational' then 
            datatype =2
        end
        impdata = tlist(["listtype","Kp","Ki","Kd","Tf","Iformula","Dformula","InputDelay","OutputDelay","Ts","Name","Notes","UserData"]) //
//____________________________________________________________________________________________________________________________________________________________________________
        // find the count
        count_numb = 0
        tempCount = min(rhs,5)
            for jj = 1:tempCount
                if typeof(varargin(jj)) == 'constant' | typeof(varargin(jj)) == 'rational' | typeof(varargin(jj)) == 'hypermat' then 
                    count_numb = count_numb + 1
                else
                    break
                end    
            end
//// Vrification of the syntax
////        syntax = ["Iformula","Dformula","Name","Notes","UseData"]
//        syntax = ['Iformula','Dformula','Name','Notes','UseData']
//        disp(rhs)
//        disp(tempCount)
//        if rhs > 5 then
//            for ii = 6:rhs
//                disp(ii)
//                printf('\n****')
//                if ii> rhs then
//                    error(msprintf(gettext("Please enter data in right formate")))
//                end
//                numbCount = 0
//                tempData = varargin(ii)
//                for jj = 1:5
//                    compString = isequal(syntax(jj),tempData)
//                    if compString == %T then
//                        ii = ii + 1
//                        printf('\n shaktimaan')
//                    elseif compString == %F then
//                        numbCount = numbCount+1
//                    end
//                end
//                disp(ii)
//                printf('\n---')
//                disp(numbCount)
//            end
//            
//        end
//            
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        if count_numb <= 4 | count_numb == 5 & varargin(5) == 0 then
            inctd = 1                                                           // in countinuous time domain
            indtd = 0                                                           // in discrete time domain
        elseif count_numb == 5 & varargin(5) ~= 0 then
            indtd = 1
            inctd = 0    
        end
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        kk = 0
        qq = 0
        //disp(count_numb)
        for ii=1:count_numb
            // check dimension
            //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if size(varargin(ii))==[1 1] then
                qq=qq+1
                if qq == count_numb then
                    current_length =1
                    current_size = [1 1]
                end
            else
                if kk==0 then 
                    temp_size = size(varargin(ii))
                    temp_length = length(temp_size)
                    kk=kk+1
                end
                    current_size = size(varargin(ii))
                    current_length = length(current_size)
                if temp_length ~= current_length then
                        error(msprintf(gettext("Please enter a right dimensions")))
                else
                    max_size = max(current_length,temp_length)
                    for jj =1:max_size
                         if current_size(jj)~=temp_size(jj) then
                             error(msprintf(gettext("Please enter a right dimensions")))
                         end
                    end
                end
            end
           //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    //******************************************************************
        //find the total number of array elements
        
        total_num_element = 1
        for ii=1:current_length
            total_num_element = total_num_element*current_size(ii)
        end
        // data transfer
        varargin_data=varargin
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        if datatype == 1 & inctd == 1 & indtd == 0 then                                                             //continuous time domain count_numb <= 4 &
            printf('\n Continuous Domain')
            sysTs = 'c'
            for ii =1:total_num_element
            //##########################################################
                for jj =1:4
                    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if jj<=count_numb then
                        if typeof(varargin(jj))=='hypermat' then
                            if jj<=rhs then
                                pid_data(jj,1)= varargin_data(jj).entries(ii)
                            elseif jj>rhs then
                                pid_data(jj,1)=0
                            end
                        elseif typeof(varargin_data(jj))=='constant' then
                            if size(varargin_data(jj))== [1 1] then
                                pid_data(jj,1)=varargin_data(jj)
                            else
                                temp_data = varargin_data(jj)
                                pid_data(jj,1)= temp_data(ii)
                            end
                        end
                    elseif jj>count_numb then
                            pid_data(jj,1)= 0
                    end
                //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end//end of data extraction
                // checking the imaginary and negative values
                if phasemag(pid_data(1,1),'c') ~= 0 & phasemag(pid_data(1,1),'c') ~= 180 | pid_data(1,1) == %inf then
                     error(msprintf(gettext("Proportional gain (Kp) must be a finite real number.")))
                elseif phasemag(pid_data(2,1),'c') ~= 0 & phasemag(pid_data(2,1),'c') ~= 180 | pid_data(2,1) == %inf then
                    error(msprintf(gettext("Integral gain (Ki) must be a finite positive real number.")))
                elseif phasemag(pid_data(3,1),'c') ~= 0 & phasemag(pid_data(3,1),'c') ~= 180 | pid_data(3,1) == %inf then
                    error(msprintf(gettext("Derivative gain (Kd) must be a finite positive real number.")))
                elseif phasemag(pid_data(4,1),'c') ~= 0 & phasemag(pid_data(4,1),'c') ~= 180 | pid_data(4,1) == %inf then
                    error(msprintf(gettext("Time constant of the first-order derivative filter(Tf) must be a finite positive real number.")))
                end            
//                disp(pid_data(1,1))
//                disp(pid_data(2,1))
//                disp(pid_data(3,1))
//                disp(pid_data(4,1))
//                disp(pid_data(3,1)+pid_data(1,1)*pid_data(4,1))
//                disp(pid_data(1,1)+pid_data(2,1)*pid_data(4,1))
                tempNum = (pid_data(3,1)+pid_data(1,1)*pid_data(4,1))*s^2+(pid_data(1,1)+pid_data(2,1)*pid_data(4,1))*s+pid_data(2,1)   // numerator of the parallel PID 
                tempDen = pid_data(4,1)*s^2+s                                                                                           // denominator of the parallel PID
                simp_mode(%t)
                tempTF = tempNum/tempDen
                TFdata(ii,1) = syslin('c',tempTF)                               //formation of the continuous time transfer function
                simp_mode(%f)
                temp_Kp_data(ii,1) = pid_data(1,1)
                temp_Ki_data(ii,1) = pid_data(2,1)
                temp_Kd_data(ii,1) = pid_data(3,1)
                temp_Tf_data(ii,1) = pid_data(4,1)
            end//end of continuous time PID
                // storing the parallel PID data in hyper matrix
                if current_length<=2 then
                    output_data = hypermat([current_size],TFdata(:,1))
                elseif current_length > 2 then
                    output_data = hypermat([current_size])
                    output_data(:,:,:,:)=TFdata
                end
                output = output_data
//=======================================================================================================================================================================
        elseif datatype == 1 & inctd == 0 & indtd == 1 then                                                                        //discrete time domain
                printf('\n Discrete time system\n')
                if size(varargin_data(5))~=[1 1] then
                    error(msprintf(gettext("Please enter valid sampling timing. " )))
                end
//**********************************************************************************************************************************************************************
                //now find the type of discrete time integral and derivative formula
                Iformula_index = 0
                Dformula_index = 0
                if rhs >5 then
                    for ii = 6:rhs
                        if varargin_data(ii)=='Iformula' then
                            Iformula_index = ii
                        elseif varargin_data(ii)=='Dformula' then
                            Dformula_index = ii
                        end
                    end
                end
                //validity of sampling formula
                if Iformula_index~= 0 then
                    if Iformula_index+1 > rhs then
                        error(msprintf(gettext(" Please select a proper sampling formula")))
                    elseif varargin_data(Iformula_index+1)~='F'& varargin_data(Iformula_index+1)~= 'B'& varargin_data(Iformula_index+1)~= 'T' then
                        error(msprintf(gettext("not valid sampling method, Please select \n F for Forward Euler Sampling \n B for Backward Euler Sampling \n T for Trapezoidal Sampling ")))
                    end
                end
                if Dformula_index~= 0 then
                    if Dformula_index+1 > rhs then
                        error(msprintf(gettext(" Please select a proper sampling formula")))
                    elseif varargin_data(Dformula_index+1)~='F' & varargin_data(Dformula_index+1)~='B' & varargin_data(Dformula_index+1)~='T' then
                        error(msprintf(gettext("not valid sampling method, Please select \n F for Forward Euler Sampling \n B for Backward Euler Sampling \n T for Trapezoidal Sampling ")))
                    end
                end
//**********************************************************************************************************************************************************************
                 //section of integral formula
                if Iformula_index==0 | varargin_data(Iformula_index+1)=='F' then
                    Iformula = varargin_data(5)/(z-1)
                    Iformula_method = 'F'
                elseif varargin_data(Iformula_index+1)=='B' then
                    Iformula = varargin_data(5)*z/(z-1)
                    Iformula_method = 'B'
                elseif varargin_data(Iformula_index+1)=='T' then
                    Iformula = varargin_data(5)*(z+1)/(2*(z-1))
                    Iformula_method = 'T'
                else//if Iformula_index+1>rhs|varargin_data(Iformula_index+1)~='T'|varargin_data(Iformula_index+1)~='B'|varargin_data(Iformula_index+1)~='F' then
                    Iformula = varargin_data(5)/(z-1)
                    Iformula_method = 'F'
                end
                
                 //selection of derivative formula
                if Dformula_index==0 |varargin_data(Dformula_index+1)=='F' then
                    Dformula_method = 'F'
                elseif varargin_data(Dformula_index+1)=='B' then
                    Dformula = varargin_data(5)*z/(z-1)
                    Dformula_method = 'B'
                elseif varargin_data(Dformula_index+1)=='T' then
                    Dformula = varargin_data(5)*(z+1)/(2*(z-1))
                    Dformula_method = 'T'
                else
                    Dformula_method = 'F'
                end
                
                // if Dformula is Farword Euler type then checking the following parameter
                if Dformula_method == 'F' then
                    temp_Tf = varargin_data(4)
                    temp_Kd = varargin_data(3)
                    size_Tf = size(temp_Tf)
                    if temp_Kd ~= 0 & temp_Tf ~= 0 then
                        if size_Tf == [1 1] then
                            temp_index = 1
                        else
                            temp_index = total_num_element 
                        end
                        for ii =1:temp_index
                            if varargin_data(5)> 2*temp_Tf(ii) then
                                printf('Index number = %d ,ii')
                                error(msprintf(gettext("Sample timing must be greater then two times of filter timing Ts>2*Tf")))
                            end
                        end
                     end
                    Dformula = varargin_data(5)/(z-1)
                end    
//*************************************************************************************************************************************************************************
        //now find the type of discrete time integral and derivative formula
            for ii =1:total_num_element
            //##########################################################
                for jj =1:4
                    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if jj<=count_numb then
                        if typeof(varargin(jj))=='hypermat' then
                            if jj<=rhs then
                                pid_data(jj,1)= varargin_data(jj).entries(ii)
                            elseif jj>rhs then
                                pid_data(jj,1)=0
                            end
                        elseif typeof(varargin_data(jj))=='constant' then
                            if size(varargin_data(jj))== [1 1] then
                                pid_data(jj,1)=varargin_data(jj)
                            else
                                temp_data = varargin_data(jj)
                                pid_data(jj,1)= temp_data(ii)
                            end
                        end
                    elseif jj>count_numb then
                        if jj == 4 then 
                            pid_data(jj,1)=%inf
                        else
                            pid_data(jj,1)= 0
                        end
                    end
                //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
                // checking the imaginary and negative values
                if phasemag(pid_data(1,1),'c') ~= 0 & phasemag(pid_data(1,1),'c') ~= 180 | pid_data(1,1) == %inf then
                     error(msprintf(gettext("Proportional gain (Kp) must be a finite real number.")))
                elseif phasemag(pid_data(2,1),'c') ~= 0 & phasemag(pid_data(2,1),'c') ~= 180 | pid_data(2,1) == %inf then
                    error(msprintf(gettext("Integral gain (Ki) must be a finite real number.")))
                elseif phasemag(pid_data(3,1),'c') ~= 0 & phasemag(pid_data(3,1),'c') ~= 180 | pid_data(3,1) == %inf then
                    error(msprintf(gettext("Derivative gain (Kd) must be a finite real number.")))
                elseif phasemag(pid_data(4,1),'c') ~= 0 & phasemag(pid_data(4,1),'c') ~= 180 | pid_data(4,1) == %inf then
                    error(msprintf(gettext("Derivative filter divisor(Tf) must be a finite real number.")))
                elseif phasemag(varargin_data(5)) ~= 0 | varargin_data(5) == %inf then
                    error(msprintf(gettext("Time constant of the first-order derivative filter (Ts) must be a finite positive real number.")))   
                end
                tempKi = Iformula*pid_data(2,1)
                tempKd = pid_data(3,1)/(pid_data(4,1)+Dformula)
                simp_mode(%t)
                tempTF = pid_data(1,1)+tempKi+tempKd
                TFdata(ii,1) = syslin('d',tempTF)                               //formation of the discrete time transfer function
                simp_mode(%f)
                temp_Kp_data(ii,1) = pid_data(1,1)
                temp_Ki_data(ii,1) = pid_data(2,1)
                temp_Kd_data(ii,1) = pid_data(3,1)
                temp_Tf_data(ii,1) = pid_data(4,1)                
            end
                if current_length<=2 then
                    output_data = hypermat([current_size], TFdata(:,1))
                elseif current_length > 2 then
                    output_data = hypermat([current_size])
                    output_data(:,:,:,:)= TFdata
                end
                output_data.dt = varargin_data(5)
                sysTs = output_data.dt
                output = output_data
//=================================================================================================================================================================================
            elseif datatype == 2 then
                printf('\n rational data type')
// storing the dimension of the transfer matrix in current_size
                    current_size = size(varargin_data(1))
                    current_length = length(current_size)
// total number of elements
                    total_num_element = 1
                    for ii=1:current_length
                        total_num_element = total_num_element*current_size(ii)
                    end
                    systf_num = varargin_data(1).num
                    systf_den = varargin_data(1).den
                    sysTs = varargin_data(1).dt
                    disp(sysTs)
                if varargin_data(1).dt == 'c' then                                   // continuous domain
                    for ii = 1:total_num_element
//------------------------------------------------------------------------------
// finding out the transfer function 
                    //trans_func = systf(ii)
                    num_coeff = coeff(systf_num(ii))
                    den_coeff = coeff(systf_den(ii))
                    
                    length_num_coeff = length(num_coeff)
                    length_den_coeff = length(den_coeff)
                    last_num_coeff = num_coeff(length_num_coeff)
                    last_den_coeff = den_coeff(length_den_coeff)
                    num_coeff = num_coeff/last_num_coeff
                    den_coeff = den_coeff/last_den_coeff
                    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                                      //----------------------------------------------------
                        if length_num_coeff == 1 & length_den_coeff == 1 then
                            temp_Kp_data(ii,1) = last_num_coeff*num_coeff(1)/last_den_coeff
                            temp_Ki_data(ii,1) = 0
                            temp_Kd_data(ii,1) = 0
                            temp_Tf_data(ii,1) = 0
                        //I type of controller
                        elseif length_num_coeff == 1 & length_den_coeff == 2 & den_coeff(1) == 0 then
                            temp_Kp_data(ii,1) = 0
                            temp_Ki_data(ii,1) = last_num_coeff*num_coeff(1)/last_den_coeff
                            temp_Kd_data(ii,1) = 0
                            temp_Tf_data(ii,1) = 0
                        // D type of controller(without filter)
                        elseif length_num_coeff == 2 & length_den_coeff == 1 & num_coeff(1) == 0 then     
                            temp_Kp_data (ii,1) = 0
                            temp_Ki_data (ii,1) = 0
                            temp_Kd_data (ii,1) = last_num_coeff
                            temp_Tf_data(ii,1) = 0
                        // D type controller with filter
                        elseif length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1) == 0 then
                            temp_Tf_data(ii,1) = 1/den_coeff(1)
                            temp_Kd_data(ii,1) = temp_Tf_data(ii,1)*last_num_coeff/last_den_coeff
                            temp_Kp_data(ii,1) = 0
                            temp_Ki_data(ii,1) = 0
                        // PI type of controller
                        elseif length_num_coeff == 2 & length_den_coeff == 2 & den_coeff(1) == 0 then
                            temp_Kp_data(ii,1) = last_num_coeff/last_den_coeff
                            temp_Ki_data(ii,1) = temp_Kp_data(ii,1)*num_coeff(1)
                            temp_Kd_data(ii,1) = 0
                            temp_Tf_data(ii,1) = 0
                        // PD type of controller
                        elseif length_num_coeff == 2 & length_den_coeff == 1 then 
                            temp_Ki_data(ii,1) = 0
                            temp_Kd_data(ii,1) = last_num_coeff/last_den_coeff
                            temp_Kp_data(ii,1) = temp_Kd_data(ii,1)*num_coeff(1)
                            temp_Tf_data(ii,1) = 0
                        // PID type of controller
                        elseif length_num_coeff == 3 & length_den_coeff == 2 & den_coeff(1) == 0 then
                            temp_Kd_data(ii,1) = last_num_coeff/last_den_coeff
                            temp_Kp_data(ii,1) = num_coeff(2)*temp_Kd_data(ii,1)
                            temp_Ki_data(ii,1) = num_coeff(1)*temp_Kd_data(ii,1)
                            temp_Tf_data(ii,1) = 0
                        elseif length_num_coeff == 3 & length_den_coeff == 3 & den_coeff(1) == 0 then
                            temp_Tf_data(ii,1) = 1/den_coeff(2)
                            //disp(temp_Tf_data(ii,1))
                            A = [temp_Tf_data(ii,1) 0 1;
                                (1- temp_Tf_data(ii,1)*num_coeff(2)) temp_Tf_data(ii,1) -num_coeff(2);
                                -temp_Tf_data(ii,1)*num_coeff(1) 1 -num_coeff(1)]
                            b = [temp_Tf_data(ii,1)*last_num_coeff/last_den_coeff 0 0]
                            temp_pid = A\b'
                            temp_Kp_data(ii,1) = temp_pid(1)
                            temp_Ki_data(ii,1) = temp_pid(2)
                            temp_Kd_data(ii,1) = temp_pid(3)
                        else
                            disp(varargin_data_index(ii,1))
                            error(msprintf(gettext("Above equation is not in parallel PID form")))
                        end
                    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                    end
                    // end of continuous domain
                elseif varargin_data(1).dt == 'd' | typeof(varargin_data(1).dt)=='constant' then                //discrete time domain
                    if varargin_data(1).dt == 'd' then
                        varargin_data(1).dt = 1                                 // If the sampling time is not given then the code will take 1 sec as sampling timing 
                    end
                    
//------------------------------------------------------------------------------
        // find the Iformula and Dformula index
                    Iformula_index = 0
                    Dformula_index = 0
                    for ii = 1:rhs
                        temp_data = varargin_data(ii)
                        I_index = isequal(temp_data,'Iformula')
                        D_index= isequal(temp_data,'Dformula')
                        if I_index== %T then
                            Iformula_index = ii
                        elseif D_index == %T then
                            Dformula_index = ii
                        end
                    end
            // end of finding of Iformula and Dformula
//------------------------------------------------------------------------------
                //validity of sampling formula
                if Iformula_index~= 0 then
                    if Iformula_index+1 > rhs then
                        error(msprintf(gettext(" Please select a proper sampling formula")))
                    elseif varargin_data(Iformula_index+1)~='F'& varargin_data(Iformula_index+1)~= 'B'& varargin_data(Iformula_index+1)~= 'T' then
                        error(msprintf(gettext("not valid sampling method, Please select \n F for Forward Euler Sampling \n B for Backward Euler Sampling \n T for Trapezoidal Sampling ")))
                    end
                end
                if Dformula_index~= 0 then
                    if Dformula_index+1 > rhs then
                        error(msprintf(gettext(" Please select a proper sampling formula")))
                    elseif varargin_data(Dformula_index+1)~='F' & varargin_data(Dformula_index+1)~='B' & varargin_data(Dformula_index+1)~='T' then
                        error(msprintf(gettext("not valid sampling method, Please select \n F for Forward Euler Sampling \n B for Backward Euler Sampling \n T for Trapezoidal Sampling ")))
                    end
                end
                
                 //section of integral formula
                if Iformula_index==0 | varargin_data(Iformula_index+1)=='F' then
                    Iformula = 'F'
                elseif varargin_data(Iformula_index+1)=='B' then
                    Iformula = 'B'
                elseif varargin_data(Iformula_index+1)=='T' then
                    Iformula = 'T'
                end
                
                 //selection of derivative formula
                if Dformula_index==0 |varargin_data(Dformula_index+1)=='F' then
                    Dformula = 'F'
                elseif varargin_data(Dformula_index+1)=='B' then
                    Dformula = 'B'
                elseif varargin_data(Dformula_index+1)=='T' then
                    Dformula = 'T'
                end                
                                    
                    Iformula_method = Iformula
                    Dformula_method = Dformula
//extract the element from the given trnsfer function
//____________________________________________________________________________________________________________________________________________________________________________
                    for ii = 1:total_num_element
//------------------------------------------------------------------------------
// finding out the transfer function 
                    //trans_func = systf(ii)
                    num_coeff = coeff(systf_num(ii))
                    den_coeff = coeff(systf_den(ii))
                    length_num_coeff = length(num_coeff)
                    length_den_coeff = length(den_coeff)
                    last_num_coeff = num_coeff(length_num_coeff)
                    last_den_coeff = den_coeff(length_den_coeff)
                    num_coeff = num_coeff/last_num_coeff
                    den_coeff = den_coeff/last_den_coeff
//------------------------------------------------------------------------------
                            //type of PID equation 
                            //(Iformula Dformula)
                            //(F F),(F B) (F T)
                            //(B F),(B B) (B T)
                            //(T F),(T B) (T T)
                            if Iformula =='F' & Dformula == 'F' then
                    //******************************************************************************
                    //printf('\n do ***')
                    // when Ifomula = Forward Euler and Dformula = Forward Euler
                            // For P type controller
                                if length_num_coeff == 1 & length_den_coeff == 1 then
                                    temp_Kp_data(ii,1) = last_num_coeff/last_den_coeff
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for I type Controller 
                            //           Ts
                            //      Ki*-----
                            //          z-1
                                elseif length_num_coeff == 1 & length_den_coeff == 2 then
                                        temp_Kp_data(ii,1) = 0
                                        temp_Ki_data(ii,1) = last_num_coeff/sysTs
                                        temp_Kd_data(ii,1) = 0
                                        temp_Tf_data(ii,1) = 0
                            // for D type Controller
                            //               1
                            //       Kd*---------------
                            //                   Ts
                            //          Tf+  --------
                            //                  z-1
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1)== -1 & den_coeff(1)~=1 then
                                    temp_Tf_data(ii,1) = sysTs/(1+den_coeff(1))
                                    temp_Kd_data(ii,1) = temp_Tf_data(ii,1)*(last_num_coeff/last_den_coeff)
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kp_data(ii,1) = 0
                             // for PI type Controller
                             //           Ts
                             //  Kp+Ki*---------
                             //          z-1
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & den_coeff(1)== -1 then
                                    temp_Kp_data(ii,1) = last_num_coeff/last_den_coeff
                                    temp_Ki_data(ii,1) = temp_Kp_data(ii,1)*(1+num_coeff(1))/sysTs
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for PID type Controller
                            //              Ts         1
                            //   Kp+ Ki*-------+Kd*--------------
                            //              z-1              Ts
                            //                      Tf+---------
                            //                            z-1
                                elseif length_num_coeff == 3 &  length_den_coeff == 3  then //& () == (sysTs/(1-den_coeff(1)))
                                    temp_Tf_data(ii,1) = sysTs/(1-den_coeff(1))
                                    comparison_Tf = sysTs/(2+den_coeff(2))
                                    temp_variable = comparison_Tf - temp_Tf_data(ii,1)
                                    if round(temp_variable) == 0 then    
                                        A = [temp_Tf_data(ii,1),0,1;
                                             ((sysTs-2*temp_Tf_data(ii,1))-num_coeff(2)*temp_Tf_data(ii,1)),(sysTs*temp_Tf_data(ii,1)),-(2+num_coeff(2));
                                             ((temp_Tf_data(ii,1)-sysTs)-num_coeff(1)*temp_Tf_data(ii,1)),(sysTs*(sysTs-temp_Tf_data(ii,1))),(1-num_coeff(1))]
                                        b = [last_num_coeff*temp_Tf_data(ii,1)/last_den_coeff 0 0]
                                        temp_pid= A\b'
                                        temp_Kp_data(ii,1) = temp_pid(1)
                                        temp_Ki_data(ii,1) = temp_pid(2)
                                        temp_Kd_data(ii,1) = temp_pid(3)
                                      else
                                          disp(systf_num(ii)/(systf_den(ii)))
                                          error(msprintf(gettext("Above equation is not in parallel PID form")))
                                      end
                                      
                            // for PD type Controller
                            //               1
                            //   Kp+ Kd*--------------
                            //                Ts
                            //           Tf+---------
                            //                z-1
                                elseif length_num_coeff == 2 & length_den_coeff == 2 then
                                    temp_Tf_data(ii,1) = sysTs/(1+den_coeff(1))
                                    A = [temp_Tf_data(ii,1) 1;
                                         (sysTs-temp_Tf_data(ii,1))-num_coeff(1)*temp_Tf_data(ii,1) -(1+num_coeff(1))]
                                    b = [last_num_coeff*temp_Tf_data(ii,1)/last_den_coeff 0]
                                    temp_pid = A\b'
                                    temp_Kp_data(ii,1) = temp_pid(1)
                                    temp_Kd_data(ii,1) = temp_pid(2)
                                    temp_Ki_data(ii,1) = 0
                                else
                                    disp(systf_num(ii)/(systf_den(ii)))
                                    error(msprintf(gettext("Above equation is not in parallel PID form"))) 
                            end
                                
                    // when Ifomula = Forward Euler and Dformula = Backward Euler        
                    //******************************************************************************
                            elseif Iformula == 'F' & Dformula == 'B' then 
                    //******************************************************************************
                                // For P type controller
                                if length_num_coeff == 1 & length_den_coeff == 1 then
                                    temp_Kp_data(ii,1) = last_num_coeff/last_den_coeff
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for I type Controller 
                            //           Ts
                            //      Ki*-----
                            //          z-1
                                elseif length_num_coeff == 1 & length_den_coeff == 2 then
                                    temp_Kp_data(ii,1) = 0
                                    temp_Ki_data(ii,1) = last_num_coeff/(sysTs*last_den_coeff)
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for D type Controller
                            //               1
                            //       Kd*---------------
                            //                 z*Ts
                            //          Tf+  --------
                            //                  z-1
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1)== -1 then
                                    //printf('\n derivative type controller\n')
                                    temp_Tf_data(ii,1) = -sysTs*den_coeff(1)/(1+den_coeff(1))
                                    temp_Kd_data(ii,1) = (sysTs+temp_Tf_data(ii,1))*last_num_coeff/last_den_coeff
                                    temp_Kp_data(ii,1) = 0
                                    temp_Ki_data(ii,1) = 0
                            // for PI type Controller
                            //           Ts
                            //  Kp+Ki*---------
                            //          z-1
                                elseif length_num_coeff ==2 & length_den_coeff ==2 & den_coeff(1)==-1 then
                                    temp_Kp_data(ii,1) = last_num_coeff/last_den_coeff
                                    temp_Ki_data(ii,1) = temp_Kp_data(ii,1)*(1+num_coeff(1))/sysTs
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for PID type Controller
                            //              Ts         1
                            //   Kp+ Ki*-------+Kd*--------------
                            //              z-1           z*Ts
                            //                      Tf+---------
                            //                            z-1
                                elseif length_num_coeff == 3 & length_den_coeff ==3 then
                                    temp_Tf_data(ii,1) = sysTs*den_coeff(1)/(1-den_coeff(1))
                                    comparison_Tf = -sysTs*(1+den_coeff(2))/(2+den_coeff(2))
                                    temp_variable = comparison_Tf - temp_Tf_data(ii,1)
                                    if round(temp_variable) == 0 then 
                                        A = [(temp_Tf_data(ii,1)+sysTs) 0 1;
                                              -((2*temp_Tf_data(ii,1)+sysTs)+(temp_Tf_data(ii,1)+sysTs)*num_coeff(2)) sysTs*(temp_Tf_data(ii,1)+sysTs) -(2+num_coeff(2));
                                               temp_Tf_data(ii,1)-(temp_Tf_data(ii,1)+sysTs)*num_coeff(1) -sysTs*temp_Tf_data(ii,1) 1-num_coeff(1)]
                                        b = [(sysTs+temp_Tf_data(ii,1))*last_num_coeff/last_den_coeff; 0;0]
                                        temp_pid = A\b
                                        temp_Kp_data(ii,1) = temp_pid(1)
                                        temp_Ki_data(ii,1) = temp_pid(2)
                                        temp_Kd_data(ii,1) = temp_pid(3)
                                    else
                                        disp(systf_num(ii)/(systf_den(ii)))
                                        error(msprintf(gettext("Above equation is not in parallel PID form")))
                                    end
                                
                            // for PD type Controller
                            //              1
                            //   Kp+Kd*--------------
                            //               z*Ts
                            //         Tf+---------
                            //               z-1
                            elseif length_num_coeff == 2 & length_den_coeff ==2 then
                                    temp_Tf_data(ii,1) = -sysTs*den_coeff(1)/(1+den_coeff(1))

                                    A = [temp_Tf_data(ii,1)+sysTs 1;
                                        -(temp_Tf_data(ii,1)+(temp_Tf_data(ii,1)+sysTs)*num_coeff(1)) -(1+ num_coeff(1))]
                                    b = [(sysTs+temp_Tf_data(ii,1))*last_num_coeff/last_den_coeff; 0]
                                    temp_pid = A\b
                                    temp_Kp_data(ii,1) = temp_pid(1)
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = temp_pid(2)
                            else
                                    disp(systf_num(ii)/(systf_den(ii)))
                                    error(msprintf(gettext("Above equation is not in parallel PID form")))
                            end
                    // when Ifomula = Forward Euler and Dformula = Trapezoidal
                    //******************************************************************************
                            elseif Iformula == 'F' & Dformula == 'T' then 
                    //******************************************************************************
                                // For P type controller
                                if length_num_coeff == 1 & length_den_coeff == 1 then
                                    temp_Kp_data(ii,1) = last_num_coeff/last_den_coeff
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for I type Controller 
                            //           Ts
                            //      Ki*-----
                            //          z-1
                                elseif length_num_coeff == 1 & length_den_coeff == 2 then
                                    temp_Kp_data(ii,1) = 0
                                    temp_Ki_data(ii,1) = last_num_coeff/sysTs
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for D type Controller
                            //               1
                            //       Kd*---------------
                            //               Ts (z+1)
                            //          Tf+  --------
                            //               2  (z-1)
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1)==-1 then
                                    temp_Tf_data(ii,1) = sysTs*(1-den_coeff(1))/(2*(1+den_coeff(1)))
                                    temp_Kd_data(ii,1) = (2*temp_Tf_data(ii,1)+sysTs)*last_num_coeff/(2*last_den_coeff)
                                    temp_Kp_data(ii,1) = 0
                                    temp_Ki_data(ii,1) = 0
                            // for PI type Controller
                            //           Ts
                            //  Kp+Ki*---------
                            //          z-1
                                elseif length_num_coeff ==2 & length_den_coeff ==2 & den_coeff(1)==-1 then
                                    //
                                    temp_Kp_data(ii,1) = last_num_coeff
                                    temp_Ki_data(ii,1) = temp_Kp_data(ii,1)*(1+num_coeff(1))/sysTs
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for PID type Controller
                            //              Ts         1
                            //   Kp+ Ki*-------+Kd*--------------
                            //              z-1          Ts*(z+1)
                            //                      Tf+---------
                            //                           2*(z-1)
                                elseif length_num_coeff == 3 & length_den_coeff == 3 then
                                    temp_Tf_data(ii,1) = -sysTs*(1+den_coeff(1))/(2*(den_coeff(1)-1))
                                    comparison_Tf = -sysTs*den_coeff(2)/(4+2*den_coeff(2))
                                    temp_variable = comparison_Tf - temp_Tf_data(ii,1)
                                    if round(temp_variable) == 0 then
                                        A =[(2*temp_Tf_data(ii,1)+sysTs) 0 2;
                                            -(4*temp_Tf_data(ii,1)+(2*temp_Tf_data(ii,1)+sysTs)*num_coeff(2)) sysTs*(2*temp_Tf_data(ii,1)+sysTs) -(4+2*num_coeff(2));
                                            -((sysTs-2*temp_Tf_data(ii,1))+(2*temp_Tf_data(ii,1)+sysTs)*num_coeff(1)) sysTs*(sysTs-2*temp_Tf_data(ii,1)) (2-2*num_coeff(1))]
                                        b = [(2*temp_Tf_data(ii,1)+sysTs)*last_num_coeff/last_den_coeff;0;0]
                                        temp_pid = A\b
                                        temp_Kp_data(ii,1) = temp_pid(1)
                                        temp_Ki_data(ii,1) = temp_pid(2)
                                        temp_Kd_data(ii,1) = temp_pid(3)
                                    else
                                        disp(systf_num(ii)/(systf_den(ii)))
                                        error(msprintf(gettext("Above equation is not in parallel PID form")))
                                    end
                            // for PD type Controller
                            //              1
                            //   Kp+Kd*--------------
                            //              Ts*(z+1)
                            //         Tf+---------
                            //               2(z-1)
                                elseif length_num_coeff == 2 & length_den_coeff == 2 then
                                    temp_Tf_data(ii,1) = sysTs*(1-den_coeff(1))/(2*(1+den_coeff(1)))
                                    A = [(2*temp_Tf_data(ii,1)+sysTs) 2;
                                          (sysTs-2*temp_Tf_data(ii,1))-(2*temp_Tf_data(ii,1)+sysTs)*num_coeff(1) -2*(1+num_coeff(1))]
                                    b = [(2*temp_Tf_data(ii,1)+sysTs)*last_num_coeff/last_den_coeff;0]
                                    temp_pid = A\b
                                    temp_Kp_data(ii,1) = temp_pid(1)
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = temp_pid(2)
                                else
                                        disp(systf_num(ii)/(systf_den(ii)))
                                        error(msprintf(gettext("Above equation is not in parallel PID form")))
                                end
                    // when Ifomula = Backward Euler and Dformula = Forward Eular
                    //******************************************************************************
                            elseif Iformula == 'B' & Dformula == 'F' then 
                    //******************************************************************************
                                // For P type controller
                                if length_num_coeff == 1 & length_den_coeff == 1 then
                                    temp_Kp_data(ii,1) = last_num_coeff/last_den_coeff
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for I type Controller 
                            //           Ts*z
                            //      Ki*-------
                            //           z-1
                                elseif length_num_coeff ==2&length_den_coeff==2&num_coeff(1)==0&den_coeff(1)==-1 then
                                    temp_Ki_data(ii,1) = last_num_coeff/(sysTs*last_den_coeff)
                                    temp_Kp_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for D type Controller
                            //               1
                            //       Kd*---------------
                            //                   Ts
                            //          Tf+  --------
                            //                  z-1
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1)== -1 & den_coeff(1)~=1 then
                                    temp_Tf_data(ii,1) = sysTs/(1+den_coeff(1))
                                    temp_Kd_data(ii,1) = temp_Tf_data(ii,1)*(last_num_coeff/last_den_coeff)
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kp_data(ii,1) = 0
                            // for PI type Controller
                            //          Ts*z
                            //  Kp+Ki*---------
                            //          z-1
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & den_coeff(1)== -1 then
                                    A = [1 sysTs;
                                         1+num_coeff(1) sysTs*num_coeff(1)]
                                    b = [last_num_coeff/last_den_coeff; 0]
                                    temp_pid = A\b
                                    temp_Kp_data(ii,1) = temp_pid(1)
                                    temp_Ki_data(ii,1) = temp_pid(2)
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for PID type Controller
                            //              z*Ts         1
                            //   Kp+ Ki*-------+Kd*--------------
                            //              z-1           Ts
                            //                      Tf+---------
                            //                            z-1
                                elseif length_num_coeff == 3 & length_den_coeff ==3 then
                                    temp_Tf_data(ii,1) = sysTs/(1-den_coeff(1))
                                    comparison_Tf = sysTs/(2+den_coeff(2))
                                    temp_variable = comparison_Tf - temp_Tf_data(ii,1)
                                    if round(temp_variable) == 0 then
                                        A = [temp_Tf_data(ii,1) temp_Tf_data(ii,1)*sysTs 1;
                                             ((sysTs-2*temp_Tf_data(ii,1))-(temp_Tf_data(ii,1)*num_coeff(2))) (sysTs*(sysTs-temp_Tf_data(ii,1))-temp_Tf_data(ii,1)*sysTs*num_coeff(2)) -(2+num_coeff(2));
                                             ((temp_Tf_data(ii,1)-sysTs)-(temp_Tf_data(ii,1)*num_coeff(1))) -temp_Tf_data(ii,1)*sysTs*num_coeff(1) 1-num_coeff(1)]
                                        b = [temp_Tf_data(ii,1)*last_num_coeff/last_den_coeff;0;0]
                                        temp_pid = A\b
                                        temp_Kp_data(ii,1) = temp_pid(1)
                                        temp_Ki_data(ii,1) = temp_pid(2)
                                        temp_Kd_data(ii,1) = temp_pid(3)
                                     else
                                        disp(systf_num(ii)/(systf_den(ii)))
                                        error(msprintf(gettext("Above equation is not in parallel PID form")))
                                    end
                                    
                            // for PD type Controller
                            //               1
                            //   Kp+ Kd*--------------
                            //                Ts
                            //           Tf+---------
                            //                z-1
                                elseif length_num_coeff == 2 & length_den_coeff == 2 then
                                    temp_Tf_data(ii,1) = sysTs/(1+den_coeff(1))
                                    A = [temp_Tf_data(ii,1) 1;
                                         (sysTs-temp_Tf_data(ii,1))-num_coeff(1)*temp_Tf_data(ii,1) -(1+num_coeff(1))]
                                    b = [last_num_coeff*temp_Tf_data(ii,1) 0]
                                    temp_pid = A\b'
                                    temp_Kp_data(ii,1) = temp_pid(1)
                                    temp_Kd_data(ii,1) = temp_pid(2)
                                    temp_Ki_data(ii,1) = 0 
                                else
                                    disp(systf_num(ii)/(systf_den(ii)))
                                    error(msprintf(gettext("Above equation is not in parallel PID form")))
                                end
                    // when Ifomula = Backward Euler and Dformula = Backward Eular
                    //******************************************************************************
                            elseif Iformula == 'B' & Dformula == 'B' then 
                    //******************************************************************************
                                // For P type controller
                                if length_num_coeff == 1 & length_den_coeff == 1 then
                                    temp_Kp_data(ii,1) = last_num_coeff/last_den_coeff
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for I type Controller 
                            //           Ts*z
                            //      Ki*-------
                            //           z-1
                                elseif length_num_coeff ==2 & length_den_coeff==2&num_coeff(1)==0 & den_coeff(1)==-1 then
                                    temp_Ki_data(ii,1) = last_num_coeff/(sysTs*last_den_coeff)
                                    temp_Kp_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for D type Controller
                            //               1
                            //       Kd*---------------
                            //                 z*Ts
                            //          Tf+  --------
                            //                  z-1
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1)== -1 then
                                    //printf('\n derivative type controller\n')
                                    temp_Tf_data(ii,1) = -sysTs*den_coeff(1)/(1+den_coeff(1))
                                    temp_Kd_data(ii,1) = (sysTs+temp_Tf_data(ii,1))*last_num_coeff/last_den_coeff
                                    temp_Kp_data(ii,1) = 0
                                    temp_Ki_data(ii,1) = 0
                            // for PI type Controller
                            //          Ts*z
                            //  Kp+Ki*---------
                            //          z-1
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & den_coeff(1)== -1 then
                                    A = [1 sysTs;
                                         1+num_coeff(1) sysTs*num_coeff(1)]
                                    b = [last_num_coeff/last_den_coeff; 0]
                                    temp_pid = A\b
                                    temp_Kp_data(ii,1) = temp_pid(1)
                                    temp_Ki_data(ii,1) = temp_pid(2)
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for PID type Controller
                            //              z*Ts         1
                            //   Kp+ Ki*-------+Kd*--------------
                            //              z-1           z*Ts
                            //                      Tf+---------
                            //                            z-1
                                elseif length_num_coeff == 3 & length_den_coeff ==3 then
                                    temp_Tf_data(ii,1) = sysTs*den_coeff(1)/(1-den_coeff(1))
                                    comparison_Tf = sysTs*(den_coeff(2)-1)/(2-den_coeff(2))
                                    temp_variable = comparison_Tf - temp_Tf_data(ii,1)
                                    if round(temp_variable) == 0 then
                                        A =[(temp_Tf_data(ii,1)+sysTs) sysTs*(temp_Tf_data(ii,1)+sysTs) 1;
                                            ((2*temp_Tf_data(ii,1)+sysTs)+(temp_Tf_data(ii,1)+sysTs)*num_coeff(2)) (temp_Tf_data(ii,1)*sysTs+sysTs*(sysTs+temp_Tf_data(ii,1))*num_coeff(2)) (2+num_coeff(2));
                                            (temp_Tf_data(ii,1)-(temp_Tf_data(ii,1)+sysTs)*num_coeff(1)) -sysTs*(temp_Tf_data(ii,1)+sysTs)*num_coeff(1) (1-num_coeff(1))]
                                        b = [(temp_Tf_data(ii,1)+sysTs)*last_num_coeff/last_den_coeff;0;0]
                                        temp_pid = A\b
                                        temp_Kp_data(ii,1) = temp_pid(1)
                                        temp_Ki_data(ii,1) = temp_pid(2)
                                        temp_Kd_data(ii,1) = temp_pid(3)
                                    else
                                        disp(systf_num(ii)/(systf_den(ii)))
                                        error(msprintf(gettext("Above equation is not in parallel PID form")))
                                    end
                            // for PD type Controller
                            //              1
                            //   Kp+Kd*--------------
                            //               z*Ts
                            //         Tf+---------
                            //               z-1
                                elseif length_num_coeff == 2 & length_den_coeff ==2 then
                                    temp_Tf_data(ii,1) = -sysTs*den_coeff(1)/(1+den_coeff(1))
                                    A = [temp_Tf_data(ii,1)+sysTs 1;
                                        -(temp_Tf_data(ii,1)+(temp_Tf_data(ii,1)+sysTs)*num_coeff(1)) -(1+ num_coeff(1))]
                                    b = [(sysTs+temp_Tf)*last_num_coeff/last_den_coeff; 0]
                                    temp_pid = A\b
                                    temp_Kp_data(ii,1) = temp_pid(1)
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = temp_pid(2)
                                else
                                    disp(systf_num(ii)/(systf_den(ii)))
                                    error(msprintf(gettext("Above equation is not in parallel PID form")))
                                end
                    // when Ifomula = Backward Euler and Dformula = Trapezoidal Eular
                    //******************************************************************************
                            elseif Iformula == 'B' & Dformula == 'T' then 
                    //******************************************************************************
                                // For P type controller
                                if length_num_coeff == 1 & length_den_coeff == 1 then
                                    temp_Kp_data(ii,1) = last_num_coeff/last_den_coeff
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for I type Controller 
                            //           Ts*z
                            //      Ki*-------
                            //           z-1
                                elseif length_num_coeff ==2 & length_den_coeff==2 & num_coeff(1)==0 & den_coeff(1)==-1 then
                                    temp_Ki_data(ii,1) = last_num_coeff/(sysTs*last_den_coeff)
                                    temp_Kp_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for D type Controller
                            //               1
                            //       Kd*---------------
                            //               Ts (z+1)
                            //          Tf+  --------
                            //               2  (z-1)
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1)==-1 then
                                    temp_Tf_data(ii,1) = sysTs*(1-den_coeff(1))/(2*(1+den_coeff(1)))
                                    temp_Kd_data(ii,1) = (2*temp_Tf_data(ii,1)+sysTs)*last_num_coeff/(2*last_den_coeff)
                                    temp_Kp_data(ii,1) = 0
                                    temp_Ki_data(ii,1) = 0
                            // for PI type Controller
                            //          Ts*z
                            //  Kp+Ki*---------
                            //          z-1
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & den_coeff(1)== -1 then
                                    A = [1 sysTs;
                                         1+num_coeff(1) sysTs*num_coeff(1)]
                                    b = [last_num_coeff/last_den_coeff; 0]
                                    temp_pid = A\b
                                    temp_Kp_data(ii,1) = temp_pid(1)
                                    temp_Ki_data(ii,1) = temp_pid(2)
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for PID type Controller
                            //              z*Ts         1
                            //   Kp+ Ki*-------+Kd*--------------
                            //              z-1          Ts*(z+1)
                            //                      Tf+---------
                            //                           2*(z-1)
                                elseif length_num_coeff == 3 & length_den_coeff == 3 then
                                    temp_Tf_data(ii,1) = sysTs*(1+den_coeff(1))/(2*(1-den_coeff(1)))
                                    comparison_Tf = -sysTs*den_coeff(2)/(4+2*den_coeff(2))
                                    temp_variable = comparison_Tf - temp_Tf_data(ii,1)
                                    if round(temp_variable) == 0 then
                                        A = [(2*temp_Tf_data(ii,1)+sysTs) sysTs*(sysTs+2*temp_Tf_data(ii,1)) 2;
                                             -(4*temp_Tf_data(ii,1)+(2*temp_Tf_data(ii,1)+sysTs)*num_coeff(2)) sysTs*(sysTs-2*temp_Tf_data(ii,1))-sysTs*(sysTs+2*temp_Tf_data(ii,1))*num_coeff(2) -(4+2*num_coeff(2));
                                             ((2*temp_Tf_data(ii,1)-sysTs)-(2*temp_Tf_data(ii,1)+sysTs)*num_coeff(1)) -sysTs*(sysTs+2*temp_Tf_data(ii,1))*num_coeff(1) (2-2*num_coeff(1))]
                                        b = [(2*temp_Tf_data(ii,1)+sysTs)*last_num_coeff/last_den_coeff;0;0]
                                        temp_pid = A\b
                                        temp_Kp_data(ii,1) = temp_pid(1)
                                        temp_Ki_data(ii,1) = temp_pid(2)
                                        temp_Kd_data(ii,1) = temp_pid(3)
                                    else
                                    disp(systf_num(ii)/(systf_den(ii)))
                                    error(msprintf(gettext("Above equation is not in parallel PID form")))
                                    end
                            // for PD type Controller
                            //              1
                            //   Kp+Kd*--------------
                            //              Ts*(z+1)
                            //         Tf+---------
                            //               2(z-1)
                                elseif length_num_coeff == 2 & length_den_coeff == 2 then
                                    temp_Tf_data(ii,1) = sysTs*(1-den_coeff(1))/(2*(1+den_coeff(1)))
                                    A = [(2*temp_Tf_data(ii,1)+sysTs) 2;
                                          (sysTs-2*temp_Tf_data(ii,1))-(2*temp_Tf_data(ii,1)+sysTs)*num_coeff(1) -2*(1+num_coeff(1))]
                                    b = [(2*temp_Tf+sysTs)*last_num_coeff/last_den_coeff;0]
                                    temp_pid = A\b
                                    temp_Kp_data(ii,1) = temp_pid(1)
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = temp_pid(2)
                                else
                                    disp(systf_num(ii)/(systf_den(ii)))
                                    error(msprintf(gettext("Above equation is not in parallel PID form")))
                                end
                    // when Ifomula = Trapezoidal Euler and Dformula = Forward Euler
                    //******************************************************************************
                            elseif Iformula == 'T' & Dformula == 'F' then 
                    //******************************************************************************
                                // For P type controller
                                if length_num_coeff == 1 & length_den_coeff == 1 then
                                    temp_Kp_data(ii,1) = last_num_coeff/last_den_coeff
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for I type Controller 
                            //           Ts*(z+1)
                            //      Ki*------------
                            //            2*(z-1)            
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1) == 1 & den_coeff(1) == -1 then
                                    temp_Ki_data(ii,1) = 2*last_num_coeff/(sysTs*last_den_coeff)
                                    temp_Kp_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for D type Controller
                            //               1
                            //       Kd*---------------
                            //                   Ts
                            //          Tf+  --------
                            //                  z-1
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1)== -1 & den_coeff(1)~=1 then
                                    temp_Tf_data(ii,1) = sysTs/(1+den_coeff(1))
                                    temp_Kd_data(ii,1) = temp_Tf_data(ii,1)*(last_num_coeff/last_den_coeff)
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kp_data(ii,1) = 0
                            // for PI type Controller
                            //          Ts*(z+1)
                            //  Kp+Ki*---------
                            //          2*(z-1)
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & den_coeff(1)==-1 then
                                    A = [2 sysTs;
                                         -2*(1+num_coeff(1)) sysTs*(1-num_coeff(1))]
                                    b = [2*last_num_coeff/last_den_coeff;0]
                                    temp_pid = A\b
                                    temp_Kp_data(ii,1) = temp_pid(1)
                                    temp_Ki_data(ii,1) = temp_pid(2)
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for PID type Controller
                            //             Ts*(z+1)         1

                            //   Kp+ Ki*-------------+Kd*--------------
                            //              2*(z-1)             Ts
                            //                           Tf+---------
                            //                                 z-1
                                elseif length_num_coeff == 3 & length_den_coeff == 3 then
                                    temp_Tf_data(ii,1) = sysTs/(1-den_coeff(1))
                                    comparison_Tf = sysTs/(2+den_coeff(2))
                                    temp_variable = comparison_Tf - temp_Tf_data(ii,1)
                                    if round(temp_variable) == 0 then
                                        A = [2*temp_Tf_data(ii,1) temp_Tf_data(ii,1)*sysTs 2;
                                            ((2*sysTs-4*temp_Tf_data(ii,1))-2*temp_Tf_data(ii,1)*num_coeff(2)) (sysTs^2-temp_Tf_data(ii,1)*sysTs*num_coeff(2)) -(4+2*num_coeff(2));
                                            (2*(temp_Tf_data(ii,1)-sysTs)-2*temp_Tf_data(ii,1)*num_coeff(1)) (sysTs*(sysTs-temp_Tf_data(ii,1))-temp_Tf_data(ii,1)*sysTs*num_coeff(1)) (2-2*num_coeff(1))]
                                        b = [2*temp_Tf_data(ii,1)*last_num_coeff/last_den_coeff;0;0]
                                        temp_pid = A\b
                                        temp_Kp_data(ii,1) = temp_pid(1)
                                        temp_Ki_data(ii,1) = temp_pid(2)
                                        temp_Kd_data(ii,1) = temp_pid(3)
                                    else
                                        disp(systf_num(ii)/(systf_den(ii)))
                                        error(msprintf(gettext("Above equation is not in parallel PID form")))
                                    end
                            // for PD type Controller
                            //               1
                            //   Kp+ Kd*--------------
                            //                Ts
                            //           Tf+---------
                            //                z-1
                                elseif length_num_coeff == 2 & length_den_coeff == 2 then
                                    temp_Tf_data(ii,1) = sysTs/(1+den_coeff(1))
                                    A = [temp_Tf_data(ii,1) 1;
                                         (sysTs-temp_Tf_data(ii,1))-num_coeff(1)*temp_Tf_data(ii,1) -(1+num_coeff(1))]
                                    b = [last_num_coeff*temp_Tf_data(ii,1) 0]
                                    temp_pid = A\b'
                                    temp_Kp_data(ii,1) = temp_pid(1)
                                    temp_Kd_data(ii,1) = temp_pid(2)
                                    temp_Ki_data(ii,1) = 0
                                else
                                    disp(systf_num(ii)/(systf_den(ii)))
                                    error(msprintf(gettext("Above equation is not in parallel PID form")))
                                end
                    // when Ifomula = Trapezoidal Euler and Dformula = Backward Euler
                    //******************************************************************************
                            elseif Iformula == 'T' & Dformula == 'B' then 
                    //******************************************************************************
                                // For P type controller
                                if length_num_coeff == 1 & length_den_coeff == 1 then
                                    temp_Kp_data(ii,1) = last_num_coeff/last_den_coeff
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for I type Controller 
                            //           Ts*(z+1)
                            //      Ki*------------
                            //            2*(z-1)            
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1) == 1 & den_coeff(1) == -1 then
                                    temp_Ki_data(ii,1) = 2*last_num_coeff/(sysTs*last_den_coeff)
                                    temp_Kp_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for D type Controller
                            //               1
                            //       Kd*---------------
                            //                 z*Ts
                            //          Tf+  --------
                            //                  z-1
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1)== -1 then
                                    //printf('\n derivative type controller\n')
                                    temp_Tf_data(ii,1) = -sysTs*den_coeff(1)/(1+den_coeff(1))
                                    temp_Kd_data(ii,1) = (sysTs+temp_Tf_data(ii,1))*last_num_coeff/last_den_coeff
                                    temp_Kp_data(ii,1) = 0
                                    temp_Ki_data(ii,1) = 0
                            // for PI type Controller
                            //          Ts*(z+1)
                            //  Kp+Ki*---------
                            //          2*(z-1)
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & den_coeff(1)==-1 then
                                    A = [2 sysTs;
                                         -2*(1+num_coeff(1)) sysTs*(1-num_coeff(1))]
                                    b = [2*last_num_coeff/last_den_coeff;0]
                                    temp_pid = A\b
                                    temp_Kp_data(ii,1) = temp_pid(1)
                                    temp_Ki_data(ii,1) = temp_pid(2)
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for PID type Controller
                            //             Ts*(z+1)         1
                            //   Kp+ Ki*-------------+Kd*--------------
                            //              2*(z-1)             Ts
                            //                           Tf+---------
                            //                                 z-1
                                elseif length_num_coeff == 3 & length_den_coeff == 3 then
                                    temp_Tf_data(ii,1) = sysTs*den_coeff(1)/(1-den_coeff(1))
                                    comparison_Tf = -sysTs*(den_coeff(2)+1)/(2+den_coeff(2))
                                    temp_variable = comparison_Tf - temp_Tf_data(ii,1)
                                    if round(temp_variable) == 0 then
                                        A = [2*(temp_Tf_data(ii,1)+sysTs) sysTs*(temp_Tf_data(ii,1)+sysTs) 2;
                                             -((2*sysTs+4*temp_Tf_data(ii,1))+2*(temp_Tf_data(ii,1)+sysTs)*num_coeff(2)) (sysTs^2 -sysTs*(sysTs+temp_Tf_data(ii,1))*num_coeff(2)) -(4+2*num_coeff(2));
                                             (2*temp_Tf_data(ii,1)-2*(temp_Tf_data(ii,1)+sysTs)*num_coeff(1)) -(temp_Tf_data(ii,1)*sysTs+sysTs*(sysTs+temp_Tf_data(ii,1))*num_coeff(1)) 2*(1-num_coeff(1))]
                                        b = [2*(temp_Tf_data(ii,1)+sysTs)*last_num_coeff/last_den_coeff;0;0]
                                        temp_pid = A\b
                                        temp_Kp_data(ii,1) = temp_pid(1)
                                        temp_Ki_data(ii,1) = temp_pid(2)
                                        temp_Kd_data(ii,1) = temp_pid(3)
                                    else
                                        disp(systf_num(ii)/(systf_den(ii)))
                                        error(msprintf(gettext("Above equation is not in parallel PID form")))
                                    end
                            // for PD type Controller
                            //              1
                            //   Kp+Kd*--------------
                            //               z*Ts
                            //         Tf+---------
                            //               z-1
                                elseif length_num_coeff == 2 & length_den_coeff ==2 then
                                    temp_Tf_data(ii,1) = -sysTs*den_coeff(1)/(1+den_coeff(1))
                                    A = [temp_Tf_data(ii,1)+sysTs 1;
                                        -(temp_Tf_data(ii,1)+(temp_Tf_data(ii,1)+sysTs)*num_coeff(1)) -(1+ num_coeff(1))]
                                    b = [(sysTs+temp_Tf_data(ii,1))*last_num_coeff/last_den_coeff; 0]
                                    temp_pid = A\b
                                    temp_Kp_data(ii,1) = temp_pid(1)
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = temp_pid(2)                
                                else
                                    disp(systf_num(ii)/(systf_den(ii)))
                                    error(msprintf(gettext("Above equation is not in parallel PID form")))
                                end
                    // when Ifomula = Trapezoidal Euler and Dformula = Trapezoidal
                    //******************************************************************************
                            elseif Iformula == 'T' & Dformula == 'T' then 
                    //******************************************************************************
                                // For P type controller
                                if length_num_coeff == 1 & length_den_coeff == 1 then
                                    temp_Kp_data(ii,1) = last_num_coeff/last_den_coeff
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for I type Controller 
                            //           Ts*(z+1)
                            //      Ki*------------
                            //            2*(z-1)            
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1) == 1 & den_coeff(1) == -1 then
                                    temp_Ki_data(ii,1) = 2*last_num_coeff/(sysTs*last_den_coeff)
                                    temp_Kp_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for D type Controller
                            //               1
                            //       Kd*---------------
                            //               Ts (z+1)
                            //          Tf+  --------
                            //               2  (z-1)
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1)==-1 then
                                    temp_Tf_data(ii,1) = sysTs*(1-den_coeff(1))/(2*(1+den_coeff(1)))
                                    temp_Kd_data(ii,1) = (2*temp_Tf_data(ii,1)+sysTs)*last_num_coeff/(2*last_den_coeff)
                                    temp_Kp_data(ii,1) = 0
                                    temp_Ki_data(ii,1) = 0
                             // for PI type Controller
                            //          Ts*(z+1)
                            //  Kp+Ki*---------
                            //          2*(z-1)
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & den_coeff(1)==-1 then
                                    A = [2 sysTs;
                                         -2*(1+num_coeff(1)) sysTs*(1-num_coeff(1))]
                                    b = [2*last_num_coeff/last_den_coeff;0]
                                    temp_pid = A\b
                                    temp_Kp_data(ii,1) = temp_pid(1)
                                    temp_Ki_data(ii,1) = temp_pid(2)
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                            // for PID type Controller
                            //             Ts*(z+1)         1
                            //   Kp+ Ki*-------------+Kd*--------------
                            //              2*(z-1)          Ts*(z+1)
                            //                           Tf+---------
                            //                                2*(z-1)
                                elseif length_num_coeff == 3 & length_den_coeff == 3 then
                                    temp_Tf_data(ii,1) = sysTs*(1+den_coeff(1))/(2*(1-den_coeff(1)))
                                    comparison_Tf = -sysTs*den_coeff(2)/(4+2*den_coeff(2))
                                    temp_variable = comparison_Tf - temp_Tf_data(ii,1)
                                    if round(temp_variable) == 0 then
                                        A = [(4*temp_Tf_data(ii,1)+2*sysTs) sysTs*(2*temp_Tf_data(ii,1)+sysTs) 4;
                                             -(8*temp_Tf_data(ii,1)+(4*temp_Tf_data(ii,1)+2*sysTs)*num_coeff(2)) (2*sysTs^2-sysTs*(sysTs+2*temp_Tf_data(ii,1))*num_coeff(2)) -(8+4*num_coeff(2));
                                             ((4*temp_Tf_data(ii,1)-2*sysTs)-(4*temp_Tf_data(ii,1)+2*sysTs)*num_coeff(1)) (sysTs*(sysTs-2*temp_Tf_data(ii,1))-sysTs*(2*temp_Tf_data(ii,1)+sysTs)*num_coeff(1)) 4-4*num_coeff(1)]
                                        b = [(4*temp_Tf_data(ii,1)+2*sysTs)*last_num_coeff/last_den_coeff;0;0]
                                        temp_pid = A\b
                                        temp_Kp_data(ii,1) = temp_pid(1)
                                        temp_Ki_data(ii,1) = temp_pid(2)
                                        temp_Kd_data(ii,1) = temp_pid(3)
                                    else
                                        disp(systf_num(ii)/(systf_den(ii)))
                                        error(msprintf(gettext("Above equation is not in parallel PID form")))
                                    end
                                    
                            // for PD type Controller
                            //              1
                            //   Kp+Kd*--------------
                            //              Ts*(z+1)
                            //         Tf+---------
                            //               2(z-1)
                                elseif length_num_coeff == 2 & length_den_coeff == 2 then
                                    temp_Tf_data(ii,1) = sysTs*(1-den_coeff(1))/(2*(1+den_coeff(1)))
                                    A = [(2*temp_Tf_data(ii,1)+sysTs) 2;
                                          (sysTs-2*temp_Tf_data(ii,1))-(2*temp_Tf_data(ii,1)+sysTs)*num_coeff(1) -2*(1+num_coeff(1))]
                                    b = [(2*temp_Tf_data(ii,1)+sysTs)*last_num_coeff/last_den_coeff;0]
                                    temp_pid = A\b
                                    temp_Kp_data(ii,1) = temp_pid(1)
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = temp_pid(2)
                                else
                                    disp(systf_num(ii)/(systf_den(ii)))
                                    error(msprintf(gettext("Above equation is not in parallel PID form")))
                                end
                                
                    //******************************************************************************
                            end
                    //------------------------------------------------------------------------------
                    // when Iformula is Backward euler
                     end// extract the PID data from given transfer function
                    
//______________________________________________________________________________                    
                end// identifty that the given function is "c" type or "d" type
                output = varargin_data(1)
        end
//______________________________________________________________________________

                if current_length<=2 then
                    //output_data = hypermat([current_size], pid_function(:,1))
                    Kp = hypermat([current_size],temp_Kp_data(:,1))
                    Ki = hypermat([current_size],temp_Ki_data(:,1))
                    Kd = hypermat([current_size],temp_Kd_data(:,1))
                    Tf = hypermat([current_size],temp_Tf_data(:,1))
                elseif current_length > 2 then
                    //output_data = hypermat([current_size])
                    //output_data(:,:,:,:)= pid_function
                    Kp = hypermat([current_size])
                    Kp(:,:,:,:) = temp_Kp_data
                    Ki = hypermat([current_size])
                    Ki(:,:,:,:) = temp_Ki_data
                    Kd = hypermat([current_size])
                    Kd(:,:,:,:) = temp_Kd_data
                    Tf = hypermat([current_size])
                    Tf(:,:,:,:) = temp_Tf_data
                end
                impdata.Kp = Kp
                impdata.Ki = Ki
                impdata.Kd = Kd
                impdata.Tf = Tf
                impdata.Ts = sysTs
//------------------------------------------------------------------------------
                if sysTs == 'c' then
                    impdata.Iformula = '[]'
                else
                    if Iformula_method == 'F' then
                        impdata.Iformula = 'ForwardEuler'
                    elseif Iformula_method == 'B' then
                        impdata.Iformula = 'BackwardEuler'
                    elseif Iformula_method == 'T' then
                        impdata.Iformula = 'Trapezoidal'
                    end
                end
//------------------------------------------------------------------------------
                if sysTs == 'c' then
                    impdata.Dformula = '[]'
                else
                    if Dformula_method == 'F' then
                        impdata.Dformula = 'ForwardEuler'
                    elseif Dformula_method == 'B' then
                        impdata.Dformula = 'BackwardEuler'
                    elseif Dformula_method == 'T' then
                        impdata.Dformula = 'Trapezoidal'
                    end
                end
//------------------------------------------------------------------------------
                impdata.InputDelay = 0
                impdata.OutputDelay = 0
//------------------------------------------------------------------------------
                Namedata = 0
                Notesdata = 0
                Datauser = 0
                //Namedata = find( varargin == "Name" )
                //disp(Namedata)
                if count_numb < rhs then
                    for ii = count_numb:rhs
                        if typeof(varargin(ii)) == 'hypermat' | typeof(varargin(ii)) == 'rational' then
                        else
                            if varargin(ii) == 'Name' then
                                Namedata = ii
                            elseif varargin(ii) == 'Notes' then
                                Notesdata = ii
                            elseif varargin(ii) == 'UserData' then
                                Datauser = ii
                            end
                        end
                        
                    end
                end
//------------------------------------------------------------------------------                
                if Namedata == 0 then
                    impdata.Name = []
                elseif Namedata + 1 > rhs then
                    impdata.Name = []
                else
                    impdata.Name = varargin(Namedata+1)
                end
//------------------------------------------------------------------------------                
                if Notesdata == 0 then
                    impdata.Notes = []
                elseif Notesdata + 1 > rhs then
                    impdata.Notes = []
                else
                    impdata.Notes = varargin(Notesdata+1)
                end
//------------------------------------------------------------------------------                
                if Datauser == 0 then
                    impdata.UserData = []
                elseif Datauser + 1 > rhs then
                    impdata.UserData = []
                else
                    impdata.UserData = varargin(datauser)
                end
//------------------------------------------------------------------------------                
//                disp(temp_Kp_data)
//                disp(temp_Ki_data)
//                disp(temp_Kd_data)
//                disp(temp_Tf_data)
                
                                
        
endfunction
