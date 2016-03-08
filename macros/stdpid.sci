

function [output, extdata ] =  stdpid(varargin)
//Generate standard form of PID controller, verify the standard form of PID controller
//Calling Seqence
//output = pid(Kp)
//output = pid(Kp,Ti)
//output = pid(Kp,Ti,Td)
//output = pid(Kp,Ti,Td,N)
//output = pid(Kp,Ti,Td,N,Ts)
//output = pid(sys)
//[output impdata] = pid(Kp)
//[output impdata] = pid(.....,'Notes','...')
//Parameters
//Kp : a matrix of double, represents the  Proportional gain of PID controller. Default value is zero.
//Ti : a matrix of double, represents real and positive Integral time of PID controller. Default value is infinity
//Td : a matrix of double, represents real and positive Derivative time of PID controller. Default value is zero
//N : a matrix of double, represents real and positive devisor of derivative filter. Default value is infinity
//Ts: a double, represents the Sampling timing
//sys: represents the valid form of PID controller   
//Name : a string represent  name of  Control system
//UserData : a string represent UserData of control system
//Notes : a string represent Notes about the control system
//
// Description
//output = pid(Kp,Ti,Td,N) generates the standard form of PID in continuous domain. The mathematical equation is given below
//
//   <latex>
	//    \begin{eqnarray}
    //      \[
    //          output=K_{p}\left(1+\frac{1}{T_{i}}\frac{1}{s}+\frac{T_{d}s}{\frac{T_{d}}{N}s+1}\right)
    //      \]
	//    \end{eqnarray}
	//   </latex>
//
//All the PID parameters (Kp, Ti,Td & N) must be real number and output is the transfer 	function obtained by the given PID parameters data. 
//
//output = pid(Kp,Ti,Td,N,Ts)  generates the standard form of PID in discrete domain and Ts is the sampling timing of the controller. The mathematical equation is given below
//   <latex>
	//    \begin{eqnarray}
    //      \[
    //          output=K_{p}\left(1+\frac{1}{T_{i}}Iformula+\frac{T_{d}s}{\frac{T_{d}}{N}+Dformula}\right)
    //      \]
	//    \end{eqnarray}
	//   </latex>
//
//There are numbers of way to select the sampling formula for discrete time integration and derivation (Iformula,dformula respectively). By default value of the Iformula and Dformula is forward Euler. Other sampling formula is Backward Euler and Trapezoidal. If Dformula is forward Euler and first order derivative filter time constant Tf~=0 then it must satisfy the Td/N >Ts/2 for proper derivative filter pole. The following way to select the sample method
//
//F - Forward Euler
//
//B - Backward Euler
//
//T - Trapezoidal   
//
//output = stdpid(Kp) generates proportional controller in continuous domain
//
//output = stdpid(Kp,Ti) generates proportional+integral controller in continuous domain
//
//output = stdpid(Kp,%inf,Td) generates proportional + derivative controller without first order filter in continuous domain
//
//output = stdpid(Kp,%inf,Td,N) generates proportional + derivative controller with first order filter in continuous domain
//
//output = stdpid(Kp,Ti,Td,%inf) generates proportional +integral+ derivative controller without first order filter in continuous domain
//
//output = stdpid(Kp,Ti,Td,N) generates proportional +integral+ derivative controller with first order filter in continuous domain
//
//output = stdpid(Kp,%inf,0,%inf,Ts) generates proportional controller in discrete domain
//
//output = stdpid(Kp,Ki,0,%inf,Ts) generates proportiona+Integral  controller in discrete domain
//
//output = stdpid(Kp,%inf,Td,%inf,Ts) generates proportional+derivative controller without first order filter in discrete domain
//
//output = stdpid(Kp,%inf,Td,N,Ts) generates proportional+derivative controller with first order filter in discrete domain
//
//output = stdpid(Kp,Ti,Td,%inf,Ts) generates proportional+integral+derivative controller without first order filter in discrete domain
//
//output = stdpid(Kp,Ti,Td,N,Ts) generates proportional+integral+derivative controller with first order filter in discrete domain
//
//output = stdpid(.......,Ts,'Iformula','__','Dformula','__') generates discrete domain controller with selection of sampling formula for Iformula and Dformula  
//
//output  = stdpid(sys) checks the given transfer function and verify it
//
//[output impdata] =stdpid(.....,'Notes','...','UserData','.....','Notes','....') it adds extra descriptive information about the controller in impdata 
// Examples
// output = stdpid(2)
// output = stdpid(2,3)
// output = stdpid(2,3,4)
// output = stdpid(2,3,4,5)
// output = stdpid(2,%inf,4)
// output = stdpid(2,%inf,4,5)
// output = stdpid(2,%inf,0,%inf,0.1)
// output = stdpid(2,3,0,%inf,0.1,'Iformula','B')
// output = stdpid(2,3,4,%inf,0.1,'Iformula','B','Dformula','T')
// output = stdpid(2,3,4,5,0.1,'Iformula','T','Dformula','B')
// s = poly(0,'s'); sys = syslin('c',3*(s+1)*(s+2)/s);
// [output impdata] = stdpid(sys)
//z = poly(0,'z'); sys = syslin(0.1,(19.146667 - 38.793333*z + 19.66*z^2)/(  1.6 - 3.4*z + 1.8*z^2))
// [output impdata]= stdpid(sys,'Iformula','T','Dformula','B')
// Authors
//
// Ashutosh Kumar Bhargava
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
        else
            error(msprintf(gettext("Please enter proper PID data")))
        end
        impdata = tlist(["listtype","Kp","Ki","Kd","N","Iformula","Dformula","InputDelay","OutputDelay","Ts","Name","Notes","UserData"]) //
//____________________________________________________________________________________________________________________________________________________________________________
        // find the count
        count_numb = 0
        tempCount = min(rhs,5)
            for jj = 1:tempCount
                if size(varargin(jj)) == [0 0] then
                    error(msprintf(gettext("Please enter positive real data")))
                elseif typeof(varargin(jj)) == 'constant' | typeof(varargin(jj)) == 'rational' | typeof(varargin(jj)) == 'hypermat' then //then//,, ())  
                    count_numb = count_numb + 1
                else
                    break
                end    
            end
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
            printf('\n Continuous Domain \n')
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
                        if jj == 4 | jj == 2 then 
                            pid_data(jj,1)=%inf
                        else
                            pid_data(jj,1)= 0
                        end
                    end
                //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
                // checking the imaginary and negative values
                if phasemag(pid_data(1,1),'c') ~= 0 then//& pid_data(1,1) < 0 then
                     error(msprintf(gettext("Proportional gain (Kp) must be a finite positive real number.")))
                elseif  pid_data(2,1) == 0 | phasemag(pid_data(2,1)) ~= 0 then
                    error(msprintf(gettext("Integral time (Ti) must be a finite positive real number.")))
                elseif phasemag(pid_data(3,1)) ~= 0 then
                    error(msprintf(gettext("Derivative time (Ti) must be a finite positive real number.")))
                elseif pid_data(4,1) == 0 | phasemag(pid_data(4,1)) ~= 0 then
                    error(msprintf(gettext("Derivative filter divisor(N) must be a finite positive real number.")))
                end
                temp_num = ((pid_data(1,1)*pid_data(3,1)/pid_data(4,1))+pid_data(1,1)*pid_data(3,1))*s^2+((pid_data(1,1)*pid_data(3,1))/(pid_data(2,1)*pid_data(4,1))+pid_data(1,1))*s+pid_data(1,1)/pid_data(2,1)
                temp_den = (pid_data(3,1)/pid_data(4,1))*s^2+s
                simp_mode(%t)
                if pid_data(1,1) == 0 then
                    tf_data(ii,1)=0
                else
                    temp_tf = temp_num/temp_den
                    tf_data(ii,1) = syslin('c',temp_tf)                   //formation of the continuous time transfer function
                end
                simp_mode(%f)
                temp_Kp_data(ii,1) = pid_data(1,1)
                temp_Ki_data(ii,1) = pid_data(2,1)
                temp_Kd_data(ii,1) = pid_data(3,1)
                temp_N_data(ii,1) = pid_data(4,1)
            end//end of continuous domain
//            disp(tf_data)
                if current_length<=2 then
                    output_data = hypermat([current_size],tf_data(:,1))
                elseif current_length > 2 then
                    output_data = hypermat([current_size])
                    output_data(:,:,:,:)=tf_data
                end
                output = output_data
//=================================================================================================================================================================================
        elseif datatype == 1 & inctd == 0 & indtd == 1 then                                                          //discrete time domain 
            printf('\n Discrete Domain \n')
        //now find the type of discrete time integral and derivative formula
                if size(varargin_data(5))~=[1 1] then
                    error(msprintf(gettext("Please enter valid sampling timing. " )))
                end
//**********************************************************************************************************************************************************************
                //now find the type of discrete time integral and derivative formula
                Iformula_index = 0
                Dformula_index = 0
                sysTs = varargin_data(5)
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
                    temp_N = varargin_data(4)
                    temp_Kd = varargin_data(3)
                    size_N = size(temp_N)
                    if temp_Kd ~= 0 & temp_N ~= 0 then
                            if size_N == [1 1] then
                                temp_index = 1
                            else
                                temp_index = total_num_element 
                            end
                            for ii =1:temp_index
                                  if size(temp_N) == [1 1] then
                                    temp_N_data = temp_N
                                else
                                    temp_N_data = temp_N(ii)
                                end
                                if size(temp_Kd)== [1 1] then
                                    temp_Kd_data = temp_Kd
                                else
                                    temp_Kd_data = temp_Kd(ii)
                                end
                                if temp_N_data ~= %inf then
                                    if (temp_Kd_data/temp_N_data) < varargin_data(5)/2 then//varargin_data(5)>=2*temp_Tf(ii) then
                                        error(msprintf(gettext("Sample timing must be greater then two times of filter timing Td/N>Ts/2")))
                                    end
                                end
                            end
                     end
                    Dformula = varargin_data(5)/(z-1)
                end 
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

                if phasemag(pid_data(1,1)) ~= 0 then    //isequal(p,0) == 'T' then          //& pid_data(1,1) < 0 then
                     error(msprintf(gettext("Proportional gain (Kp) must be a finite positive real number.")))
                elseif  pid_data(2,1) == 0 | phasemag(pid_data(2,1)) ~= 0 then
                    error(msprintf(gettext("Integral time (Ti) must be a finite positive real number.")))
                elseif phasemag(pid_data(3,1)) ~= 0 then
                    error(msprintf(gettext("Derivative time (Ti) must be a finite positive real number.")))
                elseif pid_data(4,1) == 0 | phasemag(pid_data(4,1)) ~= 0 then
                    error(msprintf(gettext("Derivative filter divisor(N) must be a finite positive real number.")))
                elseif phasemag(varargin_data(5)) ~= 0 then
                    error(msprintf(gettext("Sampling time (Ts) must be a finite positive real number.")))   
                end
                temp_Ti = (Iformula)*(1/pid_data(2,1))
                temp_Td = (1/((pid_data(3,1)/pid_data(4,1))+Dformula))*(pid_data(3,1))
                temp_pid = (1+temp_Ti+temp_Td)
                temp_pid = pid_data(1,1)*temp_pid
                simp_mode(%t)
                tf_data(ii,1) = simp(syslin('d',temp_pid))                            //formation of the discrete time transfer function 
                simp_mode(%f)
                temp_Kp_data(ii,1) = pid_data(1,1)
                temp_Ki_data(ii,1) = pid_data(2,1)
                temp_Kd_data(ii,1) = pid_data(3,1)
                temp_N_data(ii,1) = pid_data(4,1)
            end//end of continuous domain
                if current_length<=2 then
                    output_data = hypermat([current_size],tf_data(:,1))
                elseif current_length > 2 then
                    output_data = hypermat([current_size])
                    output_data(:,:,:,:)=tf_data
                end
                output_data.dt = varargin_data(5)
                output = output_data
         // end of discrete domain
//=====================================================================================================================================================================================
        elseif datatype == 2 then                                                                           //transfer function domain
            printf('\n rational data \n')
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
                    //sysTf_data = systf_num/systf_den
                    sysTs = varargin_data(1).dt
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
//___________________________________________________________________________________________________________________________________________________________________________
                        if length_num_coeff == 1 & length_den_coeff == 1 then
                            temp_Kp_data(ii,1) = last_num_coeff*num_coeff(1)/last_den_coeff
                            temp_Ki_data(ii,1) = 0
                            temp_Kd_data(ii,1) = 0
                            temp_N_data(ii,1) = 0
                        //Kp*(1+1/(Ti*s))
                        elseif length_num_coeff == 2 & length_den_coeff == 2 & den_coeff(1) == 0 then
                            temp_Kp_data(ii,1) = last_num_coeff/last_den_coeff
                            temp_Ki_data(ii,1) = 1/num_coeff(1)
                            temp_Kd_data(ii,1) = 0
                            temp_N_data(ii,1) = %inf
                        //Kp*(1+1/(Ti*s)+Td*s/((Td/N)*s+1))
                        elseif length_den_coeff == 3 &length_num_coeff == 3 & den_coeff(1) == 0 then
                            A = [1/den_coeff(2) -num_coeff(2);
                                 1 -num_coeff(1)]
                            b = [num_coeff(2)/den_coeff(2) - 1;num_coeff(1)/den_coeff(2)]
                            temp_Ki_n_Td = A\b
                            temp_Ki_data(ii,1) = 1/temp_Ki_n_Td(1)
                            temp_Kd_data(ii,1) = temp_Ki_n_Td(2)
                            temp_N_data(ii,1) = temp_Kd_data(ii,1)*den_coeff(2)
                            temp_Kp_data(ii,1) = (last_num_coeff/last_den_coeff)/(temp_N_data(ii,1)+1)
                        //Kp*(1+Td*s/((Td/N)*s+1))
                        elseif length_num_coeff == 2 & length_den_coeff == 2 then
                            temp_Kd_data(ii,1) = 1/num_coeff(1) - 1/den_coeff(1)
                            temp_N_data(ii,1) = temp_Kd_data(ii,1)*den_coeff(1)
                            temp_Kp_data(ii,1) = last_num_coeff/(last_den_coeff*(1+temp_N_data(ii,1)))
                            temp_Ki_data(ii,1) = %inf
                        //Kp*(1+Td*s)
                        elseif length_num_coeff == 2 & length_den_coeff == 1 then
                            temp_Kd_data(ii,1) = 1/num_coeff(1)
                            temp_Kp_data(ii,1) = last_num_coeff*num_coeff(1)/last_den_coeff
                            temp_Ki_data(ii,1) = %inf
                            temp_N_data(ii,1) = %inf
                       elseif length_num_coeff == 3 & length_den_coeff == 2 & den_coeff(1) == 0 then
                            temp_Kd_data(ii,1) = 1/num_coeff(2)
                            temp_Ki_data(ii,1) = num_coeff(2)/num_coeff(1)
                            temp_Kp_data(ii,1) = num_coeff(2)*last_num_coeff/last_den_coeff
                            temp_N_data(ii,1) = %inf    
                       elseif (length_num_coeff == 1 & length_den_coeff == 2 & den_coeff(1) == 0)|(length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1) == 0)| (length_num_coeff == 2 & length_den_coeff == 1 & num_coeff(1) == 0) then
                            error(msprintf(gettext("The given control tranfer function is in Parallel PID form. Use pid function.")))
                        else
                            disp(systf_num(ii)/systf_den(ii))
                            error(msprintf(gettext("Above equation is not in parallel PID form")))
                        end
                        output = varargin_data(1)
//____________________________________________________________________________________________________________________________________________________________________________
                        end
                    elseif varargin_data(1).dt == 'd' | typeof(varargin_data(1).dt)=='constant' then                //discrete time domain
                        if varargin_data(1).dt == 'd' then
                            varargin_data(1).dt = 1
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
                    //end of type of Iformula and Dformula
                    Iformula_method = Iformula
                    Dformula_method = Dformula
//extract the element from the given trnsfer function
//____________________________________________________________________________________________________________________________________________________________________________
                    for ii = 1:total_num_element
//------------------------------------------------------------------------------
// finding out the transfer function 
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

//***************************************************************************************************************************************************************
                            if Iformula =='F' & Dformula == 'F' then
                            // when Ifomula = Forward Euler and Dformula = Forward Euler
                                 if length_num_coeff == 1 & length_den_coeff == 1 then
                                //For P type controller
                                    temp_Kp_data(ii,1) = last_num_coeff/last_den_coeff
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                                 elseif length_num_coeff == 2 & length_den_coeff == 2 & den_coeff(1) == -1 then
                                // For PI type controller
                                     temp_Kp_data(ii,1) = last_num_coeff/last_den_coeff
                                     temp_Ki_data(ii,1) = sysTs/(1+num_coeff(1))
                                     temp_Kd_data(ii,1) = 0
                                     temp_N_data(ii,1) = 0
                                 elseif length_num_coeff == 3 & length_den_coeff == 3 & round(den_coeff(2)+den_coeff(1)) == -1 then
                                 // PID N ~= %inf
                                     A = [(2+num_coeff(2)) -sysTs;
                                          (1-num_coeff(1)) sysTs*(1+den_coeff(2))]
                                     b = [den_coeff(2)-num_coeff(2);1+den_coeff(2)+num_coeff(1)]
                                     temp_pid = A\b
                                     temp_N_data(ii,1) = temp_pid(1)
                                     temp_Ki_data(ii,1) = 1/temp_pid(2)
                                     temp_Kp_data(ii,1) = last_num_coeff/(last_den_coeff*(1+temp_N_data(ii,1)))
                                     temp_Kd_data(ii,1) = temp_N_data(ii,1)*sysTs/(2+den_coeff(2))
                                 elseif length_num_coeff == 3 & length_den_coeff == 2 & den_coeff(1) == -1 then
                                  // PID N == %inf
                                     temp_Kd_data(ii,1) = sysTs/(2+num_coeff(2))
                                     temp_Ki_data(ii,1) = sysTs^2/(temp_Kd_data(ii,1)*(num_coeff(1)-1)+sysTs)
                                     temp_Kp_data(ii,1) = sysTs*last_num_coeff/(last_den_coeff*temp_Kd_data(ii,1))
                                     temp_N_data(ii,1) = %inf
                                 elseif length_num_coeff == 2 & length_den_coeff == 2 then
                                  //PD N~= %inf
                                     temp_N_data(ii,1) = (den_coeff(1)-num_coeff(1))/(1+num_coeff(1))
                                     temp_Kd_data(ii,1) = temp_N_data(ii,1)*sysTs/(1+den_coeff(1))
                                     temp_Kp_data(ii,1) = last_num_coeff/((1+temp_N_data(ii,1))*last_den_coeff)
                                     temp_Ki_data(ii,1) = %inf
                                 elseif length_num_coeff == 2 & length_den_coeff == 1 then
                                  // PD N == %inf
                                     temp_Kd_data(ii,1) = sysTs/(1+num_coeff(1))
                                     temp_Kp_data(ii,1) = sysTs*last_num_coeff/(temp_Kd_data(ii,1)*last_den_coeff)
                                     temp_N_data(ii,1) = %inf
                                     temp_Ki_data(ii,1) = %inf
                                 elseif (length_num_coeff == 1 & length_den_coeff == 2) | (length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1) == -1) | (length_num_coeff == 2 & length_den_coeff == 1 & num_coeff(1) == -1) then
                                  // Parallel form of controller I type controller, D type of controller
                                     error(msprintf(gettext("The given control tranfer function is in Parallel PID form. Use pid function for validation.")))
                                 else
                                     disp(systf_num(ii)/systf_den(ii))
                                     error(msprintf(gettext("Above equation is not a valid standard PID controller")))
                                 end // end of when Ifomula = Forward Euler and Dformula = Forward Euler
//***********************************************************************************************************************************************************************
                            elseif Iformula =='F' & Dformula == 'B' then
                            // when Ifomula = Forward Euler and Dformula = Backward Euler
                                 if length_num_coeff == 1 & length_den_coeff == 1 then
                                //For P type controller
                                    temp_Kp_data(ii,1) = last_num_coeff/last_den_coeff
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                                 elseif length_num_coeff == 2 & length_den_coeff == 2 & den_coeff(1) == -1 then
                                // For PI type controller
                                     temp_Kp_data(ii,1) = last_num_coeff/last_den_coeff
                                     temp_Ki_data(ii,1) = sysTs/(1+num_coeff(1))
                                     temp_Kd_data(ii,1) = 0
                                     temp_N_data(ii,1) = 0
                                 elseif length_num_coeff == 3 & length_den_coeff == 3 & den_coeff(2) == -1 & den_coeff(1) == 0 then
                                  // PID N == %inf
                                      temp_Kd_data(ii,1) = sysTs*num_coeff(1)/(1-num_coeff(1))
                                      temp_Ki_data(ii,1) = sysTs^2/(temp_Kd_data(ii,1)*(2+num_coeff(2))+sysTs*(1+num_coeff(2)))
                                      temp_Kp_data(ii,1) = sysTs*last_num_coeff/((temp_Kd_data(ii,1)+sysTs)*last_den_coeff)
                                      temp_N_data(ii,1) = %inf
                                 elseif length_num_coeff == 3 & length_den_coeff == 3 & round((2+den_coeff(2))*den_coeff(1)) == round((1-den_coeff(1))*(1+den_coeff(2))) then
                                 // PID N ~= %inf
                                     printf('\n N is not inf')
                                     temp_data = (1-den_coeff(1))/den_coeff(1)
                                     A = [sysTs*(1+temp_data) -(2+num_coeff(2));
                                         -sysTs 1-num_coeff(1)]
                                     b = [temp_data+(1+temp_data)*num_coeff(2)+2;(1+temp_data)*num_coeff(1)-1]
                                     temp_pid = A\b
                                     temp_Ki_data(ii,1) = 1/temp_pid(1)
                                     temp_N_data(ii,1) = temp_pid(2)
                                     temp_Kd_data(ii,1) = temp_N_data(ii,1)*sysTs/temp_data
                                     temp_Kp_data(ii,1) = (1+temp_data)*last_num_coeff/((1+temp_N_data(ii,1)+temp_data)*last_den_coeff)
                                     
                                 
                                 elseif (length_num_coeff == 1 & length_den_coeff == 2) | (length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1) == -1 & den_coeff(1) < 0 ) | (length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1) == -1 & den_coeff(1) == 0) then
                                     // Parallel form of controller I type controller, D type of controller
                                     error(msprintf(gettext("The given control tranfer function is in Parallel PID form. Use pid function for validation.")))
                                     
                                     
                                 elseif length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1) < 0 & den_coeff(1) < 0 then
                                 // PD controller N~=%inf
                                     temp_data = (1+den_coeff(1))/den_coeff(1)
                                     temp_N_data(ii,1) = -(1+(1-temp_data)*num_coeff(1))/(1+num_coeff(1))
                                     temp_Kd_data(ii,1) = -temp_N_data(ii,1)*sysTs/temp_data
                                     temp_Kp_data(ii,1) = (temp_Kd_data(ii,1)+temp_N_data(ii,1)*sysTs)*last_num_coeff/(last_den_coeff*(temp_Kd_data(ii,1)+temp_Kd_data(ii,1)*temp_N_data(ii,1)+temp_N_data(ii,1)*sysTs))
                                     temp_Ki_data(ii,1) = %inf
                                     
                                 elseif length_num_coeff == 2 & length_den_coeff == 2 & den_coeff(1) == 0 & num_coeff(1) < 0 then
                                 // PD controller N ~= %inf 
                                     temp_Kd_data(ii,1) = -sysTs*num_coeff(1)/(1+num_coeff(1))
                                     temp_Kp_data(ii,1) = sysTs*last_num_coeff/((temp_Kd_data(ii,1)+sysTs)*last_den_coeff)
                                     temp_N_data(ii,1) = %inf
                                     temp_Ki_data(ii,1) = %inf
                                     
                                 else
                                     disp(systf_num(ii)/systf_den(ii))
                                     error(msprintf(gettext("Above equation is not a valid standard PID controller")))
                                 end // end of when Ifomula = Forward Euler and Dformula = Backward Euler
//***********************************************************************************************************************************************************************                 
                            elseif Iformula =='F' & Dformula == 'T' then
                            // when Ifomula = Forward Euler and Dformula = Trapezoidal
                                 if length_num_coeff == 1 & length_den_coeff == 1 then
                                //For P type controller
                                     temp_Kp_data(ii,1) = last_num_coeff/last_den_coeff
                                     temp_Ki_data(ii,1) = 0
                                     temp_Kd_data(ii,1) = 0
                                     temp_Tf_data(ii,1) = 0
                                 elseif length_num_coeff == 2 & length_den_coeff == 2 & den_coeff(1) == -1 then
                                // For PI type controller
                                     temp_Kp_data(ii,1) = last_num_coeff/last_den_coeff
                                     temp_Ki_data(ii,1) = sysTs/(1+num_coeff(1))
                                     temp_Kd_data(ii,1) = 0
                                     temp_N_data(ii,1) = 0
                                 elseif length_num_coeff == 3 & length_den_coeff == 3 & round((2+den_coeff(2))*(den_coeff(1)+1)) == round((den_coeff(1)-1)*den_coeff(2)) then
                                 // For PID type controller N ~= %inf    
                                     temp_data = 2*(1-den_coeff(1))/(1+den_coeff(1))
                                     A = [(sysTs*temp_data+2*sysTs) -(4+2*num_coeff(2));
                                          (sysTs*temp_data-2*sysTs) 2*(1-num_coeff(1))]
                                     b = [4+(2+temp_data)*num_coeff(2);(2+temp_data)*num_coeff(1)+temp_data-2]
                                     temp_pid = A\b
                                     temp_Ki_data(ii,1) = 1/temp_pid(1)
                                     temp_N_data(ii,1) = temp_pid(2)
                                     temp_Kd_data(ii,1) = (sysTs*temp_N_data(ii,1))/temp_data
                                     temp_Kp_data(ii,1) = ((2+temp_data)*last_num_coeff)/((2+2*temp_N_data(ii,1)+temp_data)*last_den_coeff)
                                     
                                 elseif length_den_coeff == 3 & length_num_coeff == 3 & den_coeff(2) == 0 & den_coeff(1) == -1 then
                                 // For PID type controller N == %inf
                                     A = [sysTs^2 -(4+2*num_coeff(2));
                                          sysTs^2 2*(1-num_coeff(1))]
                                     b = [sysTs*num_coeff(2);sysTs*(1+num_coeff(1))]
                                     temp_pid = A\b
                                     temp_Ki_data(ii,1) = 1/temp_pid(1)
                                     temp_Kd_data(ii,1) = temp_pid(2)
                                     temp_Kp_data(ii,1) = sysTs*last_num_coeff/((2*temp_Kd_data(ii,1)+sysTs)*last_den_coeff)
                                     temp_N_data(ii,1) = %inf
                                     
                                 elseif (length_num_coeff == 1 & length_den_coeff == 2) | (length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1) == -1) | (length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1) == -1 & den_coeff(1) == 1) then
                                     // Parallel form of controller I type controller, D type of controller
                                     error(msprintf(gettext("The given control tranfer function is in Parallel PID form. Use pid function for validation.")))
                                        
                                 elseif length_num_coeff == 2 & length_den_coeff == 2 & den_coeff(1) < 0 & num_coeff(1) < 0 then
                                  // PD type controller N ~= %inf    
                                      temp_data = 2*(1+den_coeff(1))/(1-den_coeff(1))
                                      temp_N_data(ii,1) = -((2-temp_data)+(2+temp_data)*num_coeff(1))/(2*(1+num_coeff(1)))
                                      temp_Kd_data(ii,1) = temp_N_data(ii,1)*sysTs/temp_data
                                      temp_Kp_data(ii,1) = ((2+temp_data)*last_num_coeff)/((2+2*temp_N_data(ii,1)+temp_data)*last_den_coeff)
                                      temp_Ki_data(ii,1) = %inf
                                  
                                 elseif length_num_coeff == 2 & length_den_coeff == 2 then // & den_coeff(1) == 1      
                                  // PD type controller N == %inf    
                                      temp_Kd_data(ii,1) = sysTs*(1-num_coeff(1))/(2*(1+num_coeff(1)))
                                      temp_Kp_data(ii,1) = sysTs*last_num_coeff/((2*temp_Kd_data(ii,1)+sysTs)*last_den_coeff)
                                      temp_Ki_data(ii,1) = %inf
                                      temp_N_data(ii,1) = %inf
                                 else
                                     disp(systf_num(ii)/systf_den(ii))
                                     error(msprintf(gettext("Above equation is not a valid standard PID controller")))
                                 end // end of when Iformula = Forward Euler and Dformula = Trapezoidal
                                                             
//************************************************************************************************************************************************************************
                            elseif Iformula =='B' & Dformula == 'F' then
                            // When Iformula = Backward Euler and Dformula = Forward Euler    
                                 if length_num_coeff == 1 & length_den_coeff == 1 then
                                 //For P type controller
                                     temp_Kp_data(ii,1) = last_num_coeff/last_den_coeff
                                     temp_Ki_data(ii,1) = 0
                                     temp_Kd_data(ii,1) = 0
                                     temp_Tf_data(ii,1) = 0
                                     
                                 elseif length_num_coeff == 2 & length_den_coeff == 2 & den_coeff(1) == -1 & num_coeff(1) < 0 then
                                 // For PI type controller     
                                     temp_Ki_data(ii,1) = -sysTs/(1+1/num_coeff(1))
                                     temp_Kp_data(ii,1) = -num_coeff(1)*last_num_coeff/last_den_coeff
                                     temp_Kd_data(ii,1) = 0
                                     temp_N_data(ii,1) = %inf
                                     
                                 elseif length_num_coeff == 3 & length_den_coeff == 3 & round(den_coeff(1)+den_coeff(2)) == -1 then    
                                 // PID type controller N~= %inf
                                     temp_data = 2+den_coeff(2)
                                     A = [sysTs*(temp_data-1-num_coeff(2)) -(2+num_coeff(2));
                                          -sysTs*num_coeff(1) 1-num_coeff(1)]
                                     b = [num_coeff(2)+2-temp_data; temp_data-1+num_coeff(1)]
                                     temp_pid = A\b
                                     temp_Ki_data(ii,1) = 1/temp_pid(1)
                                     temp_N_data(ii,1) = temp_pid(2)
                                     temp_Kd_data(ii,1) = temp_N_data(ii,1)*sysTs/temp_data
                                     temp_Kp_data(ii,1) = last_num_coeff/((1+temp_N_data(ii,1)+sysTs/temp_Ki_data(ii,1))*last_den_coeff)
                                     
                                 elseif length_num_coeff == 3 & length_den_coeff == 2 & den_coeff(1) == -1 then
                                 // PID type controller N == %inf    
                                     temp_Kd_data(ii,1) = sysTs/(1-num_coeff(1))
                                     temp_Ki_data(ii,1) = sysTs^2/(temp_Kd_data(ii,1)*(2+num_coeff(2))-sysTs)
                                     temp_Kp_data(ii,1) = sysTs*last_num_coeff/(temp_Kd_data(ii,1)*last_den_coeff)
                                     temp_N_data(ii,1) = %inf     
                                     
                                 elseif (length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1) == 0 & den_coeff(1) == -1) | (length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1) == -1) | (length_num_coeff == 2 & length_den_coeff == 1 & num_coeff(1) == -1) then
                                     // Parallel form of controller I type controller, D type of controller
                                     error(msprintf(gettext("The given control tranfer function is in Parallel PID form. Use pid function for validation.")))
                                          
                                 elseif length_num_coeff == 2 & length_den_coeff == 2 then
                                 // PD type controller N ~= %inf
                                     temp_data = 1+den_coeff(1)
                                     temp_N_data(ii,1) = (temp_data-1-num_coeff(1))/(1+num_coeff(1))
                                     temp_Kd_data(ii,1) = temp_N_data(ii,1)*sysTs/temp_data
                                     temp_Kp_data(ii,1) = last_num_coeff/((1+temp_N_data(ii,1))*last_den_coeff)
                                     temp_Ki_data(ii,1) = %inf
                                     
                                 elseif length_num_coeff == 2 & length_den_coeff == 1 then
                                 // PD type controller N == %inf
                                     temp_Kd_data(ii,1) = sysTs/(1+num_coeff(1))
                                     temp_Kp_data(ii,1) = (1+num_coeff(1))*last_num_coeff/last_den_coeff
                                     temp_Ki_data(ii,1) = %inf
                                     temp_N_data(ii,1) = %inf
                                     
                                 else
                                     disp(systf_num(ii)/systf_den(ii))
                                     error(msprintf(gettext("Above equation is not a valid standard PID controller")))
                                 end // end of when Iformula = Backward Euler and Dformula = Forward Euler    
//************************************************************************************************************************************************************************             
                            elseif Iformula =='B' & Dformula == 'B' then
                            //when Iformula = Backward Euler and Dformula = Backward Euler
                                if length_num_coeff == 1 & length_den_coeff == 1 then
                                //For P type controller
                                     temp_Kp_data(ii,1) = last_num_coeff/last_den_coeff
                                     temp_Ki_data(ii,1) = 0
                                     temp_Kd_data(ii,1) = 0
                                     temp_Tf_data(ii,1) = 0
                                     
                                 elseif length_num_coeff == 2 & length_den_coeff == 2 & den_coeff(1) == -1 & num_coeff(1) < 0 then
                                 // For PI type controller     
                                     temp_Ki_data(ii,1) = -sysTs/(1+1/num_coeff(1))
                                     temp_Kp_data(ii,1) = -num_coeff(1)*last_num_coeff/last_den_coeff
                                     temp_Kd_data(ii,1) = 0
                                     temp_N_data(ii,1) = %inf
                                     
                                     
                                 elseif length_num_coeff == 3 & length_den_coeff == 3 & den_coeff(1) == 0 & den_coeff(2) == -1 then 
                                 // PID type controller N== %inf
                                     printf('\n hi ')
                                     A = [-(sysTs^2*num_coeff(2)) -(2+num_coeff(2));
                                          -sysTs^2*num_coeff(1) (1-num_coeff(1))]
                                     b = [sysTs+sysTs*num_coeff(2);sysTs*num_coeff(1)]
                                     temp_pid = A\b
                                     temp_Ki_data(ii,1) = 1/temp_pid(1)
                                     temp_Kd_data(ii,1) = temp_pid(2)
                                     temp_Kp_data(ii,1) = sysTs*last_num_coeff/((temp_pid(1)*sysTs^2+sysTs+temp_Kd_data(ii,1))*last_den_coeff)
                                     temp_N_data(ii,1) = %inf
                                 
                                 elseif length_num_coeff == 3 & length_den_coeff == 3 & round((2+den_coeff(2))*den_coeff(1)) == round((den_coeff(1)-1)*(den_coeff(2)+1)) then
                                 // PID type controller N~= %inf
                                     temp_data = (1-den_coeff(1))/den_coeff(1)
                                     A = [-sysTs*(1+(1+temp_data)*num_coeff(2)) -(2+num_coeff(2));
                                          -sysTs*(1+temp_data)*num_coeff(1) 1-num_coeff(1)]
                                     b = [(1+temp_data)*num_coeff(2)+(2+temp_data);(1+temp_data)*num_coeff(1)-1]            
                                     temp_pid = A\b
                                     temp_Ki_data(ii,1) = 1/temp_pid(1)
                                     temp_N_data(ii,1) = temp_pid(2)
                                     temp_Kd_data(ii,1) = temp_N_data(ii,1)*sysTs/temp_data
                                     temp_Kp_data(ii,1) = (1+temp_data)*last_num_coeff/((1+temp_N_data(ii,1)+temp_data+temp_pid(1)*sysTs*(1+temp_data))*last_den_coeff)
                                          
                                 elseif (length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1) == 0 & den_coeff(1) == -1) | (length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1) == -1 & den_coeff(1) < 0 ) | (length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1) == -1 & den_coeff(1) == 0) then
                                 // Parallel form of controller I type controller, D type of controller
                                     error(msprintf(gettext("The given control tranfer function is in Parallel PID form. Use pid function for validation.")))    
                                     
                                 elseif length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1) < 0 & den_coeff(1) < 0 then
                                 // PD controller N~=%inf
                                     temp_data = (1+den_coeff(1))/den_coeff(1)
                                     temp_N_data(ii,1) = -(1+(1-temp_data)*num_coeff(1))/(1+num_coeff(1))
                                     temp_Kd_data(ii,1) = -temp_N_data(ii,1)*sysTs/temp_data
                                     temp_Kp_data(ii,1) = (temp_Kd_data(ii,1)+temp_N_data(ii,1)*sysTs)*last_num_coeff/(last_den_coeff*(temp_Kd_data(ii,1)+temp_Kd_data(ii,1)*temp_N_data(ii,1)+temp_N_data(ii,1)*sysTs))
                                     temp_Ki_data(ii,1) = %inf
                                     
                                 elseif length_num_coeff == 2 & length_den_coeff == 2 & den_coeff(1) == 0 & num_coeff(1) < 0 then
                                 // PD controller N ~= %inf 
                                     temp_Kd_data(ii,1) = -sysTs*num_coeff(1)/(1+num_coeff(1))
                                     temp_Kp_data(ii,1) = sysTs*last_num_coeff/((temp_Kd_data(ii,1)+sysTs)*last_den_coeff)
                                     temp_N_data(ii,1) = %inf
                                     temp_Ki_data(ii,1) = %inf
                                     
                                 else
                                     disp(systf_num(ii)/systf_den(ii))
                                     error(msprintf(gettext("Above equation is not a valid standard PID controller")))         
                                 end // end of when Iformula = Backward Euler and Dformula = Backward Euler
                            
//************************************************************************************************************************************************************************             
                            elseif Iformula =='B' & Dformula == 'T' then
                            //when Iformula = Backward Euler and Dformula = Trapezoidal
                                if length_num_coeff == 1 & length_den_coeff == 1 then
                                //For P type controller
                                     temp_Kp_data(ii,1) = last_num_coeff/last_den_coeff
                                     temp_Ki_data(ii,1) = 0
                                     temp_Kd_data(ii,1) = 0
                                     temp_Tf_data(ii,1) = 0
                                     
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & den_coeff(1) == -1 & num_coeff(1) < 0 then
                                 // For PI type controller     
                                     temp_Ki_data(ii,1) = -sysTs/(1+1/num_coeff(1))
                                     temp_Kp_data(ii,1) = -num_coeff(1)*last_num_coeff/last_den_coeff
                                     temp_Kd_data(ii,1) = 0
                                     temp_N_data(ii,1) = %inf
                                     
                                elseif length_num_coeff == 3 & length_den_coeff == 3 & den_coeff(2) == 0 & den_coeff(1) == -1 then 
                                // For PID type controller N == %inf
                                     A = [sysTs^2*(1-num_coeff(2)) -(4+2*num_coeff(2));
                                          -sysTs^2*num_coeff(1) 2*(1-num_coeff(1))]
                                     b = [sysTs*num_coeff(2);sysTs*(1+num_coeff(1))]
                                     temp_pid = A\b
                                     temp_Ki_data(ii,1) = 1/temp_pid(1)
                                     temp_Kd_data(ii,1) = temp_pid(2)
                                     temp_Kp_data(ii,1) = sysTs*last_num_coeff/((temp_pid(1)*sysTs^2+sysTs+2*temp_Kd_data(ii,1))*last_den_coeff)
                                     temp_N_data(ii,1) = %inf
                                          
                                elseif length_num_coeff == 3 & length_den_coeff == 3 & round((2+den_coeff(2))*(1+den_coeff(1)))==round((den_coeff(1)-1)*den_coeff(2)) then
                                // For PID type controller N ~= %inf
                                     temp_data = 2*(1-den_coeff(1))/(1+den_coeff(1))
                                     A = [sysTs*((temp_data-2)-(2+temp_data)*num_coeff(2)) -(4+2*num_coeff(2));
                                          -sysTs*(2+temp_data)*num_coeff(1) 2*(1-num_coeff(1))]
                                     b = [4+(2+temp_data)*num_coeff(2);(2+temp_data)*num_coeff(1)+temp_data-2]
                                     temp_pid = A\b
                                     temp_Ki_data(ii,1) = 1/temp_pid(1)
                                     temp_N_data(ii,1) = temp_pid(2)
                                     temp_Kd_data(ii,1) = temp_N_data(ii,1)*sysTs/temp_data
                                     temp_Kp_data(ii,1) = (2+temp_data)*last_num_coeff/((2+temp_data+2*temp_N_data(ii,1)+temp_pid(1)*(2+temp_data)*sysTs)*last_den_coeff)
                                     
                                elseif (length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1) == 0 & den_coeff(1) == -1) | (length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1) == -1) | (length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1) == -1 & den_coeff(1) == 1) then
                                // Parallel form of controller I type controller, D type of controller
                                     error(msprintf(gettext("The given control tranfer function is in Parallel PID form. Use pid function for validation.")))
                                        
                                 elseif length_num_coeff == 2 & length_den_coeff == 2 & den_coeff(1) < 0 & num_coeff(1) < 0 then
                                  // PD type controller N ~= %inf    
                                      temp_data = 2*(1+den_coeff(1))/(1-den_coeff(1))
                                      temp_N_data(ii,1) = -((2-temp_data)+(2+temp_data)*num_coeff(1))/(2*(1+num_coeff(1)))
                                      temp_Kd_data(ii,1) = temp_N_data(ii,1)*sysTs/temp_data
                                      temp_Kp_data(ii,1) = ((2+temp_data)*last_num_coeff)/((2+2*temp_N_data(ii,1)+temp_data)*last_den_coeff)
                                      temp_Ki_data(ii,1) = %inf
                                  
                                 elseif length_num_coeff == 2 & length_den_coeff == 2 then // & den_coeff(1) == 1      
                                  // PD type controller N == %inf    
                                      temp_Kd_data(ii,1) = sysTs*(1-num_coeff(1))/(2*(1+num_coeff(1)))
                                      temp_Kp_data(ii,1) = sysTs*last_num_coeff/((2*temp_Kd_data(ii,1)+sysTs)*last_den_coeff)
                                      temp_Ki_data(ii,1) = %inf
                                      temp_N_data(ii,1) = %inf
                                 else
                                     disp(systf_num(ii)/systf_den(ii))
                                     error(msprintf(gettext("Above equation is not a valid standard PID controller")))          
                                end // end of when Iformula = Backward Euler and Dformula = Trapezoidal
//************************************************************************************************************************************************************************
                            elseif Iformula =='T' & Dformula == 'F' then
                            // when Iformula = Trapezoidal and Dformula = Forward Euler
                                if length_num_coeff == 1 & length_den_coeff == 1 then
                                //For P type controller
                                    temp_Kp_data(ii,1) = last_num_coeff/last_den_coeff
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                                     
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & den_coeff(1) == -1 then
                                // For PI type controller
                                    temp_Ki_data(ii,1) = sysTs*(1-num_coeff(1))/(2*(1+num_coeff(1)))
                                    temp_Kp_data(ii,1) = 2*last_num_coeff/((2+ sysTs/temp_Ki_data(ii,1))*last_den_coeff)
                                    temp_Kd_data(ii,1) = 0
                                    temp_N_data(ii,1) = %inf
                                                                             
                                elseif length_num_coeff == 3 & length_den_coeff == 2 & den_coeff(1) == -1 then
                                // For PID type controller  N==%inf
                                    A = [sysTs^2 -(4+2*num_coeff(2));
                                         sysTs^2 2*(1-num_coeff(1))]
                                    b = [-2*sysTs;2*sysTs]
                                    temp_pid = A\b
                                    temp_Ki_data(ii,1) = 1/temp_pid(1)
                                    temp_Kd_data(ii,1) = temp_pid(2)
                                    temp_Kp_data(ii,1) = sysTs*last_num_coeff/(temp_Kd_data(ii,1)*last_den_coeff)
                                    temp_N_data(ii,1) = %inf
                                         
                                elseif length_num_coeff == 3 & length_den_coeff == 3 & round(den_coeff(1)+den_coeff(2)) == -1 then
                                // For PID type controller N ~=%inf
                                     temp_data = 1-den_coeff(1)
                                     A = [sysTs*(temp_data-num_coeff(2)) -2*(2+num_coeff(2));
                                          sysTs*(temp_data-1-num_coeff(1)) 2*(1-num_coeff(1))]
                                     b = [(4-2*temp_data)+2*num_coeff(2); 2*(temp_data-1)+2*num_coeff(1)]
                                     temp_pid = A\b
                                     temp_Ki_data(ii,1) = 1/temp_pid(1)
                                     temp_N_data(ii,1) = temp_pid(2)
                                     temp_Kd_data(ii,1) = (temp_N_data(ii,1)*sysTs)/temp_data
                                     temp_Kp_data(ii,1) = 2*last_num_coeff/((2+2*temp_N_data(ii,1)+temp_pid(1)*sysTs)*last_den_coeff)
                                elseif (length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1) == 1 & den_coeff(1) == -1) | (length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1) == -1) | (length_num_coeff == 2 & length_den_coeff == 1 & num_coeff(1) == -1) then
                                     // Parallel form of controller I type controller, D type of controller
                                     error(msprintf(gettext("The given control tranfer function is in Parallel PID form. Use pid function for validation.")))
                                          
                                 elseif length_num_coeff == 2 & length_den_coeff == 2 then
                                 // PD type controller N ~= %inf
                                     temp_data = 1+den_coeff(1)
                                     temp_N_data(ii,1) = (temp_data-1-num_coeff(1))/(1+num_coeff(1))
                                     temp_Kd_data(ii,1) = temp_N_data(ii,1)*sysTs/temp_data
                                     temp_Kp_data(ii,1) = last_num_coeff/((1+temp_N_data(ii,1))*last_den_coeff)
                                     temp_Ki_data(ii,1) = %inf
                                     
                                 elseif length_num_coeff == 2 & length_den_coeff == 1 then
                                 // PD type controller N == %inf
                                     temp_Kd_data(ii,1) = sysTs/(1+num_coeff(1))
                                     temp_Kp_data(ii,1) = (1+num_coeff(1))*last_num_coeff/last_den_coeff
                                     temp_Ki_data(ii,1) = %inf
                                     temp_N_data(ii,1) = %inf
                                     
                                 else
                                     disp(systf_num(ii)/systf_den(ii))
                                     error(msprintf(gettext("Above equation is not a valid standard PID controller")))          
                                end// end of when Iformula = Trapezoidal and Dformula = Forward Euler
                                
//************************************************************************************************************************************************************************
                            elseif Iformula =='T' & Dformula == 'B' then
                            // when Iformula = Trapezoidal and Dformula = Backward Euler
                                if length_num_coeff == 1 & length_den_coeff == 1 then
                                //For P type controller
                                    temp_Kp_data(ii,1) = last_num_coeff/last_den_coeff
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                                     
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & den_coeff(1) == -1 then
                                // For PI type controller
                                    temp_Ki_data(ii,1) = sysTs*(1-num_coeff(1))/(2*(1+num_coeff(1)))
                                    temp_Kp_data(ii,1) = 2*last_num_coeff/((2+ sysTs/temp_Ki_data(ii,1))*last_den_coeff)
                                    temp_Kd_data(ii,1) = 0
                                    temp_N_data(ii,1) = %inf
                                    
                                elseif length_num_coeff == 3 & length_den_coeff == 3 & den_coeff(1) == 0 & den_coeff(2) == -1 then
                                // For PID type controller N~=%inf
                                    A = [sysTs^2*(1-num_coeff(2)) -(4+2*num_coeff(2));
                                         sysTs^2*num_coeff(1) 2*(num_coeff(1)-1)]
                                    b = [2*sysTs*(1+num_coeff(2));2*(-sysTs*num_coeff(1))]
                                    temp_pid = A\b
                                    temp_Ki_data(ii,1) = 1/temp_pid(1)
                                    temp_Kd_data(ii,1) = temp_pid(2)
                                    temp_Kp_data(ii,1) = 2*sysTs*last_num_coeff/((temp_pid(1)*sysTs^2+2*sysTs+2*temp_Kd_data(ii,1))*last_den_coeff)
                                    temp_N_data(ii,1) = %inf
                                    
                                elseif length_num_coeff == 3 & length_den_coeff == 3 & round((2+den_coeff(2))*(den_coeff(1))) == round((den_coeff(1)-1)*(den_coeff(2)+1)) then
                                // For PID type controller N ~= %inf
                                    temp_data = (1-den_coeff(1))/den_coeff(1)
                                    A = [sysTs*(temp_data-(1+temp_data)*num_coeff(2)) -(4+2*num_coeff(2));
                                         -sysTs*(1+(1+temp_data)*num_coeff(1)) 2*(1-num_coeff(1))]
                                    b = [2*(temp_data+2)+2*(1+temp_data)*num_coeff(2);2*(1+temp_data)*num_coeff(1)-2]
                                    temp_pid = A\b
                                    temp_Ki_data(ii,1) = 1/temp_pid(1)
                                    temp_N_data(ii,1) = temp_pid(2)
                                    temp_Kd_data(ii,1) = (temp_N_data(ii,1)*sysTs)/temp_data
                                    temp_Kp_data(ii,1) = 2*(1+temp_data)*last_num_coeff/((temp_pid(1)*sysTs*(1+temp_data)+2*temp_N_data(ii,1)+2*(1+temp_data))*last_den_coeff)
                                    
                                elseif (length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1) == 1 & den_coeff(1) == -1) | (length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1) == -1 & den_coeff(1) < 0 ) | (length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1) == -1 & den_coeff(1) == 0) then
                                 // Parallel form of controller I type controller, D type of controller
                                     error(msprintf(gettext("The given control tranfer function is in Parallel PID form. Use pid function for validation.")))    
                                     
                                 elseif length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1) < 0 & den_coeff(1) < 0 then
                                 // PD controller N~=%inf
                                     temp_data = (1+den_coeff(1))/den_coeff(1)
                                     temp_N_data(ii,1) = -(1+(1-temp_data)*num_coeff(1))/(1+num_coeff(1))
                                     temp_Kd_data(ii,1) = -temp_N_data(ii,1)*sysTs/temp_data
                                     temp_Kp_data(ii,1) = (temp_Kd_data(ii,1)+temp_N_data(ii,1)*sysTs)*last_num_coeff/(last_den_coeff*(temp_Kd_data(ii,1)+temp_Kd_data(ii,1)*temp_N_data(ii,1)+temp_N_data(ii,1)*sysTs))
                                     temp_Ki_data(ii,1) = %inf
                                     
                                 elseif length_num_coeff == 2 & length_den_coeff == 2 & den_coeff(1) == 0 & num_coeff(1) < 0 then
                                 // PD controller N ~= %inf 
                                     temp_Kd_data(ii,1) = -sysTs*num_coeff(1)/(1+num_coeff(1))
                                     temp_Kp_data(ii,1) = sysTs*last_num_coeff/((temp_Kd_data(ii,1)+sysTs)*last_den_coeff)
                                     temp_N_data(ii,1) = %inf
                                     temp_Ki_data(ii,1) = %inf
                                     
                                 else
                                     disp(systf_num(ii)/systf_den(ii))
                                     error(msprintf(gettext("Above equation is not a valid standard PID controller")))             
                                end// end of when Iformula = Trapezoidal and Dformula = Backward Euler
                                
//************************************************************************************************************************************************************************

                            elseif Iformula =='T' & Dformula == 'T' then
                            // when Iformula = trapezoidal and Dformula = trapezoidal
                                if length_num_coeff == 1 & length_den_coeff == 1 then
                                //For P type controller
                                    temp_Kp_data(ii,1) = last_num_coeff/last_den_coeff
                                    temp_Ki_data(ii,1) = 0
                                    temp_Kd_data(ii,1) = 0
                                    temp_Tf_data(ii,1) = 0
                                     
                                elseif length_num_coeff == 2 & length_den_coeff == 2 & den_coeff(1) == -1 then
                                // For PI type controller
                                    temp_Ki_data(ii,1) = sysTs*(1-num_coeff(1))/(2*(1+num_coeff(1)))
                                    temp_Kp_data(ii,1) = 2*last_num_coeff/((2+ sysTs/temp_Ki_data(ii,1))*last_den_coeff)
                                    temp_Kd_data(ii,1) = 0
                                    temp_N_data(ii,1) = %inf
                                
                                elseif length_num_coeff == 3 & length_den_coeff == 3 & den_coeff(1) == -1 & den_coeff(2) == 0 then
                                // For PID type controller N == %inf
                                    A =[sysTs^2*(2-num_coeff(2)) -4*(2+num_coeff(2));
                                        sysTs^2*(1-num_coeff(1)) 4*(1-num_coeff(1))]
                                    b = [2*sysTs*num_coeff(2);2*sysTs*(1+num_coeff(1))]
                                    temp_pid = A\b
                                    temp_Ki_data(ii,1) = 1/temp_pid(1)
                                    temp_Kd_data(ii,1) = temp_pid(2)
                                    temp_Kp_data(ii,1) = 2*sysTs*last_num_coeff/((temp_pid(1)*sysTs^2+2*sysTs+4*temp_Kd_data(ii,1))*last_den_coeff)
                                    temp_N_data(ii,1) = %inf
                                    
                                elseif length_num_coeff == 3 & length_den_coeff == 3 & round((2+den_coeff(2))*(den_coeff(1)+1)) == round((den_coeff(1)-1)*(den_coeff(2))) then
                                // For PID type controller N ~= %inf
                                    temp_data = 2*(1-den_coeff(1))/(1+den_coeff(1))
                                    A = [sysTs*(2*temp_data-(2+temp_data)*num_coeff(2)) -4*(2+num_coeff(2));
                                         sysTs*((temp_data-2)-(temp_data+2)*num_coeff(1)) 4*(1-num_coeff(1))]
                                    b = [8+(4+2*temp_data)*num_coeff(2);(2*temp_data-4)+(2*temp_data+4)*num_coeff(1)]
                                    temp_pid = A\b
                                    temp_Ki_data(ii,1) = 1/temp_pid(1)
                                    temp_N_data(ii,1) = temp_pid(2)
                                    temp_Kd_data(ii,1) = temp_N_data(ii,1)*sysTs/temp_data
                                    temp_Kp_data(ii,1) = (4+2*temp_data)*last_num_coeff/(((4+2*temp_data)+4*temp_N_data(ii,1)+temp_pid(1)*sysTs*(2+temp_data))*last_den_coeff)
                                    
                                elseif (length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1) == 1 & den_coeff(1) == -1) |  (length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1) == -1) | (length_num_coeff == 2 & length_den_coeff == 2 & num_coeff(1) == -1 & den_coeff(1) == 1) then
                                // Parallel form of controller I type controller, D type of controller
                                     error(msprintf(gettext("The given control tranfer function is in Parallel PID form. Use pid function for validation.")))
                                        
                                 elseif length_num_coeff == 2 & length_den_coeff == 2 & den_coeff(1) < 0 & num_coeff(1) < 0 then
                                  // PD type controller N ~= %inf    
                                      temp_data = 2*(1+den_coeff(1))/(1-den_coeff(1))
                                      temp_N_data(ii,1) = -((2-temp_data)+(2+temp_data)*num_coeff(1))/(2*(1+num_coeff(1)))
                                      temp_Kd_data(ii,1) = temp_N_data(ii,1)*sysTs/temp_data
                                      temp_Kp_data(ii,1) = ((2+temp_data)*last_num_coeff)/((2+2*temp_N_data(ii,1)+temp_data)*last_den_coeff)
                                      temp_Ki_data(ii,1) = %inf
                                  
                                 elseif length_num_coeff == 2 & length_den_coeff == 2 then // & den_coeff(1) == 1      
                                  // PD type controller N == %inf    
                                      temp_Kd_data(ii,1) = sysTs*(1-num_coeff(1))/(2*(1+num_coeff(1)))
                                      temp_Kp_data(ii,1) = sysTs*last_num_coeff/((2*temp_Kd_data(ii,1)+sysTs)*last_den_coeff)
                                      temp_Ki_data(ii,1) = %inf
                                      temp_N_data(ii,1) = %inf
                                 else
                                     disp(systf_num(ii)/systf_den(ii))
                                     error(msprintf(gettext("Above equation is not a valid standard PID controller")))                          
                                end // end of when Irormula = Trapezoidal and dformula = trapezoidal                            
//************************************************************************************************************************************************************************
                            end// end of the Iformula and Dformula definder
                    end//end of for loop
        end
        output = varargin_data(1)
    end//end of extraction of the transfer function
    if current_length<=2 then
        //output_data = hypermat([current_size], pid_function(:,1))
        Kp = hypermat([current_size],temp_Kp_data(:,1))
        Ki = hypermat([current_size],temp_Ki_data(:,1))
        Kd = hypermat([current_size],temp_Kd_data(:,1))
        N = hypermat([current_size],temp_N_data(:,1))
    elseif current_length > 2 then
        //output_data = hypermat([current_size])
        //output_data(:,:,:,:)= pid_function
        Kp = hypermat([current_size])
        Kp(:,:,:,:) = temp_Kp_data
        Ki = hypermat([current_size])
        Ki(:,:,:,:) = temp_Ki_data
        Kd = hypermat([current_size])
        Kd(:,:,:,:) = temp_Kd_data
        N = hypermat([current_size])
        N(:,:,:,:) = temp_N_data
    end
//                        disp(Kp)
//                        disp(Ki)
//                        disp(Kd)
//                        disp(N)
        //impdata = 2
                impdata.Kp = Kp
                impdata.Ki = Ki
                impdata.Kd = Kd
                impdata.N = N
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
                TimeData = 0
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
                            elseif varargin(ii) == 'TimeUnit' | varargin(ii) == 'time' then
                                TimeData = ii
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
                disp(TimeData)
                if TimeData == 0 | TimeData+1 > rhs then
                    timeUnit = 'second'
                elseif TimeData ~= 0 then
                    timeUnitArray = ["nanoseconds" "nanosecond" "microseconds" "microsecond"  "milliseconds" "millisecond" "seconds" "second" "minutes" "minute" "hours" "hour" "days" "day" "weeks" "week" "months" "month" "years" "year" ]
                    findTimeUnit = find(timeUnitArray == varargin(TimeData+1))
                    //disp(findTimeUnit)
                    if size(findTimeUnit,"r") ~= 0 then
                        timeUnit = varargin(TimeData+1)
                    else
                        error(msprintf(gettext("specified time units is nanoseconds, microseconds, milliseconds, seconds, minutes, hours, days, weeks, months, years .")))
                    end
                    
                end
//------------------------------------------------------------------------------                
                extdata("         Kp") = Kp
                extdata("         Ki") = Ki
                extdata("         Kd") = Kd
                extdata("          N") = N   
                extdata("         Ts") = sysTs
                extdata("   TimeUnit") = timeUnit                
                extdata("   Iformula") = impdata.Iformula
                extdata("   Dformula") = impdata.Dformula
                extdata(" InputDelay") = impdata.InputDelay
                extdata("OutputDelay") = impdata.OutputDelay
                extdata("       Name") = impdata.Name
                extdata("      Notes") = impdata.Notes
                extdata("   UserData") = impdata.UserData
endfunction
