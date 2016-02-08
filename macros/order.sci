
function [output]= order(h)
//finds model order
//Calling Sequence
//output = order(sys)
//
//Parameters
//sys : transfer function / state space model
//
//Description
//output = order(sys) returns the order of sys i.e. if sys is a transfer function then it will give the  number of poles (all the pole must define at the left side of the imaginary axis) and if sys is state space model type then it will return order as the number of the stats in model.  
//
//Examples
//s=poly(0,'s');H1=(1+2*s)/s^2;S1bis=syslin('c',H1)
//output = order(S1bis)
//z=poly(0,'z'); h=syslin(0.1,(1-2*z)/(z^2+0.3*z+1))
//output = order(h)
// Authors
//
// Ashutosh Kumar Bhargava 
    [lhs,rhs]=argn(0)
    //disp(rhs)
    if rhs ~= 1 then
        error(msprintf("Wrong size for input argument"))
    end
// find the data type
    datatype = typeof(h)
    if datatype == 'rational' then                  //if transfred data is transfer function type
        output = length(coeff(h.den))-1             //it will show the order of the transfer function
    elseif datatype == 'state-space' 
        dim = size(h("A"))
        output = dim(1)                            // it will tranfer the number of states of the state-space equation
    else
         error(msprintf("Wrong data type"))
    end
   // output = 1
endfunction
