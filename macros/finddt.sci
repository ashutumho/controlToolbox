function output= finddt(h)
//It finds that the given transfer function is of continuous domain type.
//Calling Sequence
//output = finddt(tf)
//
//Parameters
//tf : Transfer function or hyper matrix of the transfer functions
//
//Description
//output : The output is Boolean type. It returns "T"(True) when the model or the hyper matrix of the model is discrete domain type and return "F"(False) if the model is continuous domain type. 
//
//Examples
//s=poly(0,'s');H1=(1+2*s)/s^2;S1bis=syslin('c',H1)
//output = finddt(S1bis)
//z=poly(0,'z'); h=syslin(0.1,(1-2*z)/(z^2+0.3*z+1))
//output = finddt(h)
// Authors
//
// Ashutosh Kumar Bhargava 
    [lhs,rhs]=argn(0)
    if typeof(h)<>'rational' then
        error(msprintf(gettext("input argument is not rational/ transfer function type")))
    end
    if h.dt=='c'| h.dt == 0 then
        output = %F
    else
        output = %T
    end
endfunction
