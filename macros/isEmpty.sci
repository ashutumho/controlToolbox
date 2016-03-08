
function output =isEmpty(varargin)
//Determines that given state space, transfer function, matrix, array and hyper-matrix is empty
//
//Calling Sequence
//bool = isempty(sys)
//bool = isempty(mat)
//
//Parameters
//sys : represents the valid form of state space model or transfer function
//mat : a matrix of double, represents the array, matrix or hyper-matrix
//bool : a boolean, represents the output    
//Description 
//bool = isempty(sys) , the output will be true(T) when there is no input, output, or both of a system. The output is false(F) for the transfer functions and state space equation with proper input/output.
//bool = isempty(mat), the output will be true(T) when the matrix is null otherwise output will be true.
//Examples
//a = rand(2,2,2); bool = isempty(a);a(:,:,:) = [];bool = isempty(a)
//a = [1 2 3]; bool = isempty(a);a =[]; bool = isempty(a)
//A=diag([1,2,3]);B=[2;2;2];C=[1 2 3];D=0; sys= syslin('c',A,B,C,D); bool = isempty(sys)
//A=diag([1,2,3]);B=[2;2;2];C=[];D=0; sys= syslin('c',A,B,C,D); bool = isempty(sys)
//A=diag([1,2,3]);B=[2;2;2];C=[];D=[]; sys= syslin('c',A,B,C,D); bool = isempty(sys)
//A=diag([1,2,3]);B=[];C=[1 2 3];D=0; sys= syslin('c',A,B,C,D); bool = isempty(sys)
//A=diag([1,2,3]);B=[];C=[1 2 3];D=[]; sys= syslin('c',A,B,C,D); bool = isempty(sys)
//A=diag([1,2,3]);B=[];C=[];D=0; sys= syslin('c',A,B,C,D); bool = isempty(sys)
//A=diag([1,2,3]);B=[];C=[];D=[]; sys= syslin('c',A,B,C,D); bool = isempty(sys)
//Author
//Ashutosh Kumar Bhargava

    [lhs,rhs]=argn(0)
    // check the input elements
    if rhs == 0 then
        error(msprintf(gettext("Function has no input argument..")))
    elseif rhs >=2 then
        error(msprintf(gettext("Incorrect number of input arguments.")))
    end
    // storing the data
    sysData = varargin($)
    select typeof(sysData)
    case "constant" then
    case "hypermat" then
    case "st" then
    case "state-space" then
    case "rational" then
        sysData = tf2ss(sysData)
    else
        error(msprintf(gettext("Incompatible input argument.")))
    end
    //data analysis
    if typeof(sysData) == 'rational' | typeof(sysData) == 'state-space' then
        [A1,B1,C1,D1]=abcd(sysData)
        if ((size(C1,"r")==0)&(size(D1,"r")==0)) | ((size(B1,"r")==0)&(size(D1,"r")==0)) | ((size(B1,"r")==0)&(size(C1,"r")==0)&(size(D1,"r")==0)) then
            output = %T
        else
            output = %F
        end
    elseif typeof(sysData) == 'constant' | typeof(sysData) == 'hypermat' then
        if size(sysData,"r") == 0 then
            output = %T
        else
            output = %F    
        end
    elseif typeof(sysData) == 'st' then
        if size(sysData.freq,"r") == 0 then
            output = %T
        else
            output = %F
        end
    end

endfunction
