mode(1)
//
// Demo of isEmpty.sci
//

a = rand(2,2,2); bool = isempty(a);a(:,:,:) = [];bool = isempty(a)
a = [1 2 3]; bool = isempty(a);a =[]; bool = isempty(a)
A=diag([1,2,3]);B=[2;2;2];C=[1 2 3];D=0; sys= syslin('c',A,B,C,D); bool = isempty(sys)
A=diag([1,2,3]);B=[2;2;2];C=[];D=0; sys= syslin('c',A,B,C,D); bool = isempty(sys)
A=diag([1,2,3]);B=[2;2;2];C=[];D=[]; sys= syslin('c',A,B,C,D); bool = isempty(sys)
A=diag([1,2,3]);B=[];C=[1 2 3];D=0; sys= syslin('c',A,B,C,D); bool = isempty(sys)
A=diag([1,2,3]);B=[];C=[1 2 3];D=[]; sys= syslin('c',A,B,C,D); bool = isempty(sys)
A=diag([1,2,3]);B=[];C=[];D=0; sys= syslin('c',A,B,C,D); bool = isempty(sys)
A=diag([1,2,3]);B=[];C=[];D=[]; sys= syslin('c',A,B,C,D); bool = isempty(sys)
Author
Ashutosh Kumar Bhargava
//========= E N D === O F === D E M O =========//
