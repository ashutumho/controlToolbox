mode(1)
//
// Demo of isproper.sci
//

s = poly(0,'s'); sys = syslin('c',1/(s+1));bool = isproper(sys)
sys = pid(rand(2,2),0,3);
bool = isproper(sys)
bool = isproper(sys,'elem')
sys =  ssrand(2,3,8,list('st',2,5,5)); bool = isproper(sys)
halt()   // Press return to continue
 
//========= E N D === O F === D E M O =========//
