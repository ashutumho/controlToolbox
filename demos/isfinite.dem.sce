mode(1)
//
// Demo of isfinite.sci
//

mat = [1 2 %inf 0 -%inf]; bool = isfinite(mat)
s = poly(0,'s') ; sys = syslin('c',1/(s+1)); bool = isfinite(sys)
sys = ssrand(2,3,8,list('st',2,5,5)); bool = isfinite(sys)
bool = isfinite(sys,'elem')
halt()   // Press return to continue
 
//========= E N D === O F === D E M O =========//
