mode(1)
//
// Demo of chgTimeUnit.sci
//

s = poly(0,'s'); z = poly(0,'z');
sys = syslin('c',(s+2)/(s^2+3*s+1))
systime = chgTimeUnit(sys,'hour')
systime1 = chgTimeUnit(sys,'nanosecond')
halt()   // Press return to continue
 
//========= E N D === O F === D E M O =========//
