mode(1)
//
// Demo of chgFreqUnit.sci
//

sys = frd(1:10,1:10)
sysfreq = chgFreqUnit(sys,'Khz')
sysfreq1 = chgFreqUnit(sys,'rpm')
halt()   // Press return to continue
 
//========= E N D === O F === D E M O =========//
