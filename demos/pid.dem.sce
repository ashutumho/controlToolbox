mode(1)
//
// Demo of pid.sci
//

output = pid(2)
output = pid(2,3)
output = pid(2,3,4)
output = pid(2,3,4,5)
output = pid(2,0,4)
output = pid(2,0,4,5)
output = pid(2,0,0,0,0.1)
output = pid(2,3,0,0,0.1,'Iformula','B')
output = pid(2,3,4,0,0.1,'Iformula','B','Dformula','T')
output = pid(2,3,4,5,0.1,'Iformula','T','Dformula','B')
s = poly(0,'s'); sys = syslin('c',3*(s+1)*(s+2)/s);
[output impdata]= pid(sys)
z = poly(0,'z'); sys = syslin(0.1,(19.146667 - 38.793333*z + 19.66*z^2)/(  1.6 - 3.4*z + 1.8*z^2))
[output impdata] = pid(sys,'Iformula','T','Dformula','B')
//========= E N D === O F === D E M O =========//
