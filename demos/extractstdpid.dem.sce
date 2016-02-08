mode(1)
//
// Demo of extractstdpid.sci
//

s = poly(0,'s'); sys = syslin('c',3*(s+1)*(s+2)/s);
output = extractstdpid(sys)
[Kp,Ki,Kd,Tf,Ts] = extractstdpid(sys)
z = poly(0,'z'); sys = syslin(0.1,(7*z-6.9)/(4*z-3.9));
[Kp,Ki,Kd,Tf,Ts] = extractstdpid(sys)
sys = stdpid(rand(2,3,4),2,4,5);
[Kp,Ki,Kd,Tf,Ts] = extractstdpid(sys,[1 3 2])
[Kp,Ki,Kd,Tf,Ts] = extractstdpid(sys,5)
//========= E N D === O F === D E M O =========//
