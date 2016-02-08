// gain cross over frequency
// given gain cross over frequency=gg_cof()
//funcprot(0)
function[frq]= gg_cof(h,gain)
    [lhs,rhs]=argn(0);
    
//    if argn(2) < 1 then
//        error(msprintf(_("%s: Wrong number of input argument(s): %d expected.\n"),"g_cof",1));
//    end

    select typeof(h)
    case "rational" then
    case "state-space" then
        h=ss2tf(h)
    else
        msprintf("\n")
        error(97,1),
    end;
    if or(size(h)<>[1 1]) then
        error(msprintf(_("\n %s: Wrong size for input argument #%d: Single input, single output system expected.\n"),"g_cof",1))
    end
    //
    eps=1.e-7;// threshold used for testing if complex numbers are real or pure imaginary
        if h.dt=="c" then  //continuous time case
        w=poly(0,"w");
        niw = horner(h.num,%i*w)*conj(horner(h.num,%i*w));
        diw = horner(h.den,%i*w)*conj(horner(h.den,%i*w))
        
         w=roots(niw-gain^2*diw,"e");
         ws=real(w(find((abs(imag(w))<eps)&(real(w)>0)))); //frequency points with unitary modulus
         if ws ==[]  then
             frq = [];
             return;
         end
         frq = ws;
     else 
        if h.dt=="d" then
            dt=1;
        else
            dt=h.dt;
        end
        // |h(e^(i*w*dt))|=1 <-- h(e^(i*w*dt))*h(e^(-i*w*dt))
        z=poly(0,varn(h.den));
        sm=simp_mode();
        simp_mode(%f);
        hh=h*horner(h,1/z)-gain^2;
        simp_mode(sm);
        //find the numerator roots
        z=roots(hh.num,"e");
        z(abs(abs(z)-1)>eps)=[];// retain only roots with modulus equal to 1
        w=log(z)/(%i*dt);
        ws=real(w(abs(imag(w))<eps&real(w)>0)); //frequency points with unitary modulus
        if ws==[] then
            frq=[];
            return
        end
        frq = ws;
        end
        
endfunction

