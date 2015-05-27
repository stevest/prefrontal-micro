: Fast Na+ channel

NEURON {
	THREADSAFE
	SUFFIX Naf
	USEION na READ ena WRITE ina
	RANGE gnafbar, ina, gna
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
	v (mV)
	dt (ms)
	gnafbar= 0.086 (mho/cm2) <0,1e9>
	:ena = 55 (mV) : WILL BE IGNORED AND SET BY NEURON

}
STATE {
	m h
}
ASSIGNED {
	ina (mA/cm2)
	minf hinf 
	mtau (ms)
	htau (ms)
	gna (mho/cm2)
	
	ena	
}



INITIAL {
	rate(v)
	m = minf
	h = hinf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gna = gnafbar*m*m*m*h
	ina = gna*(v-55)
VERBATIM
	//if( ( t - (500 * (int)(t/500))) == 0 ) {	
	if( (t > 0.0)  ) {	
		if(   (int)(t*10)%1000 == 0 ) {	
			//printf("MODULO IS: %10.10f\n", t - (1.0 * (int)(t/1.0)) );
			//printf("mod is equal with zero float: %d\n", ( t - (1.0 * (int)(t/1.0))) == 0 );
			//printf("mod is equal with zero : %d\n", ( t - (1 * (int)(t/1))) == 0 );
			//printf("mod is equal with : %f\n",   (t*10)   );
			//printf("mod is equal with : %d\n",   (int)(t*10)   );
			//printf("mod is equal with : %f\n",   (((int)(t*10))/10)   );
			//printf("mod is equal with float: %d\n",  (int)(t*10)%100  );
			printf("@t = %10.10f\n",t);
		}
	}
ENDVERBATIM
}

DERIVATIVE states {
	rate(v)
	m' = (minf-m)/mtau
	h' = (hinf-h)/htau
}

UNITSOFF

FUNCTION malf( v){ :LOCAL va 
	:va=(v+28)
	if (fabs((v+28))<1e-04){
	   malf= -0.2816*(-9.3 + (v+28)*0.5)
	}else{
	   malf = -0.2816*(v+28)/(-1+exp(-(v+28)/9.3))
	}
}


FUNCTION mbet(v(mV))(/ms) { :LOCAL vb 
	:vb=(v+1)
	if (fabs((v+1))<1e-04){
	   mbet = 0.2464*(6 + (v+1)*0.5)
	}else{
	   mbet = 0.2464*(v+1)/(-1+exp((v+1)/6))
}
}	


FUNCTION half(v(mV))(/ms) { :LOCAL vc 
	:vc=(v+43.1)
	if (fabs((v+43.1))<1e-04){
	   half=0.098*(20 + (v+43.1)*0.5)
	}else{
	   half=0.098/exp((v+23.1)/20) :was 43.1
}
}


FUNCTION hbet(v(mV))(/ms) { :LOCAL vd
	:vd=(v+13.1)
	if (fabs((v+13.1))<1e-04){
	   hbet=1.4*(10 + (v+13.1)*0.5)
	}else{
	   hbet=1.4/(1+exp(-(v+25.1)/10)) :was 13.1 changed july 30, 2007
} 
}




PROCEDURE rate(v (mV)) { :LOCAL q10, msum, hsum, ma, mb, ha, hb
	

	:ma=malf(v)
	:mb=mbet(v)
	:ha=half(v)
	:hb=hbet(v)
	
	:msum = malf(v)+mbet(v)
	minf = malf(v)/(malf(v)+mbet(v))
	mtau = 1/(malf(v)+mbet(v))
	
	
	:hsum=half(v)+hbet(v)
	hinf=half(v)/(half(v)+hbet(v))
	htau = 1 / (half(v)+hbet(v))
	
}

	
UNITSON


