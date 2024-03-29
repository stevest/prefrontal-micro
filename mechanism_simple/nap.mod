: Persistent Na+ channel

NEURON {
	THREADSAFE
	SUFFIX nap
	USEION na READ ena WRITE ina
	RANGE gnabar, ina, gna
	RANGE DA_alphamshift,DA_betamshift
	RANGE DA_alphahfactor, DA_betahfactor
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
	v (mV)
	dt (ms)
	gnabar= 0.0022 (mho/cm2) <0,1e9>
	:ena = 55 (mV) : WILL BE IGNORED AND SET BY NEURON
	DA_alphamshift=0 : 2 for 100% DA, 0 otherwise
	DA_betamshift=0  : 5 for 100% DA,0 otherwise
	DA_alphahfactor=0: -.8e-5 for DA, 0 otherwise
	DA_betahfactor=0 : 0.014286-0.02 for DA, 0 otherwise
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
	gna = gnabar*m*h
	ina = gna*(v-55)
	
}

DERIVATIVE states {
	rate(v)
	m' = (minf-m)/mtau
	h' = (hinf-h)/htau
}

UNITSOFF

FUNCTION malf( v){ :LOCAL va 
	:va=v+12+DA_alphamshift
	if (fabs(v+12+DA_alphamshift)<1e-04){
		malf = (-0.2816*(v+12+DA_alphamshift+ 0.00001))/(-1+exp(-(v+12+DA_alphamshift+ 0.00001)/9.3))
	} else {
		malf = (-0.2816*(v+12+DA_alphamshift))/(-1+exp(-(v+12+DA_alphamshift)/9.3))
	}
	
}


FUNCTION mbet(v(mV))(/ms) { :LOCAL vb 
	:vb=v-15+DA_betamshift
	if (fabs(v-15+DA_betamshift)<1e-04){
		mbet = (0.2464*(v-15+DA_betamshift+ 0.00001))/(-1+exp((v-15+DA_betamshift+ 0.00001)/6))
	} else {
		mbet = (0.2464*(v-15+DA_betamshift))/(-1+exp((v-15+DA_betamshift)/6))
	}

	

}	


FUNCTION half(v(mV))(/ms) { :LOCAL vc 
	:vc=v+42.8477
	if (fabs(v+42.8477)<1e-04){
		half= (2.8e-5+DA_alphahfactor)*(exp(-(v+42.8477+ 0.00001)/4.0248))
	} else {
	        half= (2.8e-5+DA_alphahfactor)*(exp(-(v+42.8477)/4.0248))
	}

}


FUNCTION hbet(v(mV))(/ms) { :LOCAL vd
	:vd=v-413.9284
	if (fabs(v-413.9284)<1e-04){
	        hbet= (0.02+DA_betahfactor)/(1+exp(-(v-413.9284+ 0.00001)/148.2589))
	} else {
		hbet= (0.02+DA_betahfactor)/(1+exp(-(v-413.9284)/148.2589))
	}
 
}




PROCEDURE rate(v (mV)) { :LOCAL msum, hsum, ma, mb, ha, hb
	:ma=malf(v)
	:mb=mbet(v)
	:ha=half(v)
	:hb=hbet(v)
	
	:msum = malf(v)+mbet(v)
	minf = malf(v)/(malf(v)+mbet(v))
	mtau = 1/(malf(v)+mbet(v))
	
	
	:hsum = half(v)+hbet(v)
	hinf = half(v)/(half(v)+hbet(v))
	htau = 1/(half(v)+hbet(v))
}

	
UNITSON




