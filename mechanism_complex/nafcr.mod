: Fast Na+ channel

NEURON {
	THREADSAFE
	SUFFIX Nafcr
	USEION na READ ena WRITE ina
	RANGE gnafbar, ina, NAFCR_gna
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

:INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)} :deprecated!

PARAMETER {
	v (mV)
	dt (ms)
	gnafbar= 0.086 (mho/cm2) <0,1e9>
	:ena = 55 (mV) : WILL BE IGNORED AND SET BY NEURON
	ena (mV)
}

STATE {
	NAFCR_m NAFCR_h
}

ASSIGNED {
	ina (mA/cm2)
	NAFCR_minf NAFCR_hinf 
	NAFCR_mtau (ms)
	NAFCR_htau (ms)
	NAFCR_gna (mho/cm2)
}

INITIAL {
	rate(v)
	NAFCR_m = NAFCR_minf
	NAFCR_h = NAFCR_hinf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	NAFCR_gna = gnafbar*NAFCR_m*NAFCR_m*NAFCR_m*NAFCR_h
	ina = NAFCR_gna*(v-55)

}

DERIVATIVE states {
	rate(v)
	NAFCR_m' = (NAFCR_minf-NAFCR_m)/NAFCR_mtau
	NAFCR_h' = (NAFCR_hinf-NAFCR_h)/NAFCR_htau
}

UNITSOFF

FUNCTION NAFCR_malf( v){ LOCAL NAFCR_va 
	NAFCR_va=v+28
	if (fabs(NAFCR_va)<1e-04){
	   NAFCR_malf= -0.2816*(-9.3 + NAFCR_va*0.5)
	}else{
	   NAFCR_malf = -0.2816*(v+28)/(-1+exp(-(v+28)/9.3))
	}
}


FUNCTION NAFCR_mbet(v(mV))(/ms) { LOCAL NAFCR_vb 
	NAFCR_vb=v+1
	if (fabs(NAFCR_vb)<1e-04){
	   NAFCR_mbet = 0.2464*(6 + NAFCR_vb*0.5)
	}else{
	   NAFCR_mbet = 0.2464*(v+1)/(-1+exp((v+1)/6))
}
}	


FUNCTION NAFCR_half(v(mV))(/ms) { LOCAL NAFCR_vc 
	NAFCR_vc=v+43.1
	if (fabs(NAFCR_vc)<1e-04){
	   NAFCR_half=0.098*(20 + NAFCR_vc*0.5)
	}else{
	   NAFCR_half=0.098/exp((v+23.1)/20) :was 43.1
}
}


FUNCTION NAFCR_hbet(v(mV))(/ms) { LOCAL NAFCR_vd
	NAFCR_vd=v+13.1
	if (fabs(NAFCR_vd)<1e-04){
	   NAFCR_hbet=1.4*(10 + NAFCR_vd*0.5)
	}else{
	   NAFCR_hbet=1.4/(1+exp(-(v+25.1)/10)) :was 13.1 changed july 30, 2007
} 
}


PROCEDURE rate(v (mV)) {LOCAL NAFCR_msum, NAFCR_hsum, NAFCR_ma, NAFCR_mb, NAFCR_ha, NAFCR_hb
	
	NAFCR_ma=NAFCR_malf(v) NAFCR_mb=NAFCR_mbet(v) NAFCR_ha=NAFCR_half(v) NAFCR_hb=NAFCR_hbet(v)
	
	NAFCR_msum = NAFCR_ma+NAFCR_mb
	NAFCR_minf = NAFCR_ma/NAFCR_msum
	NAFCR_mtau = 1/(NAFCR_msum)
		
	NAFCR_hsum=NAFCR_ha+NAFCR_hb
	NAFCR_hinf=NAFCR_ha/NAFCR_hsum
	NAFCR_htau = 1 / (NAFCR_hsum)	
}

	
UNITSON


