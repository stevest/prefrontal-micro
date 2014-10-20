: Fast Na+ channel
: added the 's' attenuation system from hha2.mod
: Kiki Sidiropoulou
: September 27, 2007

NEURON {
	SUFFIX Nafx
	USEION na READ ena WRITE ina
	RANGE gnafbar, ina, gna, ar2
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
	v (mV)
	dt (ms)
	gnafbar	= 0 (mho/cm2)
	ena = 55 (mV)
	
	:PARAMETERS FOR S ATTENUATION SYSTEM
	NAFX_taumin = 30 (ms)  :min activation time for "s" attenuation system
        NAFX_vhalfr =-60 (mV)       :half potential for "s" attenuation system, -60
        NAFX_vvh = -58		(mV) 
 	NAFX_vvs = 2 (mV)
	NAFX_a0r = 0.0003 (/ms)
        NAFX_b0r = 0.0003 (/ms)
        NAFX_zetar = 12    
	NAFX_zetas = 12   
        NAFX_gmr = 0.2   
	ar2 = 1.0               :initialized parameter for location-dependent
                                :Na-conductance attenuation, "s", (ar=1 -> zero attenuation)
}
STATE {
	m h s
}
ASSIGNED {
	celsius (degC)
	ina (mA/cm2)
	NAFX_minf 
	NAFX_hinf
	NAFX_sinf 
	NAFX_mtau (ms)
	NAFX_htau (ms)
	NAFX_stau (ms)
	gna (mho/cm2)
	
}



INITIAL {
	rate(v, ar2)
	m = NAFX_minf
	h = NAFX_hinf
	s = NAFX_sinf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gna = gnafbar*m*m*m*h*s
	ina = gna*(v-55)
	
}

DERIVATIVE states {
	rate(v, ar2)
	m' = (NAFX_minf-m)/NAFX_mtau
	h' = (NAFX_hinf-h)/NAFX_htau
	s' = (NAFX_sinf-s)/NAFX_stau
}

UNITSOFF

FUNCTION NAFX_malf( v){ LOCAL va 
	va=v+28
	:va=v+28
	if (fabs(va)<1e-04){
	   NAFX_malf= -0.2816*(-9.3 + va*0.5)
	   :malf= -0.2816*(-9.3 + va*0.5)
	}else{
	   NAFX_malf = -0.2816*(v+28)/(-1+exp(-va/9.3))
	}
}


FUNCTION NAFX_mbet(v(mV))(/ms) { LOCAL vb 
	vb=v+1
	:vb=v+1
	if (fabs(vb)<1e-04){
	    NAFX_mbet = 0.2464*(6+vb*0.5)
	    :mbet = 0.2464*(6 + vb*0.5)
	}else{
	   NAFX_mbet = 0.2464*(v+1)/(-1+exp(vb/6))	  :/(-1+exp((v+1)/6))
	}
	}	


FUNCTION NAFX_half(v(mV))(/ms) { LOCAL vc 
	:vc=v+15.1
	vc=v+40.1	:changed to 40.1 by kiki
	if (fabs(vc)<1e-04){
	   NAFX_half=0.098*(20 + vc*0.5)
	}else{
	   NAFX_half=0.098/exp(vc+43.1/20)  :43.1, also spike train attenuation
}
}


FUNCTION NAFX_hbet(v(mV))(/ms) { LOCAL vd
	:vd=v+13.1
	vd=v+13.1  :decreasing it increases the peak current
	if (fabs(vd)<1e-04){
	   NAFX_hbet=1.4*(10 + vd*0.5)
	}else{
	   NAFX_hbet=1.4/(1+exp(-(vd-13.1)/10))  :13.1 increasing it, increases the spike train attenuation and increases spike width
} 
}


:FUNCTIONS FOR S 
FUNCTION NAFX_alpv(v(mV)) {
         NAFX_alpv = 1/(1+exp((v-NAFX_vvh)/NAFX_vvs))
}


FUNCTION NAFX_alpr(v(mV)) {       :used in "s" activation system tau

  NAFX_alpr = exp(1.e-3*NAFX_zetar*(v-NAFX_vhalfr)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION NAFX_betr(v(mV)) {       :used in "s" activation system tau

  NAFX_betr = exp(1.e-3*NAFX_zetar*NAFX_gmr*(v-NAFX_vhalfr)*9.648e4/(8.315*(273.16+celsius))) 
}



PROCEDURE rate(v (mV),NAFX_a2) {LOCAL NAFX_q10, NAFX_msum, NAFX_hsum, NAFX_ma, NAFX_mb, NAFX_ha, NAFX_hb,NAFX_c
	

	NAFX_ma=NAFX_malf(v) NAFX_mb=NAFX_mbet(v) NAFX_ha=NAFX_half(v) NAFX_hb=NAFX_hbet(v)
	
	NAFX_msum = NAFX_ma+NAFX_mb
	NAFX_minf = NAFX_ma/NAFX_msum
	NAFX_mtau = 1/(NAFX_msum)
	
	
	NAFX_hsum=NAFX_ha+NAFX_hb
	NAFX_hinf=NAFX_ha/NAFX_hsum
	NAFX_htau = 1 / (NAFX_hsum)

	NAFX_stau = NAFX_betr(v)/(NAFX_a0r*(1+NAFX_alpr(v))) 
	if (NAFX_stau<NAFX_taumin) {NAFX_stau=NAFX_taumin} :s activation tau
	NAFX_c = NAFX_alpv(v)
	NAFX_sinf = NAFX_c+NAFX_a2*(1-NAFX_c) 	
}

	
UNITSON


