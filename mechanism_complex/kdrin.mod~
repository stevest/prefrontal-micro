: Delayed rectifier K+ channel

NEURON {
	SUFFIX kdrin
	USEION k READ ki, ko WRITE ik
	RANGE gkdrbar, ik, gk
	
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
PARAMETER {
	v (mV)
	dt (ms)
	gkdrbar= 0.0338 (mho/cm2) <0,1e9>
	
	
}

STATE {
	n
}

ASSIGNED {
	ik (mA/cm2)
	KDRIN_inf
	KDRIN_tau (ms)
	gk (mho/cm2)
	ek (mV)
	ki (mM)
	ko (mM)

}


INITIAL {
	rate(v)
	n = KDRIN_inf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gk= gkdrbar*n*n*n*n
	ek = 25 * log(ko/ki)
	ik = gk*(v-ek)
	
}

DERIVATIVE states {
	rate(v)
	n' = (KDRIN_inf-n)/KDRIN_tau
}

UNITSOFF

FUNCTION KDRIN_alf(v){ LOCAL KDRIN_va 
	
	   KDRIN_va=v-13
	if (fabs(KDRIN_va)<1e-04){
	   KDRIN_va=KDRIN_va+0.0001
		KDRIN_alf= (-0.018*KDRIN_va)/(-1+exp(-(KDRIN_va/25)))
	} else {
	  	KDRIN_alf = (-0.018*(v-13))/(-1+exp(-((v-13)/25)))
	}
}


FUNCTION KDRIN_bet(v) { LOCAL KDRIN_vb 
	
	  KDRIN_vb=v-23
	if (fabs(KDRIN_vb)<1e-04){
	  KDRIN_vb=KDRIN_vb+0.0001
		KDRIN_bet= (0.0054*KDRIN_vb)/(-1+exp(KDRIN_vb/12))
	} else {
	  	KDRIN_bet = (0.0054*(v-23))/(-1+exp((v-23)/12))
	}
}	






PROCEDURE rate(v (mV)) {LOCAL q10, KDRIN_sum, KDRIN_aa, KDRIN_ab
	
	KDRIN_aa=KDRIN_alf(v) KDRIN_ab=KDRIN_bet(v) 
	
	KDRIN_sum = KDRIN_aa+KDRIN_ab
	KDRIN_inf = KDRIN_aa/KDRIN_sum
	KDRIN_tau = 1/(KDRIN_sum)
	
	
}

UNITSON	



