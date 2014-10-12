TITLE  H-current that uses Na ions
: Updated to use Cvode by Yiota Poirazi 12/1/2005

NEURON {
	SUFFIX hin
        RANGE  gbar,HIN_vhalf, K, HIN_taun, HIN_ninf, g, ihi   
	USEION hi READ ehi WRITE ihi VALENCE 1    
	THREADSAFE  
}

UNITS {
	(um) = (micrometer)
	(mA) = (milliamp)
	(uA) = (microamp)
	(mV) = (millivolt)
	(pmho) = (picomho)
	(mmho) = (millimho)
}

:INDEPENDENT {t FROM 0 TO 1 WITH 100 (ms)}

PARAMETER {              : parameters that can be entered when function is called in cell-setup
        ena    = 55    (mV)
        ehi     = -10   (mV)
	K      = 10.0   (mV)	:8.5
	gbar   = 0     (mho/cm2)  : initialize conductance to zero
	HIN_vhalf  = -90   (mV)       : half potential
}	


STATE {                : the unknown parameters to be solved in the DEs
	n
}

ASSIGNED {             : parameters needed to solve DE
        v 
	ihi (mA/cm2)
	HIN_ninf
	HIN_taun (ms)
	g
}

        


INITIAL {               : initialize the following parameter using states()
	rates()	
	n = HIN_ninf
	g = gbar*n
	ihi = g*(v-ehi)
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar*n  
	ihi = g*(v-ehi)  
}

DERIVATIVE states {
	rates()
        n' = (HIN_ninf - n)/HIN_taun
}

PROCEDURE rates() {  
 
 	if (v > -10) {
	   HIN_taun = 1
	} else {
           HIN_taun = 2*(1/(exp((v+145)/-17.5)+exp((v+16.8)/16.5)) + 10) :h activation tau +5

	}  
         HIN_ninf = 1 - (1 / (1 + exp((HIN_vhalf - v)/K)))                  :steady state value
}



