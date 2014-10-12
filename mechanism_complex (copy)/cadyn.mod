: simple first-order model of calcium dynamics

NEURON {
	THREADSAFE
        SUFFIX cadyn
        USEION ca READ cai,ica WRITE cai 
        RANGE ca :Nassi
	GLOBAL depth,cainf,taur :Nassi
}

UNITS {
        (molar) = (1/liter)		:Nassi  moles do not appear in units
        (mM) = (milli/liter)
	(um)	= (micron) :Nassi
        (mA) = (milliamp)
	(msM)	= (ms mM)  :Nassi
        FARADAY    = (faraday) (coul)
}

PARAMETER {
       depth	= .1	(um)		: depth of shell Nassi
        taur =  200 (ms)	: rate of calcium removal for stress conditions
	cainf	= 50e-6(mM)	:changed oct2
	cai		(mM)
}

ASSIGNED {
	ica		(mA/cm2)
	drive_channel	(mM/ms)
}


STATE {
	ca		(mM) 
}

INITIAL {
	ca = cainf
}

 
BREAKPOINT { 
	SOLVE state METHOD cnexp :EULER IS NOT THREAD SAFE!
}

DERIVATIVE state { 

	drive_channel =  - (10000) * ica / (2 * FARADAY * depth)
	if (drive_channel <= 0.) { drive_channel = 0.  }   : cannot pump inward 
         
        ca' = drive_channel/18 + (cainf -ca)/taur*11  :*11    :Kiki had it(7)
	cai = ca
}


INITIAL {
	ca = cainf
}
