: simple first-order model of calcium dynamics

NEURON {
        SUFFIX cadyn
        USEION ca READ cai,ica WRITE cai 
        RANGE ca :Nassi
	GLOBAL depth,cainf,taur :Nassi
         :RANGE CAF, tca, cai (Nassi)
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
        :tca= 70 (ms)           : decay time constant
	:tca= 300 (ms)           : decay time constant (new) Sept 6, 2007 Nassi
       : cainf= 100e-6   (mM)      : (50e-6)equilibrium ca2+ concentration Nassi
        :dep= 2e-4 (micron)     : depth of shell for ca2+ diffusion Nassi
        taur =  200 (ms)	: rate of calcium removal for stress conditions
	cainf	= 50e-6(mM)	:changed oct2
	cai		(mM)
}

ASSIGNED {
	ica		(mA/cm2)
	drive_channel	(mM/ms)
}
: ASSIGNED {
:        ica     (mA/cm2)
 :       diam    (micron)
  :      A       (/coul/cm)  
   :     CAF     ()
: } Nassi

STATE {
	ca		(mM) 
}
: STATE { cai (mM) } Nassi


: BREAKPOINT { 
 :       SOLVE states METHOD cnexp
: }

 : DERIVATIVE states {
         :cai'= -A*ica  + (cai-cainf)/tca*7
	: cai'= -A*CAF*ica - (cai-cainf)/tca :old Nassi
: } Nassi

: INITIAL {
 :       A =(1e4)/(F*dep)
: } Nassi
 
BREAKPOINT {
	SOLVE state METHOD euler
}

DERIVATIVE state { 

	drive_channel =  - (10000) * ica / (2 * FARADAY * depth)
	if (drive_channel <= 0.) { drive_channel = 0.  }   : cannot pump inward 
         
	:ca' = drive_channel + (cainf-ca)/taur
        ca' = drive_channel/18 + (cainf -ca)/taur*11:*11    :(7)
	cai = ca
}


INITIAL {
	ca = cainf
}
