COMMENT
//****************************//
// Created by Alon Polsky 	//
//    apmega@yahoo.com		//
//		2007			//
//****************************//
ENDCOMMENT

TITLE NMDA synapse 

NEURON {
	POINT_PROCESS nmda_segev
	USEION ca READ cai WRITE ica VALENCE 2
	NONSPECIFIC_CURRENT inmda 

	RANGE e ,gmax,inmda,tau1,tau2, gama,n

	RANGE gnmda

	GLOBAL n, gama, mg

}

UNITS {
	(nA) 	= (nanoamp)
	(mV)	= (millivolt)
	(nS) 	= (nanomho)
	(mM)    = (milli/liter)
        F	= 96480 (coul)
        R       = 8.314 (volt-coul/degC)

}

PARAMETER {
	gmax=1	(nS)
	eNMDA= 	0.0 (mV)	:0.0
	tau1=90	(ms)	: Was 90	
	tau2=5	(ms)	:was 5
	gama= 0.08	(/mV) :Scale Was 0.08 
	dt (ms)
	v		(mV)
	n=0.25	(/mM) :Translation was 0.25
	mg      = 1      (mM)           : external magnesium concentration
}

ASSIGNED { 
	inmda		(nA)   
	gnmda		(nS)
	ica 		(mA/cm2)
	cai		(mM)	
}
STATE {
	A (nS)
	B (nS)
}

INITIAL {
      gnmda=0 
	A=0
	B=0
}

BREAKPOINT {
	

	SOLVE state METHOD cnexp

	gnmda=(A-B)/(1+n*exp(-gama*v)) :  * mgblock(v)
	:inmda = (1e-3) * gnmda * (v-e)
	inmda = (1e-3) * gnmda * (v-eNMDA)
	ica = inmda/10
	:printf("t1: %f, t2: %f\n", tau1, tau2)



}


:FUNCTION mgblock(v(mV)) {
:        TABLE 
:        DEPEND mg
:        FROM -140 TO 80 WITH 1000
:
        : from Jahr & Stevens
:
:      
	 :mgblock = 1 / (1 + exp(0.072 (/mV) * -v) * (mg / 3.57 (mM)))  :was 0.062, changed to 0.072 to get a better voltage-dependence of NMDA currents, july 2008, kiki
:	
:}


DERIVATIVE state {
	
	A'=-A/tau1
	B'=-B/tau2
	
}

NET_RECEIVE (weight) {
	gmax=weight
	state_discontinuity( A, A+ gmax)
	state_discontinuity( B, B+ gmax)

}

