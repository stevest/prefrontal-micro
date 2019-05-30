COMMENT

//****************************//

// Created by Alon Polsky 	//

//    apmega@yahoo.com		//

//		2007			//

//****************************//

ENDCOMMENT



TITLE NMDA synapse 



NEURON {
	THREADSAFE
	POINT_PROCESS nmda_nojump
	USEION ca READ cai WRITE ica VALENCE 2
	NONSPECIFIC_CURRENT inmda 
	RANGE e ,gmax,inmda
	RANGE gnmda, tau1, NLalpha
	RANGE srcgid, sid, n, cellid, myflag
	GLOBAL gama,tau2
}



UNITS {

	(nA) 	= (nanoamp)

	(mV)	= (millivolt)

	(nS) 	= (nanomho)

	(mM)    = (milli/liter)

	:F	= 96480 (coul)

	:R       = 8.314 (volt-coul/degC)

	FARADAY = (faraday) (coulomb)

	R	= (k-mole) (joule/degC)



}



PARAMETER {

	gmax=1	(nS)

	e= 0.0	(mV)

	tau1=30	(ms): Was 90	

	tau2=5	(ms)	

	gama=0.08 :0.08 	(/mV) :Was 0.08

	dt (ms)

	v		(mV)
	n=0.25		(/mM) :was 0.25
	mybeta = 0.0
	cellid = 0
	srcgid = -1
	myflag = -1
	sid = -1
	NLalpha = 0

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
	: gnmda=(A-B)/(1+n*exp(-gama*v - mybeta))
	gnmda=(A-B)
	inmda = (1e-3) * gnmda * (v-e)


	ica = inmda * 0.01
	NLalpha = A
}



DERIVATIVE state {

	

	A'=-A/tau1

	B'=-B/tau2

	

}



NET_RECEIVE (weight) {
	gmax=weight
	A = A + gmax
	B = B + gmax

VERBATIM
//		printf("@NMDA=%f cellid=%g srcgid=%g sid=%g\n",t,cellid, srcgid, sid );
ENDVERBATIM
}




