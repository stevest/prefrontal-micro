TITLE kinetic NMDA receptor model

COMMENT
-----------------------------------------------------------------------------

    Kinetic model of NMDA receptors
    ===============================

    10-state gating model:
    Kampa et al. (2004) J Physiol
  
      U -- Cl  --  O
         \   | \        \
          \  |  \      \
         UMg --  ClMg - OMg
         |    |
        D1    |
         | \    |
        D2  \    |
           \    D1Mg
            \    |
            D2Mg
-----------------------------------------------------------------------------

  Based on voltage-clamp recordings of NMDA receptor-mediated currents in 
  nucleated patches of  rat neocortical layer 5 pyramidal neurons (Kampa 2004), 
  this model was fit with AxoGraph directly to experimental recordings in 
  order to obtain the optimal values for the parameters.

-----------------------------------------------------------------------------

  This mod file does not include mechanisms for the release and time course
  of transmitter; it should to be used in conjunction with a sepearate mechanism
  to describe the release of transmitter and tiemcourse of the concentration
  of transmitter in the synaptic cleft (to be connected to pointer C here).

-----------------------------------------------------------------------------

  See details of NEURON kinetic models in:

  Destexhe, A., Mainen, Z.F. and Sejnowski, T.J.  Kinetic models of 
  synaptic transmission.  In: Methods in Neuronal Modeling (2nd edition; 
  edited by Koch, C. and Segev, I.), MIT press, Cambridge, 1996.


  Written by Bjoern Kampa in 2004 

-----------------------------------------------------------------------------
ENDCOMMENT

:INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	THREADSAFE
	POINT_PROCESS NMDA_Mg
	RANGE U, Cl, D1, D2, O, UMg, ClMg, D1Mg, D2Mg, OMg
	RANGE g, gmax, rb, rmb, rmu, rbMg,rmc1b,rmc1u,rmc2b,rmc2u
	GLOBAL Erev, mg, Rb, Ru, Rd1, Rr1, Rd2, Rr2, Ro, Rc, Rmb, Rmu
	GLOBAL RbMg, RuMg, Rd1Mg, Rr1Mg, Rd2Mg, Rr2Mg, RoMg, RcMg
	GLOBAL Rmd1b,Rmd1u,Rmd2b,Rmd2u,rmd1b,rmd1u,rmd2b,rmd2u
	GLOBAL Rmc1b,Rmc1u,Rmc2b,Rmc2u
	GLOBAL vmin, vmax, valence, memb_fraction
	GLOBAL Cdur
	NONSPECIFIC_CURRENT iNMDA
	RANGE Cmax, fac, ica
	USEION ca WRITE ica
}

UNITS {
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (pS) = (picosiemens)
    (nS) = (nanosiemens)
    (umho) = (microsiemens)
    (mM) = (milli/liter)
    (uM) = (micro/liter)
}

PARAMETER {

    Erev    = 0        (mV)    : reversal potential
    gmax    = 0.0005      (umho)    : maximal conductance
    mg    = 1      (mM)    : external magnesium concentration
    vmin     = -120    (mV)
    vmax     = 100    (mV)
    valence = -1.4		:-2         parameters of voltage-dependent Mg block  THIS NASSI 27/04/10
    memb_fraction = 0.8
    Cdur    = 1    (ms)    :transmitter duration 
	Cmax    = 1    (mM)    
: Rates

    Rb        = 10e-3       (/uM /ms)    : binding         
    Ru        = 8.6e-3	:8.6e-3    ::5.6e-3       (/ms)    : unbinding        
    Ro        = 10e-3       (/ms)    : opening
    Rc        = 273e-3       (/ms)    : closing 
    Rd1        = 2.2e-3       (/ms)    : fast desensitisation
    Rr1        = 1.6e-3       (/ms)    : fast resensitisation
    Rd2         = 0.43e-3     (/ms)    : slow desensitisation  
    Rr2         = 0.5e-3    (/ms)    : slow resensitisation
    Rmb        = 0.05e-3    (/uM /ms)    : Mg binding Open
    Rmu        = 12800e-3    (/ms)    : Mg unbinding Open
    Rmc1b        = 0.00005e-3    (/uM /ms)    : Mg binding Closed
    Rmc1u        = 2.438312e-3    (/ms)    : Mg unbinding Closed
    Rmc2b        = 0.00005e-3    (/uM /ms)    : Mg binding Closed2
    Rmc2u        = 5.041915e-3    (/ms)    : Mg unbinding Closed2
    Rmd1b        = 0.00005e-3    (/uM /ms)    : Mg binding Desens1 
    Rmd1u        = 2.98874e-3    (/ms)    : Mg unbinding Desens1    
    Rmd2b        = 0.00005e-3    (/uM /ms)    : Mg binding Desens2
    Rmd2u        = 2.953408e-3    (/ms)    : Mg unbinding Desens2
    RbMg        = 10e-3        (/uM /ms)    : binding with Mg
    RuMg        = 17.1e-3    (/ms)    : unbinding with Mg
    RoMg        = 10e-3        (/ms)    : opening with Mg  
    RcMg        = 548e-3    (/ms)    : closing with Mg
    Rd1Mg        = 2.1e-3    (/ms)    : fast desensitisation with Mg   NASSI
    Rr1Mg        = 0.87e-3    (/ms)    : fast resensitisation with Mg
    Rd2Mg        = 0.26e-3    (/ms)    : slow desensitisation with Mg  NASSI
    Rr2Mg        = 0.42e-3    (/ms)    : slow resensitisation with Mg
}

ASSIGNED {
    v        (mV)    : postsynaptic voltage
    iNMDA         (nA)    : current = g*(v - Erev)
    g         (umho)    : conductance
    C         (mM)    : pointer to glutamate concentration
    rb        (/ms)   : binding, [glu] dependent
    rmb        (/ms)    : blocking V and [Mg] dependent
    rmu        (/ms)    : unblocking V and [Mg] dependent
    rbMg        (/ms)    : binding, [glu] dependent
    rmc1b        (/ms)    : blocking V and [Mg] dependent
    rmc1u        (/ms)    : unblocking V and [Mg] dependent
    rmc2b        (/ms)    : blocking V and [Mg] dependent
    rmc2u        (/ms)    : unblocking V and [Mg] dependent
    rmd1b        (/ms)    : blocking V and [Mg] dependent
    rmd1u        (/ms)    : unblocking V and [Mg] dependent
    rmd2b        (/ms)    : blocking V and [Mg] dependent
    rmd2u        (/ms)    : unblocking V and [Mg] dependent
	
	ica (nA): calcium current through NMDA receptors
	synon
	Rinf
	Rtau
	fac : facilitation
}

STATE {
    : Channel states (all fractions)
    U        : unbound
    Cl        : closed
    D1        : desensitised 1
    D2        : desensitised 2
    O        : open
    UMg        : unbound with Mg
    ClMg        : closed with Mg
    D1Mg        : desensitised 1 with Mg
    D2Mg        : desensitised 2 with Mg
    OMg        : open with Mg
	
}

INITIAL {
	U = 1
	C =0

}

BREAKPOINT {
    SOLVE kstates METHOD sparse

    g = gmax * O 
:    * (1 + fac)
    iNMDA = g * (v - Erev)
		
	:ica= (0.015)*g*(v- 120 (mV))  :cHECK THIS FOR ca++
}

KINETIC kstates {

    rb     = Rb     * (1e3) * C
    rbMg     = RbMg     * (1e3) * C
    rmb     = Rmb     * mg * (1e3) * exp((v-40) * valence * memb_fraction /25)
    rmu     = Rmu     * exp((-1)*(v-40) * valence * (1-memb_fraction) /25)
    rmc1b     = Rmc1b * mg * (1e3) * exp((v-40) * valence * memb_fraction /25)
    rmc1u     = Rmc1u * exp((-1)*(v-40) * valence * (1-memb_fraction) /25)
    rmc2b     = Rmc2b * mg * (1e3) * exp((v-40) * valence * memb_fraction /25)
    rmc2u     = Rmc2u * exp((-1)*(v-40) * valence * (1-memb_fraction) /25)
    rmd1b     = Rmd1b * mg * (1e3) * exp((v-40) * valence * memb_fraction /25)
    rmd1u     = Rmd1u * exp((-1)*(v-40) * valence * (1-memb_fraction) /25)
    rmd2b     = Rmd2b * mg * (1e3) * exp((v-40) * valence * memb_fraction /25)
    rmd2u     = Rmd2u * exp((-1)*(v-40) * valence * (1-memb_fraction) /25)

    ~ U <-> Cl    (rb,Ru)
    ~ Cl <-> O    (Ro,Rc)
    ~ Cl <-> D1    (Rd1,Rr1)
    ~ D1 <-> D2    (Rd2,Rr2)
    ~ O <-> OMg    (rmb,rmu)
    ~ UMg <-> ClMg     (rbMg,RuMg)
    ~ ClMg <-> OMg     (RoMg,RcMg)
    ~ ClMg <-> D1Mg (Rd1Mg,Rr1Mg)
    ~ D1Mg <-> D2Mg (Rd2Mg,Rr2Mg)
    ~ U <-> UMg     (rmc1b,rmc1u)
    ~ Cl <-> ClMg    (rmc2b,rmc2u)
    ~ D1 <-> D1Mg    (rmd1b,rmd1u)
    ~ D2 <-> D2Mg    (rmd2b,rmd2u)

    CONSERVE U+Cl+D1+D2+O+UMg+ClMg+D1Mg+D2Mg+OMg = 1
}


NET_RECEIVE(weight, on, nspike, r0, t0 (ms), spiket0) {

    : flag is an implicit argument of NET_RECEIVE and  normally 0
	
	: NMDA exp(-(((x-50)/65)))
	: AMPA 
	if (flag == 0) { : a spike, so turn on if not already in a Cdur pulse
		nspike = nspike + 1
		if (!on) {
			on = 1
			state_discontinuity(C,C + Cmax)
		}
		
		fac =  exp(-(((t-spiket0)-50)/60))
		spiket0 = t
		: come again in Cdur with flag = current value of nspike
		net_send(Cdur, nspike)
	}
	
	if (flag == nspike) { : if this associated with last spike then turn off


		on = 0
		state_discontinuity(C, C - Cmax)
	}
	gmax = weight 
}



