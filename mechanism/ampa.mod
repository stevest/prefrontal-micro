TITLE simple AMPA receptors

COMMENT
-----------------------------------------------------------------------------

	Simple model for glutamate AMPA receptors
	=========================================

  - FIRST-ORDER KINETICS, FIT TO WHOLE-CELL RECORDINGS

    Whole-cell recorded postsynaptic currents mediated by AMPA/Kainate
    receptors (Xiang et al., J. Neurophysiol. 71: 2552-2556, 1994) were used
    to estimate the parameters of the present model; the fit was performed
    using a simplex algorithm (see Destexhe et al., J. Computational Neurosci.
    1: 195-230, 1994).

  - SHORT PULSES OF TRANSMITTER (0.3 ms, 0.5 mM)

    The simplified model was obtained from a detailed synaptic model that 
    included the release of transmitter in adjacent terminals, its lateral 
    diffusion and uptake, and its binding on postsynaptic receptors (Destexhe
    and Sejnowski, 1995).  Short pulses of transmitter with first-order
    kinetics were found to be the best fast alternative to represent the more
    detailed models.

  - ANALYTIC EXPRESSION

    The first-order model can be solved analytically, leading to a very fast
    mechanism for simulating synapses, since no differential equation must be
    solved (see references below).



References

   Destexhe, A., Mainen, Z.F. and Sejnowski, T.J.  An efficient method for
   computing synaptic conductances based on a kinetic model of receptor binding
   Neural Computation 6: 10-14, 1994.  

   Destexhe, A., Mainen, Z.F. and Sejnowski, T.J. Synthesis of models for
   excitable membranes, synaptic transmission and neuromodulation using a 
   common kinetic formalism, Journal of Computational Neuroscience 1: 
   195-230, 1994.

-----------------------------------------------------------------------------
ENDCOMMENT



NEURON {
    THREADSAFE
	POINT_PROCESS GLU         
	RANGE R, gmax, g, ina     
	USEION na WRITE ina   
	NONSPECIFIC_CURRENT  iglu             : i
	GLOBAL Cdur, Alpha, Beta, Erev, Rinf, Rtau
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {
        Cmax	= 1	(mM)		: max transmitter concentration
:	Cdur	= 0.3	(ms)		: transmitter duration (rising phase)
	Cdur	= 1.1	(ms)		: transmitter duration (rising phase)
:	Alpha	= 0.94	(/ms)	: forward (binding) rate
	Alpha	= 10	(/ms)	: forward (binding) rate
:	Beta	= 0.18	(/ms)		: backward (unbinding) rate
	Beta	= 0.5	(/ms)		: backward (unbinding) rate
	Erev	= 0	(mV)		:0 reversal potential
}


ASSIGNED {
	v		(mV)		: postsynaptic voltage
	iglu 		(nA)		: current = g*(v - Erev)     :i
	g 		(umho)		: conductance
	Rinf				: steady state channels open
	Rtau		(ms)		: time constant of channel binding
	synon
	gmax
	ina
	ica
}

STATE {Ron Roff}

INITIAL {
        Rinf = Cmax*Alpha / (Cmax*Alpha + Beta)
       	Rtau = 1 / ((Alpha * Cmax) + Beta)
	synon = 0
}

BREAKPOINT {
:printf("AMPA Ron %f Roff %f\n",Ron, Roff)
	SOLVE release METHOD cnexp
	g = (Ron + Roff)*1(umho)
	iglu = g*(v - Erev)  :i
	ina = 0.9*iglu
	iglu = 0.1*iglu

}

DERIVATIVE release {
	Ron' = (synon*Rinf - Ron)/Rtau
	Roff' = -Beta*Roff
}

: following supports both saturation from single input and
: summation from multiple inputs
: if spike occurs during CDur then new off time is t + CDur
: ie. transmitter concatenates but does not summate
: Note: automatic initialization of all reference args to 0 except first

NET_RECEIVE(weight, on, nspike, r0, t0 (ms)) {
INITIAL{
	nspike=0
	on = 0
}
	: flag is an implicit argument of NET_RECEIVE (equals 0 if external, else is a self send message)
        if (flag == 0) { : a spike, so turn on if not already in a Cdur pulse
		nspike = nspike + 1
		if (!on) {
			synon = synon + weight
			r0 = r0*exp(-Beta*(t - t0))
			:state_discontinuity(Ron, Ron + r0) :deprecated and not thread safe!
			Ron = Ron + r0
			:state_discontinuity(Roff, Roff - r0):deprecated and not thread safe!
			Roff = Roff - r0
			t0 = t
			on = 1
			: come again in Cdur with flag = current value of nspike
			net_send(Cdur, nspike)
		}:else{
			:net_move(t+Cdur)
		:}
		

        }
	if (flag == nspike) { : if this associated with last spike then turn off
		synon = synon - weight
		r0 = weight*Rinf + (r0 - weight*Rinf)*exp(-(t - t0)/Rtau)
		state_discontinuity(Ron, Ron - r0) :deprecated and not thread safe!
		:Ron = Ron - r0
		state_discontinuity(Roff, Roff + r0) :deprecated and not thread safe!
		:Roff = Roff + r0
		t0 = t
		on = 0
	}
gmax = weight
}

