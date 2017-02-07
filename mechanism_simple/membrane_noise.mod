COMMENT
see also "random number generation and thread safety":
http://www.neuron.yale.edu/phpbb/viewtopic.php?f=31&t=2136

Since this is an electrode current, positive values of i depolarize the cell
and in the presence of the extracellular mechanism there will be a change
in vext since i is not a transmembrane current but a current injected
directly to the inside of the cell.

Edited by Stefanos to create reproducable results accross HPC nodes, as discussed here:
http://www.neuron.yale.edu/phpBB/viewtopic.php?f=31&t=2136


ENDCOMMENT

NEURON {
	POINT_PROCESS membNoise
	RANGE del, dur, pkamp, freq, phase, bias, pkampMax
	ELECTRODE_CURRENT i :ELECTRODE_CURRENT
	THREADSAFE : only true if every instance has its own distinct Random
	POINTER randObjPtr
}

UNITS {
	(nA) = (nanoamp)
}

PARAMETER {
	del=0   (ms)
	dur=0   (ms)
	pkamp=0 (nA) :1
	pkampMax=0 (nA) :1
	freq=0  (Hz)
	phase=0 (rad)
	bias=0  (nA)
	:randomEvent = 0
	PI=3.14159265358979323846
}

ASSIGNED {
        i (nA)
	randObjPtr
}

BREAKPOINT {
	if (t < del) {
		i=0   
	}else{ 
		if (t < del+dur) {
			randGen()
			i = pkamp*sin(2*PI*freq*(t-del)/1000+phase)+bias
			:printf("I = %.10f\n",i)
		}else{ 
			i = 0
		}
	}
} 

VERBATIM
	/* necessary lines for neuron not to get stub functions from nrnnoiv.c */
	double nrn_random_pick(void* r);
	void* nrn_random_arg(int argpos);
ENDVERBATIM 

PROCEDURE randGen() {
	VERBATIM
	if (_p_randObjPtr) {
		/*
		:Supports separate independent but reproducible streams for
		: each instance. However, the corresponding hoc Random
		: distribution MUST be set to Random.uniform(0,1)
		*/
		pkamp = pkampMax * nrn_random_pick(_p_randObjPtr);
	}else{
		hoc_execerror("Random object ref not set correctly for randObjPtr"," only via hoc Random");
	}
	ENDVERBATIM
}

PROCEDURE getRandObjPtr() {
	VERBATIM
	void** pv4 = (void**)(&_p_randObjPtr);
	if (ifarg(1)) {
		*pv4 = nrn_random_arg(1);
	}else{
		*pv4 = (void*)0;
	}
	ENDVERBATIM
}
