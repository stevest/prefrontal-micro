:  Vector stream of events

NEURON {
	THREADSAFE
	ARTIFICIAL_CELL NonHomPoissonStim
	POINTER randObjPtrUniform, randObjPtrPoisson
	RANGE lambdaMax
}

ASSIGNED {
	v		(mV)
	index
	lambda 
	space
	delay
	randObjPtrUniform
	randObjPtrPoisson
	randNo
	nevents
}

PARAMETER { 
	lambdaMax
} 

INITIAL {
	nevents = 0
	index = 0
	:printf("Initial voltage is: %f\n", v)
	element()
	if (index > 0) {
		:Sent always in the next t:
		net_send(1, 1)
	}
}

BREAKPOINT {
	:printf("BLAH DE BLAH\n")
	:printf("Initial voltage is: %f\n", v)
}



NET_RECEIVE (w) {
	:printf("recieved net message\n")
	if (flag == 1) {
		:printf("net message has flag == 1\n")
		element()
		:printf("lambda taken from rate vector is: %f\n", lambda)
		if (index > 0) {
			:printf("index is valid\n")
			:Generate max lambda poisson events:
			generate_poisson_events()
			:printf("number of events generated is: %f\n",nevents)
			:Thin Poisson:
			doThin()
			:printf("number of events after thinning are: %f\n",nevents)
			if ( nevents > 0 ){
				:printf("Activating mechanism at time: %f. Index is: %f\n",t,index)
				net_event(t)
			}
			:printf("Sending event for the next time step\n")
			net_send(1, 1)
		}
	}
}

VERBATIM
extern double* vector_vec();
extern int vector_capacity();
extern void* vector_arg();
/* /!\ necessary lines for neuron not to get stub functions from nrnnoiv.c */
double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);
ENDVERBATIM

PROCEDURE element() {
VERBATIM	
  { void* vv; int i, size; double* px;
	i = (int)index;
	// For positive indices of input vector
	if (i >= 0) {
		//Get vector pointer
		vv = *((void**)(&space));
		if (vv) {
			size = vector_capacity(vv);
			px = vector_vec(vv);
			if (i < size) {
				// Lambda in current t is taken from vector
				lambda = px[i];
				index += 1.;
			}else{
				index = -1.;
			}
		}else{
			index = -1.;
		}
	}
  }
ENDVERBATIM
}

PROCEDURE play() {
VERBATIM
	void** vv;
	vv = (void**)(&space);
	*vv = (void*)0;
	if (ifarg(1)) {
		//printf(".play() initialized with vector\n");
		*vv = vector_arg(1);
	}
	{
		int size;
		size = vector_capacity(*((void**)(&space)));
		//printf("Playing from vector of size: %d\n", size);
	}
	
ENDVERBATIM
}
        

PROCEDURE generate_poisson_events() {
	VERBATIM
	if (_p_randObjPtrPoisson) {
		/*
		:Supports separate independent but reproducible streams for
		: each instance. However, the corresponding hoc Random
		: distribution MUST be set to Random.uniform(0,1)
		*/
		// Generates from a (hoc defined) Poisson, using the maximum lambda:
		nevents = nrn_random_pick(_p_randObjPtrPoisson);
		//printf("Generated poisson events, success!\n");
		//printf("No of events is: %f\n",nevents);
		// nevents are the pre-thined events.
	}else{
		hoc_execerror("Random object ref not set correctly for randObjPtrPoisson"," only via hoc Random");
	}
	ENDVERBATIM
}

PROCEDURE doThin() {
	VERBATIM
	if (_p_randObjPtrUniform) {
		/*
		:Supports separate independent but reproducible streams for
		: each instance. However, the corresponding hoc Random
		: distribution MUST be set to Random.uniform(0,1)
		*/
		randNo = nrn_random_pick(_p_randObjPtrUniform);
	}else{
		hoc_execerror("Random object ref not set correctly for randObjPtrPoisson"," only via hoc Random");
	}
	if (randNo < (lambda / lambdaMax)) {
		// We have a Poisson event!
	}else{
		//No event
		nevents = 0;
	}
	ENDVERBATIM
}

PROCEDURE getRandObjPtrUniform() {
	VERBATIM
	void** pUniform = (void**)(&_p_randObjPtrUniform);
	if (ifarg(1)) {
		*pUniform = nrn_random_arg(1);
		//printf("Getting Uniform Random success!\n");
	}else{
		*pUniform = (void*)0;
	}
	ENDVERBATIM
}

PROCEDURE getRandObjPtrPoisson() {
	VERBATIM
	void** pPoisson = (void**)(&_p_randObjPtrPoisson);
	if (ifarg(1)) {
		*pPoisson = nrn_random_arg(1);
		//printf("Getting Poisson Random success!\n");
	}else{
		*pPoisson = (void*)0;
	}
	ENDVERBATIM
}
