: Slowly inactivating K+ channel

NEURON {
	SUFFIX IKs
	USEION k READ ki, ko WRITE ik
	RANGE gKsbar, ik, gk
	
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
        (mM) = (milli/liter)
	
}
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
PARAMETER {
	v (mV)
	dt (ms)
	gKsbar= 0.00014 (mho/cm2) <0,1e9>
	
}


STATE {
	a b
}


ASSIGNED {
	ik (mA/cm2)
	IKS_ainf
	IKS_binf
	IKS_atau (ms)
	IKS_btau (ms)
	gk (mho/cm2)
	ek  (mV)
	ki (mM)
	ko  (mM)
}



INITIAL {
	rate(v)
	a = IKS_ainf
	b = IKS_binf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
		
	gk = gKsbar * a * b
	ek = 25 * log(ko/ki)
	ik = gk*(v-ek)
	
}

DERIVATIVE states {
	rate(v)
	a' = (IKS_ainf-a)/IKS_atau
	b' = (IKS_binf-b)/IKS_btau
}
UNITSOFF

PROCEDURE rate(v (mV)) {LOCAL IKS_va, IKS_vb, IKS_vc, IKS_vd
	IKS_va = v + 34
	IKS_vb = v + 65
	IKS_vd = v + 63.6
	

if (fabs(IKS_va)<1e-04){ IKS_va = IKS_va+0.00001 }
	   IKS_ainf = 1/(1 + exp(-IKS_va/6.5))
	   IKS_atau = 10
	

if (fabs(IKS_vb)<1e-04){ IKS_vb = IKS_vb+0.00001 }
	   IKS_binf = 1/(1 + exp(IKS_vb/6.6))

 
if (fabs(IKS_vd)<1e-04){ IKS_vd = IKS_vd+0.00001 }
	   IKS_btau = 200 + 3200 / (1 + exp(-IKS_vd/4))
}
UNITSON







