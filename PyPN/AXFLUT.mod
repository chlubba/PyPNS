TITLE FLUT channels

: 5/16
: Carl H. Lubba
:
: Fast K+ Channel
: Iterative equations H-H notation rest = -80 mV
:
: This FLUT channel is described in detail in:
:
: McIntyre CC, Grill WM, Sherman DL, Thakor NV
: Cellular effects of deep brain stimulation: model-based analysis of activation and inhibition
: J Neurophysiol. 91(4):1457-69, 2004


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX axflut
	NONSPECIFIC_CURRENT ikf
	RANGE gkbar, ek
	RANGE n_inf
	RANGE tau_n
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {

	gkbar   = 0.08 	(mho/cm2)
	ek      = -90.0 (mV)
	celsius		(degC)
	dt              (ms)
	v               (mV)
	anA = 0.0462
	anB = 83.2
	anC = 1.1
	bnA = 0.0824
	bnB = 66
	bnC = 10.5

}

STATE {
	n
}

ASSIGNED {
	ikf      (mA/cm2)
	n_inf
	tau_n
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ikf   = gkbar * n*n*n*n * (v - ek)
}

DERIVATIVE states {   : exact Hodgkin-Huxley equations
       evaluate_fct(v)
	n'= (n_inf - n) / tau_n
}

UNITSOFF

INITIAL {

	evaluate_fct(v)
	n = n_inf

}

PROCEDURE evaluate_fct(v(mV)) { LOCAL a,b

	a = vtrap1(v)
	b = vtrap2(v)
	tau_n = 1 / (a + b)
	n_inf = a / (a + b)

}


FUNCTION vtrap1(x) {
	if (fabs((x+anB)/anC) < 1e-6) {
		vtrap1 = anA*anC
	}else{
		vtrap1 = (anA*(x+anB)) / (1 - Exp(-(x+anB)/anC))
	}
}

FUNCTION vtrap2(x) {
	if (fabs((x+bnB)/bnC) < 1e-6) {
		vtrap2 = bnA*bnC : Ted Carnevale minus sign bug fix
	}else{
		vtrap2 = (bnA*(-(x+bnB))) / (1 - Exp((x+bnB)/bnC))
	}
}

FUNCTION Exp(x) {
	if (x < -100) {
		Exp = 0
	}else{
		Exp = exp(x)
	}
}

UNITSON
