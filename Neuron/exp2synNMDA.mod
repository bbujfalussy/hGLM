COMMENT
Two state kinetic scheme synapse described by rise time tau1,
and decay time constant tau2. The normalized peak condunductance is 1.
Decay time MUST be greater than rise time.

The solution of A->G->bath with rate constants 1/tau1 and 1/tau2 is
 A = a*exp(-t/tau1) and
 G = a*tau2/(tau2-tau1)*(-exp(-t/tau1) + exp(-t/tau2))
	where tau1 < tau2

If tau2-tau1 -> 0 then we have a alphasynapse.
and if tau1 -> 0 then we have just single exponential decay.

The factor is evaluated in the
initial block such that an event of weight 1 generates a
peak conductance of 1.

Because the solution is a sum of exponentials, the
coupled equations can be solved as a pair of independent equations
by the more efficient cnexp method.

***
Modified to separate ica from i
T Branco 2012
***

ENDCOMMENT

NEURON {
	POINT_PROCESS Exp2SynNMDA
:	USEION ca READ eca WRITE ica
	RANGE tau1, tau2, e, i
	NONSPECIFIC_CURRENT i
	RANGE g, caf
}

: caf is the fraction to total current carried by calcium

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau1=.1 (ms) <1e-9,1e9>
	tau2 = 10 (ms) <1e-9,1e9>
	e=0	(mV)
	mg=1    (mM)		: external magnesium concentration
    pf = 0.03  (1)      : adjusted to give 15% ica at -60 mV
}

ASSIGNED {
	v (mV)
	i (nA)
	g (uS)
	eca (mV)
	ica (nA)
	factor
}

STATE {
	A (uS)
	B (uS)
}

INITIAL {
	LOCAL tp
	if (tau1/tau2 > .9999) {
		tau1 = .9999*tau2
	}
	A = 0
	B = 0
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = B - A
	i = g*mgblock(v)*(v - e)*(1-pf)
:    ica = g*mgblock(v)*(v - eca)*pf
}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
}

FUNCTION mgblock(v(mV)) {
	TABLE 
	DEPEND mg
	FROM -140 TO 80 WITH 1000

	: from Jahr & Stevens

	mgblock = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))
}

NET_RECEIVE(weight (uS)) {
	A = A + weight*factor
	B = B + weight*factor
}
