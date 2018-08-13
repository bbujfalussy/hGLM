:based on Uebachs et al 2010
:deleted eNa which overwrote ena
:deleted celsius which is not used
:parameterised Vh and k

: modified Konstantin Stadler 2010

TITLE nap

NEURON {
	SUFFIX nap
	USEION na READ ena WRITE ina
	RANGE  gbar, thegna, sh, Vh, k, mtau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

PARAMETER {
        gbar = 0.0052085 	(mho/cm2)
: 0.0478 mS/cm2 = 0.0478 x 10^(-3) mho/cm2
	sh = 0  		(mV)	:shift
	mtau = 1 		(ms) <1e-12, 1e9> 
	Vh = -35		(mV)    : TB was -52.3
	k = 6.8			(mV)     
	v 			(mV)
}


ASSIGNED {
	ina 		(mA/cm2)
	ena		(mV)
	thegna		(mho/cm2)
	minf 		(1)
}
 

STATE { m }


BREAKPOINT {
    	SOLVE states METHOD cnexp
	thegna = gbar*m
	ina = thegna * (v - ena)
} 

DERIVATIVE states {   
	calcMinf(v) 	
	m' = (minf-m)/mtau
}


INITIAL {
	calcMinf(v) 
	m=minf  
}

PROCEDURE calcMinf(v(mV)) {
	minf = 1 / ( 1+exp(-(v-Vh-sh)/k) )
}

