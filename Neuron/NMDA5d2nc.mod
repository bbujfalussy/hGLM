
TITLE detailed model of glutamate NMDA receptors

COMMENT
-----------------------------------------------------------------------------

Kinetic model of NMDA receptors
===============================
Based on the simplification of the 
	
	10-state gating model:
	Kampa et al. (2004) J Physiol
	
	by eliminating the dynamics of the Mg-block and hence simplifying the 
	dynamics to 5-state. 
	
	U -- Cl -- O
	     |
      	 D1
	     |
	     D2
		   
Voltage dependence of Mg2+ block:
	Jahr & Stevens 1990. J Neurosci 10: 1830.
	Jahr & Stevens 1990. J Neurosci 10: 3178.
	
-----------------------------------------------------------------------------
	
  The original model is based on voltage-clamp recordings of NMDA receptor-mediated  
  currents in nucleated patches of  rat neocortical layer 5 pyramidal neurons (Kampa 2004), 
  this model was fit with AxoGraph directly to experimental recordings in 
  order to obtain the optimal values for the parameters.
  
  -----------------------------------------------------------------------------
  
  Rates modified for near physiological temperatures with Q10 values from
  O.Cais et al 2008, Mg unbinding from Vargas-Caballero 2003, opening and
  closing from Lester and Jahr 1992.

  Tiago Branco 2010
  
  -----------------------------------------------------------------------------
  
  Modified by Balazs B Ujfalussy in 2016 to 
  
  - include a simple mechnism of transmitter release based on the exampl 10.6. 
  of the Neuron book. The synapse can be driven by NetCon object. (Note, that there 
  is some redundancy about the synaptic weight, since it is controlled by the 
  parameter gmax, Cmax and the weight of the NetCon object)
  
  - Make the dynamics of the voltage-dependent Mg binding and unbinding instantaneous 
  (based on the arguments in   
  Destexhe, A., Mainen, Z.F. and Sejnowski, T.J.  Kinetic models of 
  synaptic transmission.  In: Methods in Neuronal Modeling (2nd edition; 
  edited by Koch, C. and Segev, I.), MIT press, Cambridge, 1998, pp 1-25.)
  
  -----------------------------------------------------------------------------
  
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS NMDA5d2nc
	RANGE Cglut, Cmax, Cmin, Cdur
	RANGE U, Cl, D1, D2, O, B
	RANGE g, gmax, rb
:	GLOBAL Erev, mg, Rb, Ru, Rd, Rr, Ro, Rc
:	GLOBAL vmin, vmax
	RANGE Erev, mg, Rb, Ru, Rd1, Rd2, Rr1, Rr2, Ro, Rc
	RANGE vmin, vmax
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(umho) = (micromho)
	(mM) = (milli/liter)
	(uM) = (micro/liter)
}

PARAMETER {

	Erev	= 0    (mV)	: reversal potential
	gmax	= 500  (pS)	: maximal conductance
	mg	= 1    (mM)	: external magnesium concentration
	vmin    = -120	(mV)
	vmax    = 100	(mV)
	Cmax	= 1   (mM)     : maximal transmitter concentration at the synapse
	Cmin	= 0   (mM)     : residual transmitter concentration at the synapse
	Cdur	= 1   (ms)     : duration of transmitter release 
	
: Rates

	: : Destexhe, Mainen & Sejnowski, 1996
	: Rb	= 5e-3    (/uM /ms)	: 0.005  binding 		
	: Ru	= 12.9e-3  (/ms)	: 0.0129 unbinding		
	: Rd	= 8.4e-3   (/ms)	: 0.0084 desensitization
	: Rr	= 6.8e-3   (/ms)	: 0.0068 resensitization 
	: Ro	= 46.5e-3   (/ms)	: 0.046  opening
	: Rc	= 73.8e-3   (/ms)	: 0.0738 closing
    
    : : Kampa et al., 2004 - rates without Mg
    : Rb	= 10e-3    (/uM /ms)	: binding 		
	: Ru	= 20.16e-3  (/ms)	: unbinding
	: Ro	= 46.5e-3   (/ms)	: opening
	: Rc	= 91.6e-3   (/ms)	: closing
	
	: Rd1	= 22.66e-3   (/ms)	: fast desensitization
	: Rr1	= 7.36e-3   (/ms)	: fast resensitization 
	: Rd2	= 4.43e-3   (/ms)	: slow desensitization
	: Rr2	= 2.3e-3   (/ms)	: slow resensitization 	
	
    : Kampa et al., 2004 - rates with Mg
    Rb	= 10e-3    (/uM /ms)	: binding 		
	Ru	= 61.56e-3  (/ms)	: unbinding
	Ro	= 46.5e-3   (/ms)	: opening
	Rc	= 91.6e-3   (/ms)	: closing
	
	Rd1	= 21.66e-3   (/ms)	: fast desensitization
	Rr1	= 4.002e-3   (/ms)	: fast resensitization 
	Rd2	= 2.67e-3   (/ms)	: slow desensitization
	Rr2	= 1.93e-3   (/ms)	: slow resensitization 	
}


ASSIGNED {
	v		(mV)		: postsynaptic voltage
	i 		(nA)		: current = g*(v - Erev)
	g 		(pS)		: conductance
	Cglut 		(mM)		: pointer to glutamate concentration

	rb		(/ms)    : binding
}

STATE {
	: Channel states (all fractions)
	U		: unbound
	Cl		: closed - single, corresponds to the old, double bound
	D1		: desensitized 1
	D2		: desensitized 2
	O		: open

	B		: fraction free of Mg2+ block
}

INITIAL {
	rates(v)
	U = 1
}

BREAKPOINT {
	rates(v)
	SOLVE kstates METHOD sparse

	g = gmax * O * B
	i = (1e-6) * g * (v - Erev)
}

KINETIC kstates {
	
	rb = Rb * (1e3) * Cglut 

	~ U <-> Cl	(rb,Ru)
	~ Cl <-> D1	(Rd1,Rr1)
	~ D1 <-> D2	(Rd2,Rr2)
	~ Cl <-> O	(Ro,Rc)

	CONSERVE U+Cl+D1+D2+O = 1
}

PROCEDURE rates(v(mV)) {
	TABLE B
	DEPEND mg
	FROM vmin TO vmax WITH 200

	: from Jahr & Stevens

	B = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))
}
    
    
    
NET_RECEIVE(weight, on) {
      : on == 1 if transmitter is present ("onset" state), otherwise 0
      : flag is an implicit argument of NET_RECEIVE, normally 0
      if (flag == 0) {
        : a spike happened, so increase the transmitter concentration to Cmax
        if (!on) {
          Cglut = Cmax * weight
	  on = 1
          net_send(Cdur, 1)
        } else {
          : already in onset state, so move offset time
          net_move(t+Cdur)
        }
      }
      if (flag == 1) {
        : "turn off transmitter"
        Cglut = Cmin * weight
        on = 0
      }
}