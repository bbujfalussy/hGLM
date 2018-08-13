This library is designed to simulate and train hierarchical Generalised Linear Models (hGLMs).

It includes functions to simulate, initialise and train these models as well as other functions, necessary for simulating a hGLM model or to evaluate its derivative wrt. its parameters. The spiking is optional, so the same library handles hLNs as well.

Folders:
	Tools: main functions for hGLM simulation, evaluation etc.
	Utiles: generic functions used by most of the functions in the Tools folder
	Graphics: functions to plot the model structure for hGLMs, inut and output data (not yet implemented)
	Data: demo data from biophysical models - see the Ujfalussy et al., 2017 Fig. 2. for description of the biophysical model and the stimulation protocol

	Neuron: python and neuron files used to simulate the biophysical model
	Inputs: R scripts used to generate inputs to the biophysical model

##################################
## usage: 
examples are given in the valid_hGLN.R
fitting data from biophysical models is in Fit_hGLN.R

1. defining hGLM model:
A hGLM model is defined by its subunit structure (vector, Jc) and the inputs of the subunits (Wc.e and Wc.i).
Jc: a vector of length M, representing the subunit's parent. The root has 0. Jc <- c(0, 1, 1, 2, 2) represents the following graph:
1 - 2 - 4
  \  \
   3  5

Wc.e, Wc.i: lists of M vectors. Each vecor encodes the presynatic excitatory (Wc.e) or inhibitory (Wc.i) neurons connected to the given subunit. Repeats are allowed, but neurons must be either excitatory or inhibitory.

The input spike train (X) is an M x T matrix with elements representing the number of spikes fired by a given pre cell in a time bin. 

2. Initialising the parameters (time constant, synaptic weights and nonlinearities)
Tools/init_hGLM.R: initialize the hGLM model. Either randomly, or to match the response of a 1-layer model

3. Simulating the model.
Tools/sim_hGLM_sub.R: simulate the subthreshold response of a hGLM model to synaptic inputs
Tools/sim_hGLM_sp.R: simulate the spiking response of a hGLM model to synaptic inputs

4. Training the model: 
Tools/Train_hGLM.R: wrapper for fitting hGLM models to data. 

More examples can be found in the Tools/Tests folder, where different scripts test the functionality of the library.