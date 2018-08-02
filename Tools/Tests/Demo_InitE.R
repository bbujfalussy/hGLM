
# This Utils library contains generic functions used by most of the functions in the Tools folder
# There are two functions in the Utils_hGLM.R file that are useful for initialising the training of hGLM models.
# Here I will illustrate how this initiaalisation could be done efficiently.
#
# The two functions: 
# decouple.pars(): converts a hLN with 1 synapse/subunit to 
#		 a similar hLN with potentially many synapses per subunit.
# 		The old and the new hLN should have the SAME number of SUBUNITS
	
	# arguments: Jc: the hierarchy of the new model - should be the same for the new and the old hLN!
		# Wc.e: mapping between excitatory synapses and subunits in the new hLN
		# Wc.i: mapping between inhibitory synapses and subunits in the new hLN
		# synpars: list, containing the coupled parameters, 1 per subunit!


# reshape.pars.syn2(): a function that reshapes the 1Lsyn parameters into a hierarchical architecture
	# the old and the new hLN should have the SAME number of SYNAPSES
	
	# arguments: Jc: the hierarchy of the new model - different for the new and the old hLN!
		# Wc.e: mapping between excitatory synapses and subunits in the new hLN
		# Wc.i: mapping between inhibitory synapses and subunits in the new hLN
		# synpars: list, containing the parameters in the 1-layer model.


# The combination of the two function can be used to extend the hLN architecture. 
# We illustrate it here.

source('../Train_hGLM.R', chdir=T)
source('../../Utils/Spike_Basis.R', chdir = TRUE)
library(viridis)

############################################################
# setting up the model architectures
# We will have a 3-layer hierarchy with 10 E and 10 I inputs

L <- 10000
X <- matrix(rbinom(10*L, 1, 0.005), 10, L)
raster(X, dt=1)

# we define the reference hLN for which we want to learn the parameters
Jc3 <- c(0,1,1); 
Wc.e3 <- list(1:4, 5, 6:10); Wc.i3 <- list(NULL, NULL, NULL)
pars3 <- list(Jc=Jc3, Wc.e=Wc.e3, Wc.i=Wc.i3)
refpars <- init.hGLM(X, 1, pars3, rand.seed=1, lin.soma=T, rescale.soma=T, nSD=2, mean.v=-20, sd.v=1, alphasyn=T)
lines(refpars$v+25, col=2)


#################################
# 1. train a simple, 1N model
Jc1 <- c(0)
Wc.e1 <- list(1); Wc.i1 <- list(NULL)
initPars.1N <- list(Jc=Jc1, Wc.e=Wc.e1, Wc.i=Wc.i1)

X1 <- matrix(colSums(X[1:10,]), 1)
train2 <- list(X=X1, v=refpars$v)
pars.1N <- train.hGLM(initPars.1N, train=train2, test=train2, dt=1)

# training is quick, but the resuls are not very good....
v.1N <- sim.hGLM(X1, dt=1, pars.1N$pars.sub)
lines(v.1N + 25, col=3)

#################################
# 2.  decouple the parameters of the synapses:
Wc.e10 <- list(seq(1,10))
Wc.i10 <- list(NULL)
initpars.decoupled <- decouple.pars(Jc1, Wc.e10, Wc.i10, pars.1N$pars.sub)

# decoupling has no effect, we simulate the decoupled model:
v.1Ndecoupled <- sim.hGLM(X, dt=1, initpars.decoupled)
lines(v.1Ndecoupled + 25, col=4)

#################################
# 3.  reshape the parameters to the new architecture:
initp.reshaped <- reshape.pars.syn2(Jc=Jc3, Wc.e=Wc.e3, Wc.i=Wc.i3, initpars.decoupled) 

# this is now reshaped, but the hGLM output is very different from the original:
v.3Nreshaped <- sim.hGLM(X, dt=1, initp.reshaped)
lines(v.3Nreshaped + 25, col=4)

#################################
# 4. make the hierarchy linear to match the 1N model's output
initp.reshaped.rescaled <- init.hGLM(X, 1, initp.reshaped, lin.soma=F, rescale.soma=F, nSD = 1)
v.3Nreshaped.rescaled <- sim.hGLM(X, dt=1, initp.reshaped.rescaled)
lines(v.3Nreshaped.rescaled + 25, col=6)

# The parameter nSD controls the linearity of the subunits - larger value makes the subunits more linear.
plot(refpars$v, col=1, t='l')
lines(v.1N, col=2)
lines(v.3Nreshaped.rescaled, col=3)

#################################
# 5. train the new hLN with complex hierarchy
train <- list(X=X, v=refpars$v)
pars.3N <- train.hGLM(initp.reshaped.rescaled, train=train, test=train, dt=1, linear.som=T)

####################################
## random initialisation can also work though
pars.3N <- train.hGLM(pars3, train=train, test=train, dt=1, seed=8, linear.soma=T, maxit=1000)

