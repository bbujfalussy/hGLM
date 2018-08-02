## integrates a spike trains with an alpha kernel - parametrised by its 
## We need these functions...
source('./AlphaFiltAmpl.R')
source('./AlphaFilt.R')
source('./DexpFilt.R')

#############################################
int.spikes <- function(X, dt, Wc, Ww, Tau, delay.t=NA, grad.Tau=F, grad.delay=F, alpha.ampl=T, dTau=NA){
## This function intagrates the synaptic input arriving to a given subunit
## INPUT:
## X: NxL binary input matrix of presynaptic spikes; N: # of input neurons; L: # of timesteps
## dt: the time resolution of X in milliseconds
## Wc: a vector with the indices of the presynaptic neurons.
## Ww: either a scalar or a vector with the synaptic weights.
## Tau: time constant of the synapses - positive a scalar
## delay.t: propagation delay of the synapses - a positive scalar 
## grad.Tau: specifies whether input is filtered with the synaptic kernel or its derivative wrt. Tau. if T, alpha must be T
## alpha.ampl: filter is parametrised by keeping the amplitude (instead of the integral) independent of Tau

## OUTPUT:
## out: a vector of length L  - the total synaptic input to a given subunit
	if (length(Tau)>1){
		warning("Tau has a wrong length, only the first element will be used.")
		Tau <- Tau[1]
	}

	if (length(dTau)>1){
		warning("dTau has a wrong length, only the first element will be used.")
		dTau <- dTau[1]
	}
	
	if (length(delay.t)>1){
		warning("delay.t has wrong length, only the first element of delay.t will be used at every synapse")
		delay.t <- delay.t[1]
	}

	if (is.na(dTau)) alphasyn <- T else alphasyn <- FALSE
	if (is.na(delay.t)) syn.delay <- F else syn.delay <- T
	if (syn.delay==F) delay.t <- 0 # this is the default delay

	L <- ncol(X)
	Tmax <- L / dt
	N <- nrow(X)
	y <- rep(0, L)

	if (is.null(Wc)) out <- y else {
		n <- length(Wc)
		## test of Ww
		if (length(unique(Wc)) != length(Wc)) stop("There are repeated elements in Wc!")
		if (length(Ww)==1) Ww <- rep(Ww, n)
		if (n>nrow(X)) stop ("length of Wc must be smaller than N")
		if (max(Wc)>nrow(X)) stop ("max of Wc must be smaller than N")
		if (length(Ww) != n) stop ("length of Wc must be the same as the length of Ww")
		
		if (!alphasyn){
			if (!(grad.Tau %in% c(0, 1,2))) stop("grad.Tau must be [0, 1, 2] with DExp synapses")
		} else if (!(grad.Tau %in% c(0,1))) stop("grad.Tau must be [0, 1] with alpha synapses")

		if (!(grad.delay %in% c(T, F))) {
			stop("grad.Tau must be T or F")			
			if (!syn.delay) stop("gradiants wrt the synaptic delay can only be computed if delay is provided!")
		}
		
		# the weights of the cells
		w <- rep(0, N)
		for (i in 1:n) w[Wc[i]] <- Ww[i] 
		if (max(w) == Inf) w[w==Inf] <- 1e100		
		if (min(w) == -Inf) w[w== (-Inf)] <- -1e100		
				
		x <- w %*% X
		if (alphasyn) {
			if (alpha.ampl){
				out <- alpha.filt.ampl(Tau, x, dt, delay.t=delay.t, grad.tau=grad.Tau, grad.delay=grad.delay) 
			} else {
				out <- alpha.filt(Tau, x, dt, delay.t=delay.t, grad.tau=grad.Tau, grad.delay=grad.delay) 						
			}
		} else {
			out <- dexp.filt(Tau, dTau, x, delay.t=delay.t, dt, grad.tau=grad.Tau, grad.delay=grad.delay) 
		}
	}
	out	
}


# ## Testing this function
# X <- cbind(matrix(0, 10, 1000), matrix(rbinom(10*1000, 1, 0.005), 10, 1000))
# for (i in 1:10) X[i, 100*i] <- 1
# source('~/Programs/DendInf/10Switching/RCode/graphics/Graphics.R')
# raster(X, dt=1)

# dt <- 1
# Wc <- seq(1,10)
# Ww <- 1
# Tau <- 20

# ## basic test
# y <- int.spikes(X, 1, Wc, Ww, Tau)
# lines(y, t="l")

# ## not all cells are connected
# Wc <- seq(1,5)
# y <- int.spikes(X, 1, Wc, Ww, Tau)
# lines(y, t="l", col=2)

# ## different weights
# Ww <- seq(1,5)
# y <- int.spikes(X, 1, Wc, Ww, Tau)
# lines(y, t="l", col=3)

# ## different time constants
# Ww <- seq(10,2, by=-2)
# Tau <- seq(1,5)*4
# y <- int.spikes(X, 1, Wc, Ww, Tau)
# lines(y+2, t="l", col=4)

# ## different delays
# Ww <- seq(10,2, by=-2)
# delay.t <- seq(4,10)*4
# y <- int.spikes(X, 1, Wc, Ww, Tau, delay.t=delay.t)
# lines(y+2, t="l", col=4)

# y <- int.spikes(X, 1, Wc, Ww=6, Tau=20, delay.t=delay.t)
# lines(y+2, t="l", col=4)

#############################
## test for warnings
## Tau shorter than Wc
# raster(X, dt=1)
# Tau <- seq(1,3)*4
# y <- int.spikes(X, 1, Wc, Ww, Tau)
# lines(y, t="l", col=1)
# ## Tau shorter than Wc - with exponential synapses
# y <- int.spikes(X, 1, Wc, Ww, Tau, dTau=10*Tau)
# lines(y, t="l", col=2)

# ## element in Wc too big
# Wc <- c(1,6, 12)
# y <- int.spikes(X, 1, Wc, Ww, Tau, dTau=10*Tau)
# lines(y, t="l", col=2)

# ## Wc too long
# Wc <- rep(1, 12)
# y <- int.spikes(X, 1, Wc, Ww, Tau, dTau=10*Tau)
# lines(y, t="l", col=2)

# ## Wc too long
# Ww <- rep(1, 2); Wc <- c(1, 2); Tau <- 10
# y <- int.spikes(X, 1, Wc, Ww, Tau)
# lines(y, t="l", col=2)

# #############################
# ## testing gradients: 
# X <- cbind(matrix(0, 10, 1000), matrix(rbinom(10*1000, 1, 0.005), 10, 1000))
# for (i in 1:10) X[i, 100*i] <- 1
# raster(X, dt=1)

# dt <- 1
# Wc <- seq(1,10)
# Ww <- 1
# Tau <- 10

# ## basic test
# y <- int.spikes(X, 1, Wc, Ww, Tau, grad.Tau=1)
# lines(10*y+2, t="l")

# ## not all cells are connected
# Wc <- seq(1,5)
# y <- int.spikes(X, 1, Wc, Ww, Tau, grad.Tau=T)
# lines(10*y+2, t="l", col=2)

# ## different weights
# Ww <- seq(1,5)
# y <- int.spikes(X, 1, Wc, Ww, Tau, grad.Tau=T)
# lines(10*y+2, t="l", col=3)

# ## different time constants
# Ww <- seq(10,2, by=-2)
# Tau <- seq(1,5)*4
# y <- int.spikes(X, 1, Wc, Ww, Tau, grad.Tau=T)
# lines(10*y+2, t="l", col=4)

# y <- int.spikes(X, 1, Wc, Ww, Tau, grad.delay=T)
# lines(10*y+2, t="l", col=4)

# ## different delays
# delay.t <- seq(0,4)*4
# raster(X, dt=1)
# y <- int.spikes(X, 1, Wc, Ww, Tau=20, delay.t=delay.t, grad.Tau=F)
# lines(y+2, t="l", col=4)

# y <- int.spikes(X, 1, Wc, Ww=Ww, Tau=20, delay.t=delay.t, grad.delay=T)
# lines(y+2, t="l", col=4)


#############################################
int2.spikes <- function(X, dt, Wc, Ww, Tau, delay.t=NA, Ww.b=NA, Tau.b=NA, grad.Tau=F, grad.Taub=F, grad.Wwb=F, grad.delay=F){
## This function intagrates the synaptic input arriving to a given subunit
## this handles location dependent parameters but only alpha synapses
## INPUT:
## X: NxL binary input matrix of presynaptic spikes; N: # of input neurons; L: # of timesteps
## dt: the time resolution of X in milliseconds
## Wc: a vector with the indices of the presynaptic neurons.
## Ww: either a scalar or a vector with the synaptic weights.
## Tau: time constant of the synapses - positive a scalar
## delay.t: propagation delay of the synapses - a positive scalar
# Ww.b: location dependence of the weights
# Tau.b: location dependence of the time constants

## grad.Tau: specifies whether input is filtered with the synaptic kernel or its derivative wrt. Tau
## grad.Taub: gradient wrt Tau.b is computed
## grad.Wwb: gradient wrt Ww.b is computed
## grad.delay: gradient wrt delay is computed

## OUTPUT:
## out: a vector of length L is one of the following: 
# - the total synaptic input to the given subunit
# - its gradient wrt Tau
# - its gradient wrt Tau.b
# - its gradient wrt Ww.b
# - its gradient wrt. delay.t
# delay.t <- NA; Ww.b <- NA; Tau.b <- NA; grad.Tau <- F; grad.Taub <- F; grad.Wwb <- F; grad.delay <- F
	# cat("int2 spikes called - ")
	if (length(Tau)>1){
		warning("Tau has a wrong length, only the first element will be used.")
		Tau <- Tau[1]
	}
	
	if (length(delay.t)>1){
		warning("delay.t has wrong length, only the first element of delay.t will be used at every synapse")
		delay.t <- delay.t[1]
	}

	if (is.na(delay.t)) syn.delay <- F else syn.delay <- T
	if (syn.delay==F) delay.t <- 0 # this is the default delay

	L <- ncol(X)
	Tmax <- L / dt
	N <- nrow(X)
	y <- rep(0, L)

	# cat("initial checks done - ")

	####################################
	## checking gradients...
	if (!(grad.Tau %in% c(T,F))) stop("grad.Tau must be [0, 1] with alpha synapses")
	if (grad.Wwb & is.na(Ww.b)) {print("Ww.b is not provided for gradient"); return(y)}
	if (grad.Taub & is.na(Tau.b)) {print("Tau.b is not provided for gradient"); return(y)}
	if (!(grad.delay %in% c(T, F))) {
		stop("grad.Tau must be T or F")			
		if (!syn.delay) stop("gradiants wrt the synaptic delay can only be computed if delay is provided!")
	}
	if ((grad.Tau + grad.Taub + grad.Wwb + grad.delay) > 1) stop("gradient is cumputed wrt only one parameter!")

	# cat("gradients checked - ")
	
	if (is.null(Wc)) {
		# print("Wc is null, no integration is performed")
		out <- y 
	} else {

		n <- length(Wc)
		if (length(unique(Wc)) != length(Wc)) stop("There are repeated elements in Wc!")
		if (n>nrow(X)) stop ("length of Wc must be smaller than N")
		if (max(Wc)>nrow(X)) stop ("max of Wc must be smaller than N")
		
		d.syn <- Wc - min(Wc); if (max(d.syn)>0) d.syn <- d.syn / max(d.syn) else d.syn <- 0

		##################################
		# the weights of the cells - w is a vector of N!
		if (length(Ww)>1) {
			Ww <- Ww[1]
			warning("Ww has a wrong length, only the first element will be used.")
		}

		w <- rep(0, N)
		if (is.na(Ww.b)) { 
			w[Wc] <- Ww
			Ww.b <- 0
		} else { # if synapses are distance dependent...
			Ww <- Ww * exp(Ww.b*d.syn)
			for (i in 1:n) w[Wc[i]] <- Ww[i] 
		}
		if (max(w) == Inf) w[w==Inf] <- 1e100		
		if (min(w) == -Inf) w[w== (-Inf)] <- -1e100		

# 		cat("weights calculated - ")
				
		##################################
		# the time constants - Tau is either a CONSTANT or a VECTOR of n
		if (length(Tau)>1) {
			warning("Tau has a wrong length, only the first element will be used.")
			Tau <- Tau[1]
		}
		
		if (is.na(Tau.b)) 	Tau.b <- 0
		Tau <- Tau * exp(Tau.b*d.syn)

		# cat("time constants calculated - ")

		#####################################
		X <- X * w # MATRIX * Vector = each column of X is multiplied by w = Matrix of size X
		xx <- matrix(X[Wc,], ncol=L) ## Ww dependence is included!

		if ((grad.Tau + grad.Taub + grad.Wwb + grad.delay) == 0){ # we integrate the inputs
			for (i in 1:n) {
					if (is.na(sum((xx[i,]^2)>0))) print (w)
					if (alpha.ampl){
						y <- y + alpha.filt.ampl(Tau[i], xx[i,], dt, delay.t=delay.t)
					} else {
						y <- y + alpha.filt(Tau[i], xx[i,], dt, delay.t=delay.t)						
					}
				}
			# cat("no grad ")
		}
		
		if (alpha.ampl){
			if (grad.Tau){
				for (i in 1:n) y <- y + exp(Tau.b * d.syn[i]) * alpha.filt.ampl(Tau[i], xx[i,], dt, delay.t=delay.t, grad.tau=T)
				# cat("Tau grad ")
			}
			
			if (grad.Taub){
				for (i in 1:n) y <- y + Tau[i] *  d.syn[i] * alpha.filt.ampl(Tau[i], xx[i,], dt, delay.t=delay.t, grad.tau=T)
				# cat("Taub grad ")
			}
			
			if (grad.Wwb){
				for (i in 1:n) y <- y + d.syn[i] * alpha.filt.ampl(Tau[i], xx[i,], dt, delay.t=delay.t)
				# cat("Wwb grad ")
			}
			
			if (grad.delay){
				for (i in 1:n) y <- y + alpha.filt.ampl(Tau[i], xx[i,], dt, delay.t=delay.t, grad.delay=T)
				# cat("delay grad ")
			}
			# cat("all done! \n")
		} else {
			if (grad.Tau){
				for (i in 1:n) y <- y + exp(Tau.b * d.syn[i]) * alpha.filt(Tau[i], xx[i,], dt, delay.t=delay.t, grad.tau=T)
				# cat("Tau grad ")
			}
			
			if (grad.Taub){
				for (i in 1:n) y <- y + Tau[i] *  d.syn[i] * alpha.filt(Tau[i], xx[i,], dt, delay.t=delay.t, grad.tau=T)
				# cat("Taub grad ")
			}
			
			if (grad.Wwb){
				for (i in 1:n) y <- y + d.syn[i] * alpha.filt(Tau[i], xx[i,], dt, delay.t=delay.t)
				# cat("Wwb grad ")
			}
			
			if (grad.delay){
				for (i in 1:n) y <- y + alpha.filt(Tau[i], xx[i,], dt, delay.t=delay.t, grad.delay=T)
				# cat("delay grad ")
			}
			# cat("all done! \n")
		}
		
		out <- y
	}
	out	
}



# # ############################################
# ### testing for the location dependent synapses
# X <- cbind(matrix(0, 10, 1000), matrix(rbinom(10*1000, 1, 0.005), 10, 1000))
# for (i in 1:10) X[i, 100*i] <- 1
# raster(X, dt=1)

# dt <- 1
# Wc <- seq(1,10)
# Ww <- 2
# Tau <- 20
# Ww.b <- 1
# Tau.b <- -1

# ## basic test
# y1 <- int2.spikes(X, 1, Wc, Ww, Tau)
# y2 <- int.spikes(X, 1, Wc, Ww, Tau)
# y3 <- int2.spikes(X, 1, Wc, Ww, Tau, Ww.b=Ww.b, Tau.b=Tau.b)
# lines(y1, t="l")
# lines(y2, t="l", col=2, lty=2)
# lines(y3, t="l", col=3, lty=3)
# range(y3-y2)
# range(y3-y1)
# range(y1-y2)

# y1 <- int2.spikes(X, 1, Wc, Ww, Tau, grad.Tau=T)
# y2 <- int.spikes(X, 1, Wc, Ww, Tau, grad.Tau=T)
# y3 <- int2.spikes(X, 1, Wc, Ww, Tau, Ww.b=Ww.b, Tau.b=Tau.b, grad.Tau=T)
# lines(y1, t="l")
# lines(y2, t="l", col=2, lty=2)
# lines(y3, t="l", col=3, lty=3)
# range(y3-y2)
# range(y3-y1)

# y1 <- int2.spikes(X, 1, Wc, Ww, Tau, grad.delay=T)
# y2 <- int.spikes(X, 1, Wc, Ww, Tau, grad.delay=T)
# y3 <- int2.spikes(X, 1, Wc, Ww, Tau, Ww.b=Ww.b, Tau.b=Tau.b, grad.delay=T)
# lines(y1, t="l")
# lines(y2, t="l", col=2, lty=2)
# lines(y3, t="l", col=3, lty=3)
# range(y3-y2)
# range(y3-y1)

# ## different time constants
# y <- int2.spikes(X, 1, Wc, Ww, Tau, Ww.b=Ww.b, Tau.b=Tau.b, grad.Taub=T)
# lines(10*y+2, t="l", col=4)

# y <- int2.spikes(X, 1, Wc, Ww, Tau, Ww.b=Ww.b, Tau.b=Tau.b, grad.Wwb=T)
# lines(10*y+2, t="l", col=4)

