## A script to simulate response of a hGLM model
## AHP is conditional on the spike times - provided by the matrix X.ahp, which contains the basis functions convolved with the output spikes
## We need it to be efficient to be able to evaluate it inside the optimization
## To test this function, run 
#    source('./UnitTests_Init.R', chdir = TRUE)


source('../Utils/IntSpikes.R', chdir=T)
source('../Utils/Sigm.R')
source('./Init_hGLM.R')
source('../Graphics/Graphics.R')

#vv <- NULL; logpars <- list(Jw=T, Tau=T, W=F); regpars <- NULL; double <- F; scale.tau2 <- 2.8; verbose <- F; calc.minis <- F; X.ahp <- NULL

sim.hGLM <- function(X, dt, pars, vv=NULL, logpars=list(Jw=T, Tau=T, W=F), regpars=NULL, double=F, alphasyn=T, scale.tau2=2.8, verbose=F, calc.minis=F, X.ahp=NULL){
## function to simulate subthreshold response of a hGLM model
##
## args:
## X: NxT binary input matrix of presynaptic spikes; N: # of input neurons; T: # of timesteps
## dt: the time resolution of X in miliseconds
## pars: list with the following elements: 
	# v0: baseline of somatic voltage
	# Jc: M length connectivity vector for dendritic branches; M: # of dendritic branches
		# for each branch its parent is given. 1 is the root which has a value 0. 
		# e.g., a 3 layer deep binary tree has 7 subunits is given by: [0, 1, 1, 2, 2, 3, 3]
	# Jw: M length vector for the coupling weights associated with the dendritic branches
	# Wc.e, Wc.i: a list of M components. The component m is a vector indicating the neurons connected to branch m
	# Ww.e, (Ww.e2), Ww.i: LISTS of M components indicating the synaptic weights associated to the neurons connected to branch m
	# Th: M length vector of the thresholds for the sigmoidal nonlinearities in the branches
	# Tau.e, Tau.i: LISTS of M components of the (log) synaptic time constants in miliseconds
	# dTau.e: optional LISTS of M components of the (log) synaptic time constants for the dExp synapses in miliseconds
	# delay.t: optional vector of the synaptic delays in ms (default: 1/10 ms)
## vv: response, if provided only the error is calculated
## logpars: parameters Tau, Jw and Ww are strictly positive, so they can be defined by their log
## regpars: NULL (no regularisation) or list with the prior mean and inverse variance for the log of Ww and Tau
	# Ww.e, (Ww.e2), Ww.i, Tau.e, Tau.i, alpha.Ww, alpha.Tau
## double: logical; simple or double alpha kernels 
## scale.tau2: if double kernels are used, the scaling factor between their time constants
## verbose: regularisation details are provided
## calc.minis: also returns the calculated value of the individual synaptic amplitudes
## X.ahp: output spike train convolved with the basis functions for simulating the after-spike currents
##
## returns: either 
## - v_soma, if vv is not provided; or
## - the error between vv and its predicion, + regularisation terms

	N <- nrow(X); Ne <- length(unlist(pars$Wc.e)); Ni <- length(unlist(pars$Wc.i))
	L <- ncol(X)
	Tmax <- L / dt
	M <- length(pars$Jc)
	if (!is.null(vv)) {
		if (length(vv)!=L) stop("The training signal and the input must have the same length!")
	}

	if (!is.null(regpars) & is.null(vv)) stop('to calculate regularisation error vv must be provided')

	if (!("delay.t" %in% names(pars))) {
		if (logpars$Tau) pars$delay.t <- log(rep(1/10, M)) else pars$delay.t <- rep(1/10, M)
	}
	
	if (!is.null(X.ahp)){
		if (!("W.ahp" %in% names(pars))) {
			warning ("W.ahp must be provided for caluclating after-spike currents. X.ahp will be ignored.")
			X.ahp <- NULL
		} else {
			if (length(pars$W.ahp) != nrow(X.ahp)) {
				warning("Number of basis parameters does not match the number of basis functions.  X.ahp will be ignored.")
				X.ahp <- NULL
				pars$W.ahp <- NULL
			}
		}
	}
	
	if ("Ww.e2" %in% names(pars)) {
		if (compare.lists(pars$Ww.e2, pars$Ww.e)) { ## valid Ww.e2 found
			if (!double) {
				warning('double synapses found, double is set to TRUE')	
				double <- T
			}
		} else { ## Ww.e2 is not valid
			if (double){
				warning('double synapses are not valid, double is set to FALSE')	
				double <- F								
			}
		}
	} else {
		pars['Ww.e2'] <- list(NULL)
		double <- F
	}

	if (!alphasyn & ("dTau.e" %in% names(pars))) {
		if (compare.lists(pars$dTau.e, pars$Tau.e)) { ## valid dTau.e found
			if (double) {
				warning('double synapses found so alpha synapses will be used')	
				alphasyn <- T
			} else {
				if (!is.null(regpars)) {
					stop('error: regularization is only implemented with alpha synapses')
				}
				# warning('double exponential synapses are detected, alphasyn is set to FALSE')	
				# alphasyn <- F
			}
		} else { ## dTau.e is not valid
			warning('seems that double exponential synapses are not valid, alphasyn is set to TRUE')	
			alphasyn <- T
		}
	} else {
		pars['dTau.e'] <- list(NULL)
		alphasyn <- T
	}
	
	#######################################################################
	resp <- with(pars, {		
		# v0 <- pars$v0; Jc <- pars$Jc; Jw <- pars$Jw; Wc.e <- pars$Wc.e; Wc.i <- pars$Wc.i; Ww.e <- pars$Ww.e; Ww.e2 <- pars$Ww.e2; Ww.i <- pars$Ww.i; Th <- pars$Th; Tau.e <- pars$Tau.e; Tau.i <- pars$Tau.i; dTau.e <- pars$dTau.e; delay.t <- pars$delay.t
		
		test.pars(v0, Jc, Jw, Wc.e, Wc.i, Ww.e, Ww.e2, Ww.i, Th, Tau.e, Tau.i, dTau.e, delay.t, N, verb=F)	
	
		## we will mostly constrain Jw and Tau.e to be positive. The simplest way to achieve this is to define them by their log
		## in this case we have to transform them back to the noral form before processing
		if (logpars$Jw) Jw <- exp(Jw)
		Jw[is.na(Th)] <- 1 ## for linear subunits Jw is redundant with Ww - so we set it to 1

		if (logpars$Tau) {
			delay.t  <- exp(delay.t)
			for (m in 1:M){
				if (!is.null(Tau.e[[m]])) Tau.e[[m]] <- exp(Tau.e[[m]])
				if (!is.null(Tau.i[[m]])) Tau.i[[m]] <- exp(Tau.i[[m]])
			}
			if (!alphasyn){
				for (m in 1:M){
					if (!is.null(dTau.e[[m]])) dTau.e[[m]] <- exp(dTau.e[[m]])
				}				
			}
		}
		if (logpars$W) {
			for (m in 1:M){
				if (!is.null(Ww.e[[m]])) Ww.e[[m]] <- exp(Ww.e[[m]])
				if (!is.null(Ww.i[[m]])) Ww.i[[m]] <- (-1) * exp(Ww.i[[m]])
				if (double) {
					if (!is.null(Ww.e2[[m]])) Ww.e2[[m]] <- exp(Ww.e2[[m]])
				}
			}
		}
		
		## subunit gains and amplitudes required for regularisation
		if (!is.null(regpars) | (calc.minis)) {
			a.e <- unlist(Ww.e) * exp(-1)
			if (double) a.e2 <- unlist(Ww.e2) * exp(-1)
			a.i <- unlist(Ww.i) * exp(-1)
			gain.subunits <- rep(1, M) # c_j f_j
			baseline.subunits <- rep(0, M) 
			total.gain.subunits <- rep(1, M) # c_j f_j \prod_k c_k f_k
		}

		##################################################################
		### first, calculate the baseline for all subunits - we have to start from the leaves
		if (!is.null(regpars) | (calc.minis)) {
			subunits.done <- 0 # leaves already processed
			Jc.orig <- Jc
	
			while (length(subunits.done) < (M+1)){ 
			## we repeat this until all the subunits ar eprcessed ...
				i.leaves <- !seq(1,M) %in% Jc ## leaves are those subunits that do not receive input from other subunits
				i.remain <- !seq(1,M) %in% subunits.done ## we only need the remaining leaves!
				ii.leaves <- seq(1,M)[i.leaves & i.remain]
				if (length(ii.leaves)==0) {
					cat("No further leaves found, the graph may contain a loop among the following subunits: ", seq(1,M)[i.remain], " ")
					stop(" We will quit.")
				}
				for (ii in ii.leaves){
					ii.children <- which(Jc.orig == ii)
					if (length(ii.children) > 0){
						for (ii.child in ii.children) 	baseline.subunits[ii] <- baseline.subunits[ii] + Jw[ii.child] *sigm(baseline.subunits[ii.child], c=Th[ii.child])
					}
					Jc[ii] <- NA ## The subunit is already processed - does not give further input
					subunits.done <- c(subunits.done, ii) ## register the subunits already processed				
				}
			}
			Jc <- Jc.orig
		}


		##########################################
		## second, we calculate the synaptic input to each dendritic branch
		##########################################

		Y <- matrix(0, M, L)
		N.esyns.m <- 1
		N.isyns.m <- 1
		
		for (m in 1:M){
			# subunit gain = Jw[m] * f'(0) - evaluated in the absence of inputs
			if (!is.null(regpars) | (calc.minis)) {
				gain.subunits[m] <- Jw[m]
				if (!is.na(Th[m])) {
					ii.children <- (Jc == m)
					gain.subunits[m] <- gain.subunits[m] * sigm(baseline.subunits[m], c=Th[m], d=1)			
				}
	
				## Connectivity is defined from the root - soma - towards the leaves. 
				## So the parent has already been processed when porcessing the leaves!
				total.gain.subunits[m] <- gain.subunits[m]
				i.parent <- Jc[m]
				if (i.parent != 0) total.gain.subunits[m] <- total.gain.subunits[m] * total.gain.subunits[i.parent]
			}

			## calculate the synaptic input to each dendritic branch
			if (length(Wc.e[[m]]) > 0){
				for (i.syn in 1:length(Wc.e[[m]])) {
					# cat(m, Tau.e[[m]][i.syn], "\n")
					if (alphasyn){
						Y[m,] <- Y[m,] + int.spikes(X, dt, Wc.e[[m]][i.syn], Ww.e[[m]][i.syn], Tau.e[[m]][i.syn], delay.t[m])
						if (double) Y[m,] <- Y[m,] + int.spikes(X, dt, Wc.e[[m]][i.syn], Ww.e2[[m]][i.syn], 10.4+scale.tau2*Tau.e[[m]][i.syn], delay.t[m])
					} else {
						Y[m,] <- Y[m,] + int.spikes(X, dt, Wc.e[[m]][i.syn], Ww.e[[m]][i.syn], Tau.e[[m]][i.syn], delay.t[m], dTau=dTau.e[[m]][i.syn])
					}
				}
				
				## just take the next k where k <- length(Wc.e[[m]])!
				## Wc is indexing the presynaptic spike train and not the postsynaptic parameters / subunits!
				
				if (!is.null(regpars) | (calc.minis)) {
					n.esyns.m <- length(Wc.e[[m]])
					a.e[N.esyns.m:(N.esyns.m + n.esyns.m-1)] <- a.e[N.esyns.m:(N.esyns.m + n.esyns.m-1)] * total.gain.subunits[m]
					if (double) a.e2[N.esyns.m:(N.esyns.m + n.esyns.m-1)] <- a.e2[N.esyns.m:(N.esyns.m + n.esyns.m-1)] * total.gain.subunits[m]
					N.esyns.m <- N.esyns.m + n.esyns.m
				}
			}

			if (length(Wc.i[[m]]) > 0) {
				for (i.syn in 1:length(Wc.i[[m]])) Y[m,] <- Y[m,] + int.spikes(X, dt, Wc.i[[m]][i.syn], Ww.i[[m]][i.syn], Tau.i[[m]][i.syn], delay.t[m])
				if (!is.null(regpars) | (calc.minis)) {
					n.isyns.m <- length(Wc.i[[m]])
					a.i[N.isyns.m:(N.isyns.m + n.isyns.m-1)] <- a.i[N.isyns.m:(N.isyns.m + n.isyns.m-1)] * total.gain.subunits[m]
					N.isyns.m <- N.isyns.m + n.isyns.m
				}
			}
		}

		## Add the after-spike currents!
		if (!is.null(X.ahp)) Y[1,] <- Y[1,] + W.ahp %*% X.ahp

		##########################################
		## next, we start from the leaves and apply the nonlinearities as well as add the inputs from the children
		## 1. Where are the leaves?
		R <- matrix(0, M, L) # a matrix with the activation of the subunits
		subunits.done <- 0 # leaves already processed
		i.soma <- NULL # the index of the root subunit which is the soma
		while (length(subunits.done) < (M+1)){ 
		## we repeat this until all the subunits ar eprcessed ...
			i.leaves <- !seq(1,M) %in% Jc ## leaves are those subunits that do not receive input from other subunits
			i.remain <- !seq(1,M) %in% subunits.done ## we only need the remaining leaves!
			ii.leaves <- seq(1,M)[i.leaves & i.remain]
			if (length(ii.leaves)==0) {
				cat("No further leaves found, the graph may contain a loop among the following subunits: ", seq(1,M)[i.remain], " ")
				stop(" We will quit.")
			}
			for (ii in ii.leaves){
			# 2. apply the sigmoidal nonlinearity for every leaves
				if (is.na(Th[ii])) {
					if (ii !=1) warning(paste("Subunits are by definition nonlinear. All linear subunits should be merged into their parents. Subunit", ii, "seems to be linear."))
					R[ii,] <- Y[ii,] 
				} else R[ii,] <- sigm(Y[ii,], c=Th[ii])
	
			# 3. add the input from the child to the parent
				ii.parent <- Jc[ii]
				if (ii.parent>0){
					Y[ii.parent,] <- Y[ii.parent,] + Jw[ii] * R[ii,]
				} else {
					if (ii != 1) warning(paste("The soma should be the first subunit. Now we found that subunit ", ii, "has no parent!"))
					if (!is.null(i.soma)) stop(paste("This neuron has more than one some! Subunits", ii, "and", i.soma, "are both root subunits!"))
					i.soma <- ii
				}
				Jc[ii] <- NA ## The subunit is already processed - does not give further input
				subunits.done <- c(subunits.done, ii) ## register the subunits already processed				
			}
		}
		v <- Jw[i.soma] * R[i.soma,] + v0
		if (verbose) cat("somatic gain:", Jw[i.soma], "\n")

		if (is.null(vv)) {
			out <- v
			if (calc.minis) {
				out <- list(v=v, a.e=a.e, a.i=a.i)
				if (Ne > 0) out$logTau.e <- log(unlist(Tau.e))
				if (Ni > 0) out$logTau.i <- log(unlist(Tau.i)) else {
					out$logTau.i <- list(NULL)
					out$a.i <- list(NULL)
				}
				if (double) out$a.e2 <- a.e2
			}
		} else {
			ind.vv <- !is.na(vv) # we only check where the reference is defined
			err <- (v - vv)[ind.vv]
			out <- sum(err^2, na.rm=T) # NEGATIVE log likelihood - this is what we want to minimise!; we omit the variance of the noise here
			if (sum(is.na(v)) > 0) out <- var(vv[ind.vv]) * length(vv[ind.vv]) # if NA is predicted then we return a theoretical maximum
			
			if (!is.null(regpars)) {
				if (verbose) cat("error without regularisation:", out, "\n\n")
				# Ww - lognormal prior: the variance is set by the variance of the minis, but Ww is much smaller in vivo - so the variance is also smaller!
				# separate variance for Ww.e1, Ww.e2 and Ww.i
				out <- out + regpars$alpha.Ww * sum((log(a.e) - log(regpars$Ww.e1))^2)
				if (verbose) cat("error with a.e:", out, "\n")
				if (verbose) cat("Ww.e1:", a.e[1:5], "\n regpars WW.e1:", regpars$Ww.e1[1:5], "\n\n")
				if (double) {
					out <- out + regpars$alpha.Ww * sum((log(a.e2) - log(regpars$Ww.e2))^2)
					if (verbose) cat("error with a.e2:", out, "\n")
					if (verbose) cat("Ww.e2:", a.e2[1:5],  "\n regpars WW.e2:", regpars$Ww.e2[1:5], "\n\n")
				}
				if (Ni > 0){
					out <- out + regpars$alpha.Ww * sum((log(abs(a.i)) - log(abs(regpars$Ww.i)))^2)
					if (verbose) cat("error with a.i:", out, "\n")
					if (verbose) cat("Ww.i:", a.i[1:5], "\n regpars WW.i:", regpars$Ww.i[1:5], "\n\n")
				}
				out <- out + regpars$alpha.Tau * sum((log(unlist(Tau.e)) - regpars$logTau.e)^2)
				if (verbose) cat("error with Tau.e:", out, "\n")
				if (verbose) cat("Tau.e:", log(unlist(Tau.e)[1:5]), "\n regpars Tau.e: ", regpars$logTau.e[1:5], "\n\n")
				
				if (Ni > 0){
					out <- out + regpars$alpha.Tau * sum((log(unlist(Tau.i)) - regpars$logTau.i)^2)
					if (verbose) cat("error with Tau.i:", out, "\n")
					if (verbose) cat("Tau.i:", log(unlist(Tau.i)[1:5]), "\n regpars Tau.i: ",  regpars$logTau.i[1:5], "\n\n")
				}
			}	
		}
		out
		})
	resp
}
