## A script to calculate the gradient of an hGLM model's response wrt its parameters
## AHP is conditional on the spike times - provided by the matrix X.ahp, which contains the basis functions convolved with the output spikes
## We need it to be efficient to be able to evaluate it inside the optimization
## To test this function, run 
#    source('./UnitTests_Grad.R', chdir = TRUE)

source('./Sim_hGLM_sub.R')

# the function must handle regularisation and double-kernel synapses, gradient wrt. ahp parameters

grad.hGLM <- function(X, dt, pars, vv=NULL, logpars=list(Jw=T, Tau=T, W=F), regpars=NULL, double=F, alphasyn=T, scale.tau2=2.8, verbose=F, response=F, X.ahp=NULL){
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
	# delay.t: optional vector of the synaptic delays in ms (default: 1/10 ms) - NO GRADIENT IS CALCULATED wrt DELAY
## vv: response, if error is calculated
## logpars: parameters Tau, Jw and Ww are strictly positive, so they can be defined by their log
## regpars: NULL (no regularisation) or list with the prior mean and inverse variance for the log of Ww and Tau
	# Ww.e, (Ww.e2), Ww.i, Tau.e, Tau.i, alpha.Ww, alpha.Tau
## double: logical; simple or double alpha kernels 
## scale.tau2: if double kernels are used, the scaling factor between their time constants
## verbose: regularisation details are provided
## response: indicates whether the response should be provided
## X.ahp: output spike train convolved with the basis functions for simulating the after-spike currents

# dt <- 1; pars <-  newpars; vv <- newpars$v; double=T; logpars <- list(Jw=T, Tau=T, W=F); scale.tau2 <- 2.8; verbose <- 2; response <- F
# X <- args$inp; dt=args$dt; pars=ppars; vv=args$v; logpars=list(Jw=T, Tau=T, W=F); regpars=args$regpars; double=args$double; scale.tau2=args$scale; verbose=F; response=F; X.ahp=args$X


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
			n.basis <- length(pars$W.ahp)
		}
	}

	if ("Ww.e2" %in% names(pars)) {
		if (compare.lists(pars$Ww.e2, pars$Ww.e)) { ## valid Ww.e2 found
			if (!double) {
				warning('seems that double synapses were used, double is set to TRUE')	
				double <- T
			}
		} else { ## Ww.e2 is not valid
			if (double){
				warning('seems that double synapses are not valid, double is set to FALSE')	
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
	
	add.regularisation <- !is.null(regpars)

	if (add.regularisation & logpars$W) stop('error: regularisation is not implemented for logarithmic Ww parameters.')
	
	
	#######################################################################
	resp <- with(pars, {		
		# v0 <- pars$v0; Jc <- pars$Jc; Jw <- pars$Jw; Wc.e <- pars$Wc.e; Wc.i <- pars$Wc.i; Ww.e <- pars$Ww.e; Ww.e2 <- pars$Ww.e2; Ww.i <- pars$Ww.i; Th <- pars$Th; Tau.e <- pars$Tau.e; Tau.i <- pars$Tau.i; dTau.e <- pars$dTau.e; delay.t <- pars$delay.t
		
		test.pars(v0, Jc, Jw, Wc.e, Wc.i, Ww.e, Ww.e2, Ww.i, Th, Tau.e, Tau.i, dTau.e, delay.t, N, verb=F)	
	
		## we will mostly constrain Jw and Tau.e to be positive. The simplest way to achieve this is to define them by their log
		## in this case we have to transform them back to the normal form before processing
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
		if (add.regularisation) {
			a.e <- unlist(Ww.e) * exp(-1)
			if (double) a.e2 <- unlist(Ww.e2) * exp(-1)
			a.i <- unlist(Ww.i) * exp(-1)
			gain.subunits <- rep(1, M) # c_j f_j
			total.gain.subunits <- rep(1, M) # c_j f_j \prod_k c_k f_k
			baseline.subunits <- rep(0, M) #
			children <- list(); for (m in 1:M) children[[m]] <- m			
			dloggain.dtheta.subunits <- rep(0, M) 
		}
	
		## the algorithm does not work, if the subunits are linear!
		if (sum(is.na(Th[-1])) > 0) stop("Gradients are calculated with the assumption that subunits are nonlinear!")

		##################################################################
		### first, calculate the baseline for all subunits - we have to start from the leaves
		if (add.regularisation) {
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
		
		Y <- matrix(0, M, L) # inputs to the subunits
		esyn.subunit <- rep(NA, Ne) # a vector with the index of the subunit a particular synapse belongs to

		if (Ne > 0){
			phi.e <- matrix(0, Ne, L) # a matrix for the excitatory inputs
			d.phi.e <- matrix(0, Ne, L) # a matrix for the derivative of the excitatory inputs wrt tau
	
			if (double){
				phi.e2 <- matrix(0, Ne, L) # a matrix for the SLOW excitatory
				d.phi.e2 <- matrix(0, Ne, L) # a matrix for the derivative of the SLOW excitatory inputs wrt tau
			}
			
			if (!alphasyn) d.dphi.e <- matrix(0, Ne, L)
		}

		N.esyns.m <- 1

		if (Ni > 0){
			phi.i <- matrix(0, Ni, L) # a matrix for the inhibitory inputs to the subunits
			d.phi.i <- matrix(0, Ni, L) # a matrix for the derivative of the inhibitory inputs wrt tau
			isyn.subunit <- rep(NA, Ni)
			N.isyns.m <- 1
		}

		for (m in 1:M){
			# subunit gain = Jw[m] * f'(0) - evaluated in the absence of inputs
			if (add.regularisation) {
				gain.subunits[m] <- Jw[m]
				if (!is.na(Th[m])) {
					ii.children <- (Jc == m)
					# baseline <- sum(Jw[ii.children])/2
					gain.subunits[m] <- gain.subunits[m] * sigm(baseline.subunits[m], c=Th[m], d=1)			
				}
	
				## Connectivity is defined from the root - soma - towards the leaves. 
				## So the parent has already been processed when porcessing the leaves!
				total.gain.subunits[m] <- gain.subunits[m]
				i.parent <- Jc[m]
				if (i.parent != 0) total.gain.subunits[m] <- total.gain.subunits[m] * total.gain.subunits[i.parent]
				if (!is.na(Th[m])) {
					ii.children <- (Jc == m)
					# baseline <- sum(Jw[ii.children])/2
					dloggain.dtheta.subunits[m] <- 2*sigm(baseline.subunits[m], c=Th[m]) - 1 # r * (1-r) * (2r-1)/(r(1-r))
				}
			}

			if (length(Wc.e[[m]]) > 0){
				# print("Tau.e and gradients wrt. Tau.e  ...")
				for (i.syn in 1:length(Wc.e[[m]])) {
					ie <- Wc.e[[m]][i.syn]
					if (alphasyn){
						phi.e[ie,] <- int.spikes(X, dt, ie, Ww=1, Tau.e[[m]][i.syn], delay.t[m])
						d.phi.e[ie,] <- int.spikes(X, dt, ie, Ww=1, Tau.e[[m]][i.syn], delay.t[m], grad.Tau = T)
						Y[m,] <- Y[m,] + Ww.e[[m]][i.syn] * phi.e[ie,]
						if (double){
							phi.e2[ie,] <- int.spikes(X, dt, ie, Ww=1, 10.4 + scale.tau2*Tau.e[[m]][i.syn], delay.t[m])
							d.phi.e2[ie,] <- scale.tau2 * int.spikes(X, dt, ie, Ww=1, 10.4 + scale.tau2*Tau.e[[m]][i.syn], delay.t[m], grad.Tau = T) # we still need the derivative wrt Tau - and not scale.tau2*Tau!
							Y[m,] <- Y[m,] + Ww.e2[[m]][i.syn] * phi.e2[ie,]
						}
					} else {
						phi.e[ie,] <- int.spikes(X, dt, ie, Ww=1, Tau.e[[m]][i.syn], delay.t[m], dTau=dTau.e[[m]][i.syn])
						d.phi.e[ie,] <- int.spikes(X, dt, ie, Ww=1, Tau.e[[m]][i.syn], delay.t[m], grad.Tau = 1, dTau=dTau.e[[m]][i.syn])
						d.dphi.e[ie,] <- int.spikes(X, dt, ie, Ww=1, Tau.e[[m]][i.syn], delay.t[m], grad.Tau = 2, dTau=dTau.e[[m]][i.syn])
						Y[m,] <- Y[m,] + Ww.e[[m]][i.syn] * phi.e[ie,]
					}
					# esyn.subunit[ie] <- m
				}
				
				n.esyns.m <- length(Wc.e[[m]])
				if (add.regularisation) {
					a.e[N.esyns.m:(N.esyns.m + n.esyns.m-1)] <- a.e[N.esyns.m:(N.esyns.m + n.esyns.m-1)] * total.gain.subunits[m]
					if (double) a.e2[N.esyns.m:(N.esyns.m + n.esyns.m-1)] <- a.e2[N.esyns.m:(N.esyns.m + n.esyns.m-1)] * total.gain.subunits[m]
				}
				esyn.subunit[N.esyns.m:(N.esyns.m + n.esyns.m-1)] <- m
				N.esyns.m <- N.esyns.m + n.esyns.m
			}
			
			if (length(Wc.i[[m]]) > 0){
				for (i.syn in 1:length(Wc.i[[m]])) {
					# print("Tau.i and gradients wrt. Tau.i ...")
					ii <- Wc.i[[m]][i.syn]
					phi.i[ii-Ne,] <- int.spikes(X, dt, ii, Ww=1, Tau.i[[m]][i.syn], delay.t[m])
					d.phi.i[ii-Ne,] <- int.spikes(X, dt, ii, Ww=1, Tau.i[[m]][i.syn], delay.t[m], grad.Tau=T)
					Y[m,] <- Y[m,] + Ww.i[[m]][i.syn] * phi.i[ii-Ne,]
					# isyn.subunit[ii-Ne] <- m
				}

				n.isyns.m <- length(Wc.i[[m]])
				if (add.regularisation) {
					a.i[N.isyns.m:(N.isyns.m + n.isyns.m-1)] <- a.i[N.isyns.m:(N.isyns.m + n.isyns.m-1)] * total.gain.subunits[m]
				}
				isyn.subunit[N.isyns.m:(N.isyns.m + n.isyns.m-1)] <- m
				N.isyns.m <- N.isyns.m + n.isyns.m
			}
		}
	
		# matplot(t(phi.e), t="l")
		
		## Add the after-spike currents!
		if (!is.null(X.ahp)) Y[1,] <- Y[1,] + W.ahp %*% X.ahp

		##########################################
		## THE DOWNWARD PASS - from leaves to the root
		## next, we start from the leaves and apply the nonlinearities as well as add the inputs from the children
		Jc.orig <- Jc
		## 1. Where are the leaves?
		R <- matrix(0, M, L) # a matrix with the activation of the subunits
		subunits.done <- 0 # leaves already processed
		i.soma <- NULL # the index of the root subunit which is the soma
		while (length(subunits.done) < (M+1)){ 
		## we repeat this until all the subunits are porcessed ...
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
					if (ii !=1) stop(paste("Subunits are by definition nonlinear. All linear subunits should be merged into their parents. Subunit", ii, "seems to be linear. Gradients are not implemented with linear subunits."))
					R[ii,] <- Y[ii,]
				} else R[ii,] <- sigm(Y[ii,], c=Th[ii])
	
			# 3. ad the input from the child to the parent
				ii.parent <- Jc[ii]
				if (ii.parent>0){
					Y[ii.parent,] <- Y[ii.parent,] + Jw[ii] * R[ii,]
				} else {
					if (ii != 1) warning(paste("The soma should be the first subunit. Now we found that subunit ", ii, "has no parent!"))
					if (!is.null(i.soma)) stop(paste("This neuron has more than one soma! Subunits", ii, "and", i.soma, "are both root subunits!"))
					i.soma <- ii
				}

				if (add.regularisation) if (Jc[[ii]] > 0) 	children[[Jc[ii]]] <- union(children[[ii]], children[[Jc[ii]]])
				Jc[ii] <- NA ## The subunit is already processed - does not give further input
				subunits.done <- c(subunits.done, ii) ## register the subunits already processed				
			}
		}
		v <- Jw[i.soma] * R[i.soma,] + v0
		
		##########################################
		## THE UPWARD PASS - from the root towards the leaves
		rho <- matrix(1, M, L)
		Jc <- Jc.orig
		
		roots <- 0 # the current roots are those subunit who innervate subunit 0 - soma
		i.root <- which(Jc %in% roots) # which subunits are currently root subunits?
		## the original roots have rho=1 - it is done! So find the second layer of roots!
		roots <- i.root
		i.root <- which(Jc %in% roots) # which subunits are currently root subunits?
		
	
		subunits.done <- 1 # subunits already processed - only the soma currently
		while (length(subunits.done) < (M)){ 
			for (ii in i.root){
				i.parent <- Jc[ii]
				if ((is.na(Th[1])) & (i.parent==1)) {
					## if the soma is linear: 
					rho[ii,] <- Jw[i.parent] * rho[i.parent,]
				} else {
					rho[ii,] <- R[i.parent,] * (1-R[i.parent,]) * Jw[i.parent] * rho[i.parent,]
				}
				subunits.done <- c(subunits.done, ii)
			}
			roots <- i.root
			i.root <- which(Jc %in% roots) # which subunits are currently root subunits?
		}

		##########################################
		## Calculate the gradients of the regularisation

		if (add.regularisation){
			## the difference between the actual and the prior mean
			err.Ww.e1 <- log(a.e) - log(regpars$Ww.e1)
			if (double) err.Ww.e2 <- log(a.e2) - log(regpars$Ww.e2)
			if (Ni>0) err.Ww.i <- log(-1*(a.i)) - log(-1*(regpars$Ww.i))
			err.Tau.e <- log(unlist(Tau.e)) - regpars$logTau.e
			if (Ni>0) err.Tau.i <- log(unlist(Tau.i)) - regpars$logTau.i

			### \sum_i w_ij error_ij - the total contribution of a subunit to the amplitude
			sum.err <- rep(0, M)
			if (Ne > 0){
				for (j in 1:Ne){
					esub <- esyn.subunit[j]
					sum.err[esub] <- sum.err[esub] + err.Ww.e1[j]
					if (double) sum.err[esub] <- sum.err[esub] + err.Ww.e2[j]
				}
			}
			if (Ni>0) {
				for (j in 1:Ni){
					isub <- isyn.subunit[j]
					sum.err[isub] <- sum.err[isub] + err.Ww.i[j]
				}
			}

			## vectors to store the gradients of the regularisation terms
			de.dj.reg <- rep(0, M)
			de.dth.reg <- rep(0, M)
			## regularisation gradients wrt J and th are calculated here
			for (m in 1:M){
				ch <- children[[m]]
				for (ch.m in ch){
					de.dj.reg[m] <- de.dj.reg[m] + regpars$alpha.Ww * 2 / Jw[m] * sum.err[ch.m]
					de.dth.reg[m] <- de.dth.reg[m] + regpars$alpha.Ww * 2 * sum.err[ch.m] * dloggain.dtheta.subunits[m]
				}
			}
	
			de.dw.e1.reg <- rep(NA, Ne)
			if (double) de.dw.e2.reg <- rep(NA, Ne)
			if (Ni>0) de.dw.i.reg <- rep(NA, Ni)
			
			de.dTau.e.reg <- 2 * regpars$alpha.Tau * err.Tau.e
			if (Ni>0) de.dTau.i.reg <- 2 * regpars$alpha.Tau * err.Tau.i
		}


		
		##########################################
		## And finally we have to calculate the gradients
		dv.dj <- matrix(NA, M, L)
		dv.dth <- matrix(NA, M, L)

		dv.dw.e1 <- matrix(NA, Ne, L)
		if (double) dv.dw.e2 <- matrix(NA, Ne, L)
		dv.dtau.e <- matrix(NA, Ne, L)
		if (!alphasyn) dv.ddtau.e <- matrix(NA, Ne, L)

		if (Ni>0) {
			dv.dw.i <- matrix(NA, Ni, L)
			dv.dtau.i <- matrix(NA, Ni, L)
		}

			
		## for each subunit, the gradient wrt. J and Th
		for (j in 1:M){
			dv.dj[j,] <- R[j,] * rho[j,]
			dv.dth[j,] <- - Jw[j] * R[j,] * (1-R[j,]) * rho[j,]
		}

		if (logpars$Jw) factor.Jw <- Jw else factor.Jw <- 1
		if (logpars$Tau) {
			factor.Tau.e <- unlist(Tau.e)
			factor.Tau.i <- unlist(Tau.i)
			if (!alphasyn) factor.dTau <- unlist(dTau.e)
		} else {
			factor.Tau.e <- 1
			factor.Tau.i <- 1			
			if (!alphasyn) factor.dTau <- 1
		}
		# if (logpars$W) {
			# factor.Ww.e <- unlist(Ww.e) 
			# factor.Ww.i <- (1) * unlist(Ww.i)  # don't need the -1!
			# if (double) factor.Ww.e2 <- unlist(Ww.e2) 
		# } else {
			# factor.Ww.e <- rep(1, length(unlist(Ww.e)))
			# factor.Ww.i <- rep(1, length(unlist(Ww.i)))
			# if (double) factor.Ww.e2 <- rep(1, length(unlist(Ww.e2)))
		# }
		
		## for each synapse, the gradient wrt. Ww and Tau
		if (Ne > 0){
			for (j in 1:Ne){
				esub <- esyn.subunit[j]
				dv.dw.e1[j,] <- phi.e[j,] * Jw[esub] * R[esub,] * (1-R[esub,]) * rho[esub,]
				dv.dtau.e[j,] <- unlist(Ww.e)[j] * d.phi.e[j,] * Jw[esub] * R[esub,] * (1-R[esub,]) * rho[esub,]
				if (double) {
					dv.dw.e2[j,] <- phi.e2[j,] * Jw[esub] * R[esub,] * (1-R[esub,]) * rho[esub,]
					dv.dtau.e[j,] <- dv.dtau.e[j,] + unlist(Ww.e2)[j] * d.phi.e2[j,] * Jw[esub] * R[esub,] * (1-R[esub,]) * rho[esub,] # unlisting works, as Ww.e has the same shape as Wc.e, which was used to define j above, around line 80
				} 
				if (!alphasyn){
					dv.ddtau.e[j,] <- unlist(Ww.e)[j] * d.dphi.e[j,] * Jw[esub] * R[esub,] * (1-R[esub,]) * rho[esub,]
				}				
				if (add.regularisation){
					de.dw.e1.reg[j] <- regpars$alpha.Ww * 2 * err.Ww.e1[j] / unlist(Ww.e)[j]
					if (double) de.dw.e2.reg[j] <- regpars$alpha.Ww * 2 * err.Ww.e2[j] / unlist(Ww.e2)[j]
				}
			}
		}

		if (Ni>0) {
			for (j in 1:Ni){
				isub <- isyn.subunit[j]
				dv.dw.i[j,] <- phi.i[j,] * Jw[isub] * R[isub,] * (1-R[isub,]) * rho[isub,]
				dv.dtau.i[j,] <- unlist(Ww.i)[j] * d.phi.i[j,] * Jw[isub] * R[isub,] * (1-R[isub,]) * rho[isub,]			
				if (add.regularisation) de.dw.i.reg[j] <- regpars$alpha.Ww * 2 * err.Ww.i[j]  / unlist(Ww.i)[j]
			}
		}

########################################################	
		if (!is.null(X.ahp)){
			dv.dw.ahp <- matrix(NA, n.basis, L)
			if (is.na(Th[1])) {
				for (j in 1:n.basis) dv.dw.ahp[j,] <- Jw[1] * X.ahp[j,]
			} else {
				for (j in 1:n.basis) dv.dw.ahp[j,] <- Jw[1] * R[1,] * (1-R[1,]) * X.ahp[j,]			
			}
		}
########################################################	

		## gradients for soma...
		if (is.na(Th[1])) {
		## if the soma is linear: 
			dv.dth[1,] <- 0
			dv.dj[1,] <- 0 ## for linear subunits Jw[1] is redundant with Ww.e[1] and Jw[ch[1]] - so we set it to 1

			if (add.regularisation){
				de.dj.reg[1] <- 0
				de.dth.reg[1] <- 0
			}
			
			esyn.soma <- which(esyn.subunit==1)
			if (length(esyn.soma) > 0){
				for (j in esyn.soma){
					dv.dw.e1[j,] <- phi.e[j,] * Jw[1] * rho[1,]
					dv.dtau.e[j,] <- unlist(Ww.e)[j] * d.phi.e[j,] * Jw[1] * rho[1,]
					if (double) {
						dv.dw.e2[j,] <- phi.e2[j,] * Jw[1] * rho[1,]
						dv.dtau.e[j,] <- dv.dtau.e[j,]  + unlist(Ww.e2)[j] * d.phi.e2[j,] * Jw[1] * rho[1,]
					} 
					if (!alphasyn) {
						dv.ddtau.e[j,] <- unlist(Ww.e)[j] * d.dphi.e[j,] * Jw[1] * rho[1,]
					}
				}
			}

			if (Ni>0) {
				isyn.soma <- which(isyn.subunit==1)
				if (length(isyn.soma) > 0){
					for (j in isyn.soma){
						dv.dw.i[j,] <- phi.i[j,] * Jw[1] * rho[1,]
						dv.dtau.i[j,] <- unlist(Ww.i)[j] * d.phi.i[j,] * Jw[1] * rho[1,]			
					}
				}
			}
		}
	

		## the output	
		if (is.null(vv)) {
		## if no training singal is available, the output is the derivative of the response wrt parameters
			out <- list(v=v, dj=dv.dj, dw.e=dv.dw.e1)
			if (double) out$dw.e2 <- dv.dw.e2
			out$dth <- dv.dth
			out$dtau.e <- dv.dtau.e
			if (!alphasyn) out$ddtau.e <- dv.ddtau.e
			if (Ni>0) {
				out$dw.i <- dv.dw.i
				out$dtau.i <- dv.dtau.i
			}
			if (!is.null(X.ahp)) out$dw.ahp <- dv.dw.ahp
		} else {
		## if we have the training signal the output is the derivative of the NEGATIVE log likelihood wrt the parameters
			ind.vv <- !is.na(vv) # we only check where the reference is defined
			err <- (v - vv)[ind.vv]
			K <- length(err)
			de.dv0 <- 2 * sum(err)
			de.dj <- 2 * colSums(matrix(err * t(dv.dj[,ind.vv]), K)) * factor.Jw # the factor is due to the fact that we opimise in the log domain
			de.dw.e1 <- 2 * colSums(matrix(err * t(dv.dw.e1[,ind.vv]), K)) #* factor.Ww.e
			if (double) de.dw.e2 <- 2 * colSums(matrix(err * t(dv.dw.e2[,ind.vv]), K)) #* factor.Ww.e2
			if (Ni>0) de.dw.i <- 2 * colSums(matrix(err * t(dv.dw.i[,ind.vv]), K)) #* factor.Ww.i
			de.dth <- 2 * colSums(matrix(err * t(dv.dth[,ind.vv]), K))
			de.dtau.e <- 2 * colSums(matrix(err * t(dv.dtau.e[,ind.vv]), K)) * factor.Tau.e
			if (!alphasyn) de.ddtau.e <- 2 * colSums(matrix(err * t(dv.ddtau.e[,ind.vv]), K)) * factor.dTau
			if (Ni>0) de.dtau.i <- 2 * colSums(matrix(err * t(dv.dtau.i[,ind.vv]), K)) * factor.Tau.i
			if (!is.null(X.ahp)) de.dw.ahp <- 2 * colSums(matrix(err * t(dv.dw.ahp[,ind.vv]), K)) 

			if (add.regularisation){
				de.dj <- de.dj + de.dj.reg * factor.Jw # the factor is due to the fact that we opimise in the log domain
				de.dth <- de.dth + de.dth.reg

				de.dw.e1 <- de.dw.e1 + de.dw.e1.reg
				if (double) de.dw.e2 <- de.dw.e2 + de.dw.e2.reg
				if (Ni>0) de.dw.i <- de.dw.i + de.dw.i.reg

				de.dtau.e <- de.dtau.e + de.dTau.e.reg
				if (Ni>0) de.dtau.i <- de.dtau.i + de.dTau.i.reg
			}
			out <- list(dv0=de.dv0, dj=de.dj, dw.e=de.dw.e1)
			if (double) out$dw.e2=de.dw.e2
			out$dth <- de.dth
			out$dtau.e <- de.dtau.e
			if (!alphasyn) out$ddtau.e <- de.ddtau.e
			if (Ni>0) out$dw.i <- de.dw.i
			if (Ni>0) out$dtau.i <- de.dtau.i
			if (!is.null(X.ahp)) out$dw.ahp <- de.dw.ahp
		}
		## return the response if required
		if (response) out$v <- v		
		out
	})
	resp
}
