
##############################################
## A function to randomlly intialise the parameters Ww.e, Jw and/or Th to generate data or to initialise learning
## We initialise these parameters randomly, though we try to be close to linear in general
## Tha output is standardised - its mean is 0 and its variance is 1
## To test this function, run 
#    source('./UnitTests_Init.R', chdir = TRUE)
## initialization for the parameters (including W.ahp) is formally correct, but this does not mean that it is biologically meaningful


source('../Utils/Utils_hGLMs.R')
source('../Utils/IntSpikes.R', chdir = TRUE)
source('../Utils/Sigm.R')
source('../Graphics/Graphics.R', chdir = TRUE)
# library(viridis)

init.hGLM <- function(X, dt, pars, logpars=list(Jw=T, Tau=T, W=F), double=F, alphasyn=T, scale.tau2=2.8, lin.soma=F, rescale.soma=T, mean.v=0, sd.v=1, rand.seed=12, nSD=4, X.ahp=NULL){
	# nSD parameter controls the linearity of the cell. The larger it is the more linear the cell becomes:
	# nSD * SD of the input will be in the [-1,1] range of the sigmoid
	# if initialising an 1L model, then nSD is the variance of the output!
	# rescale.soma: somatic variance is scaled only if rescale.soma==T
	# mean.v, sd.v:  - mean and variance of the somatic output
	# 		<- this is important if we want to match a 1L/1N model with complex hieararchy
	# dt <- 1; logpars <- list(Jw=T, Tau=T, W=F); double<-F; scale.tau2<-2.8; lin.soma<-F; rescale.soma<-T; rand.seed<-12; nSD<-1; X.ahp<-NULL; alphasyn <- T

	N <- nrow(X)
	L <- ncol(X)
	Tmax <- L / dt
	M <- length(pars$Jc)
	set.seed(rand.seed)	

	if (double){
		if (alphasyn == F) {
			warning('double exponential synapses are only implemented with single kernels, double is set to FALSE')
			double <- F
		}
	} 

	if (!("v0" %in% names(pars))) pars$v0 <- -50
	if (!("Jw" %in% names(pars))) {
		if (logpars$Jw) pars$Jw <- log(rep(4, M)) else pars$Jw <- rep(4, M)  # to compensate for the gain of the sigmoid ~0.25
	}
	if (!("Ww.e" %in% names(pars))){
		if (logpars$W) pars$Ww.e <- makeParList(template=pars$Wc.e, min=log(0.5), max=log(2)) else pars$Ww.e <- makeParList(template=pars$Wc.e, min=0.5, max=2)
	}
	if(double){
		if (!("Ww.e2" %in% names(pars))) {
			if (logpars$W) pars$Ww.e2 <- makeParList(template=pars$Wc.e, min=log(0.5), max=log(2)) else pars$Ww.e2 <- makeParList(template=pars$Wc.e, min=0.5, max=2)
		}
	} else {
		pars['Ww.e2'] <- list(NULL)
		# print(names(pars))
	}
	if (!("Ww.i" %in% names(pars))) {
		if (logpars$W) pars$Ww.i <- makeParList(template=pars$Wc.i, min=log(0.1), max=log(1)) else pars$Ww.i <- makeParList(template=pars$Wc.i, min=-1, max=0)
	}
	if (!("Th" %in% names(pars))) pars$Th <- rep(1, M)
	if (!("Tau.e" %in% names(pars))){
		if (logpars$Tau) pars$Tau.e <- makeParList(template=pars$Wc.e, min=log(2), max=log(8)) else pars$Tau.e <- makeParList(template=pars$Wc.e, min=2, max=8)
	}
	if (!("Tau.i" %in% names(pars))){
		if (logpars$Tau) pars$Tau.i <- makeParList(template=pars$Wc.i, min=log(2), max=log(30)) else pars$Tau.i <- makeParList(template=pars$Wc.i, min=2, max=30)
	}

	if (!alphasyn){
		if (!("dTau.e" %in% names(pars))){
			if (logpars$Tau) pars$dTau.e <- makeParList(template=pars$Wc.e, min=log(2), max=log(32)) else pars$dTau.e <- makeParList(template=pars$Wc.e, min=2, max=32)
		}		
	} else {
		pars['dTau.e'] <- list(NULL)
	}
	
	if (!("delay.t" %in% names(pars))) {
		if (logpars$Tau) pars$delay.t <- log(rep(1/10, M)) else pars$delay.t <- rep(1/10, M)
	}
	if (lin.soma) {
		pars$Th[1] <- NA 
	}

	add.ahp <- F
	if (!is.null(X.ahp)) {
		add.ahp <- T
		if ("W.ahp" %in% names(pars)) {
			if (length(pars$W.ahp) != nrow(X.ahp)) {
				warning("length of W.ahp does not match the number of basis functions")
				add.ahp <- F
				pars['W.ahp'] <- list(NULL)
			}
		} else {
			int.w <- rowSums(X.ahp) # the integral of each basis
			n.ahp <- nrow(X.ahp)
			pars$W.ahp <- rexp(n.ahp, 1) / int.w * sample(c(-1,1), n.ahp, replace=T) # longer basis - smaller weight
			pars$W.ahp <- 2 * pars$W.ahp / sum(pars$W.ahp)
		} 
	} else {
		pars['W.ahp'] <- list(NULL)
	}
		
	# print(pars)
	#######################################################################
	newpars <- with(pars, {		
	# v0 <- pars$v0; Jc <- pars$Jc; Jw <- pars$Jw; Wc.e <- pars$Wc.e; Wc.i <- pars$Wc.i; Ww.e <- pars$Ww.e; Ww.e2 <- pars$Ww.e2; Ww.i <- pars$Ww.i; Th <- pars$Th; Tau.e <- pars$Tau.e; Tau.i <- pars$Tau.i; delay.t <- pars$delay.t; W.ahp <- pars$W.ahp
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
	
		##########################################
		## if tau is too small the input is 0, we need to get it back to a meaningful range
		for (m in 1:M){
			if (!is.null(Tau.e[[m]])) {
				for (k in 1:length(Tau.e[[m]])) if (Tau.e[[m]][k] < 0.25) Tau.e[[m]][k] <- runif(1, 1, 5)
				if (!alphasyn) {
					for (k in 1:length(dTau.e[[m]])) if (dTau.e[[m]][k] < 0.25) dTau.e[[m]][k] <- runif(1, 1, 5)				
				}
			}
			if (!is.null(Tau.i[[m]])) {
				for (k in 1:length(Tau.i[[m]])) if (Tau.i[[m]][k] < 0.25) Tau.i[[m]][k] <- runif(1, 1, 5)
			}
		}
	
		##########################################
		## calculate the synaptic input to each dendritic branch
		##########################################
		
		Y <- matrix(0, M, L)
		for (m in 1:M){
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
			}

			if (length(Wc.i[[m]]) > 0) {
				for (i.syn in 1:length(Wc.i[[m]])) Y[m,] <- Y[m,] + int.spikes(X, dt, Wc.i[[m]][i.syn], Ww.i[[m]][i.syn], Tau.i[[m]][i.syn], delay.t[m])
			}
		}
	
		## Add the after-spike currents!
		if (add.ahp) Y[1,] <- Y[1,] + W.ahp %*% X.ahp
		offset.Y <- rowMeans(Y)	

		# matplot(t(Y), t="l")
		##########################################
		## next, we start from the leaves and apply the nonlinearities as well as add the inputs from the children
		## 1. Where are the leaves?
		subunits.done <- 0; Jc.orig <- Jc # leaves already processed
		while (length(subunits.done) < (M+1)){ 
		## we repeat this until all the leaves fall down...
			i.leaves <- !seq(1,M) %in% Jc ## leaves are those subunits that do not receive input from other subunits
			i.remain <- !seq(1,M) %in% subunits.done ## we only need the remaining leaves!
			ii.leaves <- seq(1,M)[i.leaves & i.remain]
			if (length(ii.leaves)==0) {
				cat("No further leaves found, the graph may contain a loop among the following subunits: ", seq(1,M)[i.remain], " ")
				stop(" We will quit.")
			}
			for (ii in ii.leaves){
			# 2. apply the sigmoidal nonlinearity for every leaves
				if (ii !=1){ # if not the soma
					if (is.na(Th[ii])) warning(paste("Subunits are by definition nonlinear. All linear subunits should be merged into their parents. Subunit", ii, "seems to be linear."))
					##############################################################
					## This part ensures that intermediate subunits will also be reasonably nonlinear
					## if the range of Y is too small - subunit is essentially linear - we rescale the Ww.e[ii] and Jw[children] parameters...
					## note, that Y now includes input from the children!
					range.Y <- sd(Y[ii,]) # to rescale inputs to the linear range
					rescale <- 1 / (nSD * range.Y)
					Ww.e[[ii]] <- rescale * Ww.e[[ii]]			
					if (double) Ww.e2[[ii]] <- rescale * Ww.e2[[ii]]				
					Ww.i[[ii]] <- rescale * Ww.i[[ii]]				

					## rescale input from the other subunits similarly
					ii.children <- (Jc.orig == ii)
					Jw[ii.children] <- rescale * Jw[ii.children]	
					Y[ii,] <- rescale * Y[ii,]
					## compensate in the output
					Jw[ii] <- Jw[ii] / rescale					

					Th[ii] <- mean(Y[ii,])
					R <- sigm(Y[ii,], c=Th[ii])
				} else { ## soma, which can remain nonlinear!
					if (rescale.soma){ # if we force the soma's output to have variance 1
						if (lin.soma) nSD.soma <- 1/sd.v else nSD.soma <- nSD 
							# linear: nSD.soma sets the output variance
							# nonlinear: nSD.soma sets the linearity of the soma
						range.Y <- sd(Y[ii,]) # to rescale inputs to the linear range
						rescale <- 1 / (nSD.soma * range.Y)
						Ww.e[[ii]] <- rescale * Ww.e[[ii]]				
						if (double) Ww.e2[[ii]] <- rescale * Ww.e2[[ii]]				
						Ww.i[[ii]] <- rescale * Ww.i[[ii]]				
						if (add.ahp) W.ahp <- rescale * W.ahp 
								
						## rescale input from the other subunits similarly
						ii.children <- (Jc.orig == ii)
						Jw[ii.children] <- rescale * Jw[ii.children]	
						Y[ii,] <- rescale * Y[ii,]						
						R <- Y[ii,]
						
						if (!lin.soma){
							Th[ii] <- mean(Y[ii,])
							R <- sigm(Y[ii,], c=Th[ii])
							## compensate in the output
							Jw[ii] <- sd.v / sd(R)
						}							
					} else { # keep the original variance
						if (lin.soma){
							ii.children <- (Jc.orig == ii)
							R <- Y[ii,]
						} else {
							ii.children <- (Jc.orig == ii)
							Th[ii] <- Th[ii] + sum(Jw[ii.children] / 2) - sum(offset.Y[-1])  # input from the children is offset by Jw/2 - Jw*Th/4
							R <- sigm(Y[ii,], c=Th[ii])
						}
					}
				}

				# add the input from the child to the parent
				ii.parent <- Jc[ii]
				if (ii.parent>0){ ## Not the soma...
					Y[ii.parent,] <- Y[ii.parent,] + Jw[ii] * R
				} else { ## The soma...
					if (is.na(Th[ii])) { # .. is linear
						Jw[ii] <- 1
						# v0 <- v0 - sum(Jw[ii.children] / 2) # v0 should compensate for the input offsets from the children!
						v0 <- v0 - sum(Jw[ii.children] / 2) + sum(offset.Y[-1])  # input from the children is offset by Jw/2 - Jw*Th/4
					} else {
						if (rescale.soma) v0 <- v0 - Jw[ii] / 2	# v0 should compensate for the input offsets due to the nonlinearity!	
					}
					v.soma <- Jw[ii] * R + v0
					if (ii != 1) warning(paste("The soma should be the first subunit. Now we found that subunit ", ii, "has no parent!"))
					i.soma <- ii
				}
				Jc[ii] <- NA ## The subunit is already processed - does not give further input
				subunits.done <- c(subunits.done, ii) ## register the subunits already processed				
			}
		}
		v <- v.soma
		
		if (rescale.soma){
			## scaling the mean and the variance of the output
			# m.v sd.v
			m.vv <- mean(v.soma)
			sd.vv <- sd(v.soma)
			if (abs(sd.v - sd.vv) > 1e-10) warning('sd is not set properly!')			
			v0 <- v0 - m.vv + mean.v
			v <- v - m.vv + mean.v
			if (abs(mean(v) - mean.v) > 1e-10) warning('mean is not set properly!')			
		}
		
		## if we use positivity constraints on Jw and Tau.e, we have to define them by their log
		if (logpars$Jw) Jw <- log(Jw)
		if (logpars$Tau) {
			for (m in 1:M){
				if (!is.null(Tau.e[[m]])) Tau.e[[m]] <- log(Tau.e[[m]])
				if (!is.null(Tau.i[[m]])) Tau.i[[m]] <- log(Tau.i[[m]])
			}
			if (!alphasyn){
				for (m in 1:M){
					if (!is.null(dTau.e[[m]])) dTau.e[[m]] <- log(dTau.e[[m]])
				}				
			}
			delay.t  <- log(delay.t)
		}
		if (logpars$W){
			for (m in 1:M){
				if (!is.null(Ww.e[[m]])) Ww.e[[m]] <- log(Ww.e[[m]])
				if (!is.null(Ww.i[[m]])) Ww.i[[m]] <- log((-1) * Ww.i[[m]])
				if (double) {
					if (!is.null(Ww.e2[[m]])) Ww.e2[[m]] <- log(Ww.e2[[m]])
				}
			}
		}

		# the ordering is important, as we will need unlist the parameters to create a vector for optimization
		np <- list(v=v, v0=v0, Jc=Jc.orig, Wc.e=Wc.e, Wc.i=Wc.i, Jw=Jw, Ww.e=Ww.e, Ww.e2=Ww.e2, Ww.i=Ww.i, Th=Th, Tau.e=Tau.e, dTau.e=dTau.e, Tau.i=Tau.i, delay.t=delay.t, W.ahp=W.ahp)
		np
	} )
	newpars
}

