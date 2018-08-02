## A script to simulate learning in a hLN model

## We define simple wrappers that allow the evaluation of the model as well as its gradient 
## These wrapper functions can be used by the optim() function to optimise parameters of the model
source('./Grad_hGLM.R')

##########################################################
## a wrapper to calculate the error from the parmeters and the arguments
err.hGLM <- function(pars.vec, args, ret.resp=F, ret.parlist=F){
	## pars.vec: a parameter vector for the parameters of the LN model
	## args: $inp is the input spike train
	## 		$v is the target output
	## 		$dt is the time resolution of the input
	## 		$Jc is the structure of the dendrite
	## 		$Wc.e is the connectivity of the dendrite (excitatory synapses)
	## 		$Wc.i is the connectivity of the dendrite (inhibitory synapses)
	##			$linear.soma whether the soma is linear or nonlinear
	##	 The following vector contains the parameters to be optimised
	## 		$par=c(v0=T, Jw=T, Ww.e=T, Ww.e2=T, Ww.i=T, Th=T, Tau.e=T, Tau.i=T, W.ahp=T) - which parameters are in the parameter vector
	##					only parameters specified here will be optimised. All F has to be present amongst the arguments
	##		$v0 = NULL - if ($par["v0"] == F)  this will be used to evaluate the hGLM model
	## 	$Jw = NULL
	##		$Ww.e = NULL
	##		$Ww.e2 = NULL
	##		$Ww.i = NULL
	##		$Th = NULL
	##		$Tau.e = NULL
	##		$Tau.i = NULL
	##		$W.ahp = NULL
	##
	## returns:
	## default: sum of squared error between the predicted and reference membrane potential
	## if (ret.resp): the predicted voltage vector is returned
	## if (ret.parlist): the list of parameters is returned - in a form that is compatible with sim.GLM

	# test.args(pars.vec, args)
	if (ret.resp & ret.parlist) {
		warning('can not return both the response and the parameter list - only the parameter list will be returned')
		ret.resp <- F
	}

	if (!args$alphasyn & args$double) stop('error: alphasyn is FALSE and double is TRUE - double exponential synapses are only implemented in single-kernel mode.')

	M <- length(args$Jc)
	Ne <- length(unlist(args$Wc.e))
	Ni <- length(unlist(args$Wc.i))
	if (!is.null(args$X.ahp)) n.basis <- nrow(args$X.ahp)

	ppars <- list(Jc=args$Jc, Wc.e=args$Wc.e, Wc.i=args$Wc.i)

	if (args$par["v0"] == T) {
			ppars$v0 <- pars.vec[1]; pars.vec <- pars.vec[-1]
		} else ppars$v0 <- args$v0

	if (args$par["Jw"] == T) {
			ppars$Jw <- pars.vec[1:M]; pars.vec <- pars.vec[-(1:M)]
		} else ppars$Jw <- args$Jw
	if (args$par["Ww.e"] == T) {
			ppars$Ww.e <- reshapeVec2List.no.perm(args$Wc.e, pars.vec[1:Ne]); pars.vec <- pars.vec[-(1:Ne)]
	} else ppars$Ww.e <- args$Ww.e
	if (args$double){
		if (args$par["Ww.e2"] == T) {
			ppars$Ww.e2 <- reshapeVec2List.no.perm(args$Wc.e, pars.vec[1:Ne]); pars.vec <- pars.vec[-(1:Ne)]
		} else ppars$Ww.e2 <- args$Ww.e2
	} else {
		ppars['Ww.e2'] <- list(NULL)
	}
	if (args$par["Ww.i"] == T) {
			ppars$Ww.i <- reshapeVec2List.no.perm(args$Wc.i, pars.vec[1:Ni]); pars.vec <- pars.vec[-(1:Ni)]
	} else ppars$Ww.i <- args$Ww.i
		
	if (args$par["Th"] == T) {
			ppars$Th <- pars.vec[1:M]; pars.vec<-pars.vec[-(1:M)]
		} else ppars$Th <- args$Th
		
	if (args$par["Tau.e"] == T) {
			ppars$Tau.e <- reshapeVec2List.no.perm(args$Wc.e, pars.vec[1:Ne]); pars.vec <- pars.vec[-(1:Ne)]
	} else ppars$Tau.e <- args$Tau.e
	if (!args$alphasyn){
		if (args$par["dTau.e"] == T) {
				ppars$dTau.e <- reshapeVec2List.no.perm(args$Wc.e, pars.vec[1:Ne]); pars.vec <- pars.vec[-(1:Ne)]
		} else ppars$dTau.e <- args$dTau.e		
	}
	if (args$par["Tau.i"] == T) {
			ppars$Tau.i <- reshapeVec2List.no.perm(args$Wc.i, pars.vec[1:Ni]); pars.vec <- pars.vec[-(1:Ni)]
	} else ppars$Tau.i <- args$Tau.i

	if (!is.null(args$X.ahp)){
		if (args$par["W.ahp"] == T) {
			ppars$W.ahp <- pars.vec[1:n.basis]; pars.vec<-pars.vec[-(1:n.basis)]
		} else ppars$W.ahp <- args$W.ahp
	} else {
		ppars['W.ahp'] <- list(NULL)
	}

	if (!is.null(args$delay.t)) ppars$delay.t <- args$delay.t

	if (args$linear.soma) ppars$Th[1] <- NA
	
	# print(paste('alphasyn = ', alphasyn))
		
	if (ret.parlist){
		err <- ppars
	} else {
		if (ret.resp) {
			err <- sim.hGLM(X=args$inp, dt=args$dt, pars=ppars, double=args$double, scale=args$scale, X.ahp=args$X.ahp, alphasyn=args$alphasyn) 
		} else {
			err <- sim.hGLM(X=args$inp, dt=args$dt, pars=ppars, vv=args$v, regpars=args$regpars, double=args$double, scale=args$scale, X.ahp=args$X.ahp, verbose=0, alphasyn=args$alphasyn)
			# X <- args$inp; dt <- args$dt; pars <- ppars; vv <- args$v; regpars <- args$regpars; double <- args$double; scale <- args$scale; X.ahp <- args$X.ahp; verbose <- 0
		}
	}
	# if (ret.resp) err <- list(resp=err, par=ppars)
	err
}

##########################################################
## a wrapper to calculate the gradient of the error from the parmeters and the arguments
grad.err.hGLM <- function(pars.vec, args){
## pars.vec: a parameter vector for the parameters of the LN model
	## args: $inp is the input spike train
	## 		$v is the target output
	## 		$dt is the time resolution of the input
	## 		$Jc is the structure of the dendrite
	## 		$Wc.e is the connectivity of the dendrite (excitatory synapses)
	## 		$Wc.i is the connectivity of the dendrite (inhibitory synapses)
	##			$linear.soma whether the soma is linear or nonlinear
	##	 The following vector contains the parameters to be optimised
	## 		$par=c(v0=T, Jw=T, Ww.e=T, Ww.e2=T, Ww.i=T, Th=T, Tau.e=T, Tau.i=T, W.ahp=T) - which parameters are in the parameter vector
	##					only parameters specified here will be optimised. All F has to be present amongst the arguments
	##		$v0 = NULL - if ($par["v0"] == F)  this will be used to evaluate the hGLM model
	## 	$Jw = NULL
	##		$Ww.e = NULL
	##		$Ww.e2 = NULL
	##		$Ww.i = NULL
	##		$Th = NULL
	##		$Tau.e = NULL
	##		$Tau.i = NULL
	##		$W.ahp = NULL

	# test.args(pars.vec, args)

	if (!args$alphasyn & args$double) stop('error: alphasyn is FALSE and double is TRUE - double exponential synapses are only implemented in single-kernel mode.')

	M <- length(args$Jc)
	Ne <- length(unlist(args$Wc.e))
	Ni <- length(unlist(args$Wc.i))
	if (!is.null(args$X.ahp)) 	n.basis <- nrow(args$X.ahp)

	ppars <- list(Jc=args$Jc, Wc.e=args$Wc.e, Wc.i=args$Wc.i)

	if (args$par["v0"] == T) {
			ppars$v0 <- pars.vec[1]; pars.vec <- pars.vec[-1]
		} else ppars$v0 <- args$v0

	if (args$par["Jw"] == T) {
			ppars$Jw <- pars.vec[1:M]; pars.vec <- pars.vec[-(1:M)]
		} else ppars$Jw <- args$Jw
	if (args$par["Ww.e"] == T) {
			ppars$Ww.e <- reshapeVec2List.no.perm(args$Wc.e, pars.vec[1:Ne]); pars.vec <- pars.vec[-(1:Ne)]
	} else ppars$Ww.e <- args$Ww.e
	if (args$double){
		if (args$par["Ww.e2"] == T) {
			ppars$Ww.e2 <- reshapeVec2List.no.perm(args$Wc.e, pars.vec[1:Ne]); pars.vec <- pars.vec[-(1:Ne)]
		} else ppars$Ww.e2 <- args$Ww.e2
	} else {
		ppars['Ww.e2'] <- list(NULL)		
	}
	if (args$par["Ww.i"] == T) {
			ppars$Ww.i <- reshapeVec2List.no.perm(args$Wc.i, pars.vec[1:Ni]); pars.vec <- pars.vec[-(1:Ni)]
	} else ppars$Ww.i <- args$Ww.i
		
	if (args$par["Th"] == T) {
			ppars$Th <- pars.vec[1:M]; pars.vec<-pars.vec[-(1:M)]
		} else ppars$Th <- args$Th
		
	if (args$par["Tau.e"] == T) {
			ppars$Tau.e <- reshapeVec2List.no.perm(args$Wc.e, pars.vec[1:Ne]); pars.vec <- pars.vec[-(1:Ne)]
	} else ppars$Tau.e <- args$Tau.e
	if (!args$alphasyn){
		if (args$par["dTau.e"] == T) {
				ppars$dTau.e <- reshapeVec2List.no.perm(args$Wc.e, pars.vec[1:Ne]); pars.vec <- pars.vec[-(1:Ne)]
				# print(ppars$dTau.e)
		} else ppars$dTau.e <- args$dTau.e		
	}
	if (args$par["Tau.i"] == T) {
			ppars$Tau.i <- reshapeVec2List.no.perm(args$Wc.i, pars.vec[1:Ni]); pars.vec <- pars.vec[-(1:Ni)]
	} else ppars$Tau.i <- args$Tau.i

	if (!is.null(args$X.ahp)) {
		if (args$par["W.ahp"] == T) {
				ppars$W.ahp <- pars.vec[1:n.basis]; pars.vec<-pars.vec[-(1:n.basis)]
		} else ppars$W.ahp <- args$W.ahp
	} else {
		ppars['W.ahp'] <- list(NULL)
	}

	if (!is.null(args$delay.t)) ppars$delay.t <- args$delay.t

	if (args$linear.soma) ppars$Th[1] <- NA
	
	gv <- grad.hGLM(X=args$inp, dt=args$dt, pars=ppars, vv=args$v, regpars=args$regpars, double=args$double, scale=args$scale, X.ahp=args$X.ahp, alphasyn=args$alphasyn)
	
	g.err <- 0
	if (args$par["v0"] == T) g.err <- c(g.err, gv$dv0)
	if (args$par["Jw"] == T) g.err <- c(g.err, gv$dj)
	if (args$par["Ww.e"] == T) g.err <- c(g.err, gv$dw.e)
	if (args$double){
		if (args$par["Ww.e2"] == T) g.err <- c(g.err, gv$dw.e2)	
	}
	if (args$par["Ww.i"] == T) g.err <- c(g.err, gv$dw.i)
	if (args$par["Th"] == T) g.err <- c(g.err, gv$dth)
	if (args$par["Tau.e"] == T) g.err <- c(g.err, gv$dtau.e)
	if (!args$alphasyn) {
		if (args$par["dTau.e"] == T) g.err <- c(g.err, gv$ddtau.e)
	}
	if (args$par["Tau.i"] == T) g.err <- c(g.err, gv$dtau.i)
	if (!is.null(args$X.ahp)) {
		if (args$par["W.ahp"] == T) g.err <- c(g.err, gv$dw.ahp)
	}

	g.err <- g.err[-1]
	g.err
}

##########################################################
## a function that tests that the pars.vec and the args are completely specified

test.args <- function(pars.vec, args){

	## First we check the arguments - parameters we are not optimising
	if (prod(c("inp", "v", "dt", "Jc", "Wc.e", "Wc.i","par", "linear.soma") %in% names(args)) == F) stop("some argument is missing")
	## check that the length of the parameter vector is as it is expected
	M <- length(args$Jc)
	Ne <- length(unlist(args$Wc.e))
	Ni <- length(unlist(args$Wc.i))
	if (!(is.null(args$X.ahp))) n.basis <- nrow(args$X.ahp)

	## check that the length of the parameter vector is as it is expected
	l.pars <- 0
	if (args$par["v0"] == T)  l.pars <- l.pars + 1
	if (args$par["Jw"] == T)  l.pars <- l.pars + M
	if (args$par["Ww.e"] == T)   l.pars <- l.pars + Ne
	if (args$double) if (args$par["Ww.e2"] == T)   l.pars <- l.pars + Ne
	if (args$par["Ww.i"] == T)   l.pars <- l.pars + Ni
	if (args$par["Th"] == T)   l.pars <- l.pars + M		
	if (args$par["Tau.e"] == T)   l.pars <- l.pars + Ne
	if (!args$alphasyn) if (args$par["dTau.e"] == T)   l.pars <- l.pars + Ne
	if (args$par["Tau.i"] == T)   l.pars <- l.pars + Ni
	if (!(is.null(args$X.ahp))) if (args$par["W.ahp"] == T)   l.pars <- l.pars + n.basis

	if (length(pars.vec) != l.pars) {
		cat('length of the parameters is incorrect. Found:', length(pars.vec), 'expected:', l.pars, '\n')
		stop('we will stop')
	}

	## checking that parameters are either given as arguments or in the parameter vector
	if (!args$par["v0"])	{
		if ("v0" %in% names(args)){
			if (length(args$v0) != 1) stop("length of v0 must be 1")
		} else stop("v0 not specified")
	}

	if (!args$par["Jw"])	{
		if ("Jw" %in% names(args)){
			if (length(args$Jw) != M) stop("length of Jw must be M")
		} else stop("Jw not specified")
	}

	if (!args$par["Ww.e"])	{
		if ("Ww.e" %in% names(args)){
			if (length(unlist(args$Ww.e)) != Ne) stop("length of Ww.e must be Ne")
		} else stop("Ww.e not specified")
	}

	if (args$double){
		if (!args$par["Ww.e2"])	{
			if ("Ww.e2" %in% names(args)){
				if (length(unlist(args$Ww.e2)) != Ne) stop("length of Ww.e2 must be Ne")
			} else stop("Ww.e2 not specified")
		}
	}
	
	if (!args$par["Ww.i"])	{
		if ("Ww.i" %in% names(args)){
			if (length(unlist(args$Ww.i)) != Ni) stop("length of Ww.i must be Ni")
		} else stop("Ww.i not specified")
	}

	if (!args$par["Th"])	{
		if ("Th" %in% names(args)){
			if (length(args$Th) != M) stop("length of Th must be M")
		} else stop("Th not specified")
	}

	if (!args$par["Tau.e"])	{
		if ("Tau.e" %in% names(args)){
			if (length(unlist(args$Tau.e)) != Ne) stop("length of Tau.e must be Ne")
		} else stop("Tau.e not specified")
	}

	if (!args$alphasyn){
		if (!args$par["dTau.e"])	{
			if ("dTau.e" %in% names(args)){
				if (length(unlist(args$dTau.e)) != Ne) stop("length of dTau.e must be Ne")
			} else stop("dTau.e not specified")
		}
	}

	if (!args$par["Tau.i"])	{
		if ("Tau.i" %in% names(args)){
			if (length(unlist(args$Tau.i)) != Ni) stop("length of Tau.i must be Ni")
		} else stop("Tau.i not specified")
	}
}

make.args.parvec <- function(pars, v, X, dt, X.ahp=NULL, optimvec=NULL, regpars=NULL){
	# this function prepares the inputs for the err and grad functions above
	# 
	# arguments:
	# pars: list of parameters, as used by the sim.hGLM() function
	# optimvec: name of the parameters to be optimised
	# other parameters passed directly to args:
	# v: vector of the target response, to calculate the error and gradients
	# X: input spike train
	# dt: simulation timestep
	# regpars: regularisation parameters
	#
	# returns:
	# args.parvec: a list with two components:
	#		- parvec: parameter vector
	# 		- args: list of arguments

	double <- F
	if ('Ww.e2' %in% names(pars)) {
		if (!is.null(pars$Ww.e2)) {
			double <- T
		}
	}
	alphasyn <- T
	if (double==FALSE & ('dTau.e' %in% names(pars)) & is.null(regpars)) {
		if (!is.null(pars$dTau.e)) {
			alphasyn <- F
		}
	}
	
	## first, we need to reorder the pars list
	## the ordering should be the following:
	## v0, Jw, Ww.e, Ww.e2, Ww.i, Th, Tau.e, dTau.e, Tau.i, W.ahp
	if (is.na(pars$Th[1])) linear.soma <- T else linear.soma <- F
	if (is.null(optimvec)){
		newpars <- list(v0=pars$v0, Jw=pars$Jw, Ww.e=pars$Ww.e)
		if (double) newpars$Ww.e2 <- pars$Ww.e2
		newpars$Ww.i <- pars$Ww.i
		newpars$Th <- pars$Th
		if (linear.soma) newpars$Th[1] <- 0
		newpars$Tau.e <- pars$Tau.e
		if (!alphasyn) newpars$dTau.e <- pars$dTau.e
		newpars$Tau.i <- pars$Tau.i
		if ('W.ahp' %in% names(pars)) newpars$W.ahp <- pars$W.ahp
	} else {
		newpars <- list()
		if ('v0' %in% optimvec) newpars$v0 <- pars$v0
		if ('Jw' %in% optimvec) newpars$Jw <- pars$Jw
		if ('Ww.e' %in% optimvec) newpars$Ww.e <- pars$Ww.e
		if (double & ('Ww.e2' %in% optimvec)) newpars$Ww.e2 <- pars$Ww.e2
		if ('Ww.i' %in% optimvec) newpars$Ww.i <- pars$Ww.i
		if ('Th' %in% optimvec) {
			newpars$Th <- pars$Th
			if (linear.soma) newpars$Th[1] <- 0
		}
		if ('Tau.e' %in% optimvec) newpars$Tau.e <- pars$Tau.e
		if (!alphasyn & ('dTau.e' %in% optimvec)) newpars$dTau.e <- pars$dTau.e
		if ('Tau.i' %in% optimvec) newpars$Tau.i <- pars$Tau.i
		if ('W.ahp' %in% optimvec) newpars$W.ahp <- pars$W.ahp
	}


	
	## then prepare the parameter vector
	pars.vec <- unlist(newpars, use.names=F)
	
	## a vector that tells the optimiser which parameter is present
	parnames <- c(v0=F, Jw=F, Ww.e=F, Ww.e2=F, Ww.i=F, Th=F, Tau.e=F, dTau.e=F, Tau.i=F, W.ahp=F)
	for (i.newpar in names(newpars)) {
		if (length(unlist(newpars[[i.newpar]])) > 0)	{
			parnames[i.newpar] <- T
		} else {
			## remove from newpars all missing parameters
			newpars[i.newpar] <- NULL
		}
	}
	
	## prepare the arguments list
	args <- list(v=v, par=parnames, linear.soma=linear.soma, inp=X, dt=dt, double=double, alphasyn=alphasyn, scale=2.8, X.ahp=X.ahp, regpars=regpars)

	## include all parameters not in the parameter vector in the argument list (these are not optimised for)

	if ('delay.t' %in% names(pars)) args$delay.t <- pars$delay.t
	for (nn in names(pars)){
		if ((!(nn %in% names(newpars))) & (!(nn %in% names(args)))) {
			args[[nn]] <- pars[[nn]]
		} 
	}
	
	Wc.i.vec <- unlist(pars$Wc.i)
	if (is.null(Wc.i.vec)){
		args$par['Tau.i'] <- F
		args$par['Ww.i'] <- F
		args$Tau.i <- pars$Tau.i
		args$Ww.i <- pars$Ww.i
	}
	
	args.parvec <- list(args=args, parvec=pars.vec)
	args.parvec
}


