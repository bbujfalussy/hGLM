
##############################################
## testing that the parameters are formally correct
test.pars <- function(v0, Jc, Jw, Wc.e, Wc.i, Ww.e, Ww.e2, Ww.i, Th, Tau.e, Tau.i, dTau.e, delay.t, N, verb=F){
	M <- length(Jc)
	if (Jc[1] != 0 ) stop("The first subunit must be the soma, that is connected to subunit 0 - meaning no connections")
	if (max(Jc > M)) stop("Can not connect to a non-existent subunit!")
	if (sum(Jc==0) > 1) stop("There must be only one subunit connected to subunit 0 - it is the soma!")
	if (max(Jc - seq(1,M)) > (-1)) stop("subunit - Jc - connectivity must be defined in increasing order!")
	
	##########################################
	## test that all parameters have the right dimensions
	if (length(v0)>1) stop("v0 must be a scalar")

	if (length(Jw) != M) stop("Jw must be of length M")

	if (!is.list(Wc.e)) stop("Wc.e must be a list")
	if (length(Wc.e) != M) stop("Wc.e must be a list of length M")
	Wc.e.vec <- unlist(Wc.e)
	if (!is.null(Wc.e.vec)){
		Ne <- max(unlist(Wc.e))
		if (Ne > N) stop("max of Wc.e must be smaller than N")
	} else Ne <- 0

	if (!is.list(Wc.i)) stop("Wc.i must be a list")
	if (length(Wc.i) != M) stop("Wc.i must be a list of length M")
	Wc.i.vec <- unlist(Wc.i)
	if (!is.null(Wc.i.vec)){	
		if (max(Wc.i.vec) > N) stop("max of Wc.i must be smaller than N")
		if (min(Wc.i.vec) <= Ne) stop("Excitatory neurons should come first")
	}

	if (length(unique(unlist(Wc.i))) != length(unlist(Wc.i))){
		duplicated.i <- unlist(Wc.i)[duplicated(unlist(Wc.i))]
		stop(paste("input(s)", duplicated.i, "are duplicated. Currently duplicated synapses are not implemented, so duplicate inputs instead...\n"))
	} 

	if (length(unique(unlist(Wc.e))) != length(unlist(Wc.e))){
		duplicated.i <- unlist(Wc.e)[duplicated(unlist(Wc.e))]
		stop(paste("input(s)", duplicated.i, "are duplicated. Currently duplicated synapses are not implemented, so duplicate inputs instead...\n"))
	} 

	if (sum(unique(unlist(Wc.i)) %in% unique(unlist(Wc.e)))) {
		Inh <- unique(unlist(Wc.i))
		Exc <- unique(unlist(Wc.e))
		InhExc <- which(Inh %in% Exc)
		message <- paste("Neuron", Inh[InhExc], "is both excitatory and inhibitory.")
		stop(message)
	}

	if (!is.list(Ww.e)) stop("Ww.e must be a list")
	if (length(Ww.e) != M) stop("Ww.e must be a list of length M")
	if (compare.lists(Wc.e, Ww.e)==F) stop("Ww.e and Wc.e must have the same format!")

	if (!is.null(Ww.e2)){
		if (!is.list(Ww.e2)) stop("Ww.e2 must be a list")
		if (length(Ww.e2) != M) stop("Ww.e2 must be a list of length M")
		if (compare.lists(Wc.e, Ww.e2)==F) stop("Ww.e2 and Wc.e must have the same format!")
	}

	if (!is.list(Ww.i)) stop("Ww.i must be a list")
	if (length(Ww.i) != M) stop("Ww.i must be a list of length M")
	if (compare.lists(Wc.i, Ww.i)==F) stop("Ww.i and Wc.e must have the same format!")
	
	if (length(Th) >1){
		if(length(Th) != M) stop("The length of Th must be the same as the number of branches")
	}

	if (length(delay.t) >1){
		if(length(delay.t) != M) stop("The length of delay.t must be the same as the number of branches")
	}
	
	if (!is.list(Tau.e)) stop("Tau.e must be a list")
	if (length(Tau.e) != M) stop("Tau.e must be a list of length M")
	if (compare.lists(Wc.e, Tau.e)==F) stop("Tau.e and Wc.e must have the same format!")

	if (!is.null(dTau.e)){
		if (!is.list(dTau.e)) stop("dTau.e must be a list")
		if (length(dTau.e) != M) stop("dTau.e must be a list of length M")
		if (compare.lists(Tau.e, dTau.e)==F) stop("dTau.e and Tau.e must have the same format!")
	}

	if (!is.list(Tau.i)) stop("Tau.i must be a list")
	if (length(Tau.i) != M) stop("Tau.i must be a list of length M")
	if (compare.lists(Wc.i, Tau.i)==F) stop("Tau.i and Wc.i must have the same format!")

	if (verb) print("Parameters formally tested...")
}


# formally checking the parameter lists
compare.lists <- function(V, W){
	# checks the forms of the two lists - they should be identical
	# V is the reference, if V[[k]] is NULL then no test is performed 
	L.V <- length(V)
	L.W <- length(W)
	out <- T
	if (L.V != L.W) {
		out=F
		if (min(L.V, L.W) > 0) warning("list lengths do not match!")	
	} else {
		for (i in 1:L.V){
			if (!is.null(V[[i]])){
				LL.V <- length(V[[i]])
				LL.W <- length(W[[i]])
				if (LL.V != LL.W) {
					out=F
					warning("sublist lengths do not match!")	
				}
			}
		}
	}
	out
}

# ###################################################################
# ## a function that initializes parameters based on a template list and the parameter range - used in init.hGLM
makeParList <- function(templateList, min, max){
	newList <- list()
	L <- length(templateList)
	for (i in 1:L){
		LL <- length(templateList[[i]])
		if (LL > 0)	newList[[i]] <- runif(LL, min, max) else newList[i] <- list(NULL)
	}
	newList
}



# in earlier version, when each subunit had a single synapse, the parameters were vectors and not lists
# this function converts from the vector-version to the list version
reshapeParList <- function(templateList, pars){
	## take the form of the templateList and transfer the parameters from pars
	if ('v' %in% names(templateList)) templateList$v <- NULL
	K <- length(templateList)
	for (k in 1:K){
		name.par <- names(templateList)[k]
		i.par <- which(names(pars) == name.par)
		if (length(i.par) == 0) {
			if (is.null(templateList[[k]])) pars[name.par] <- list(NULL) else pars[[name.par]] <- templateList[[k]]
		} else {
			list.in.template <- is.list(templateList[[k]]) 
			list.in.pars <- is.list(pars[[i.par]]) 
			if (list.in.template & (!list.in.pars)){
				newlist <- list()
				KK <- length(templateList[[k]])
				for (kk in 1:KK) {
					newpar <- pars[[i.par]][kk]
					if (name.par %in% c('Ww.e', 'Ww.e2')) newpar <- exp(newpar)
					if (name.par == 'Ww.i') newpar <- (-1) * exp(newpar)
					newlist[[kk]] <- newpar
					}
				pars[[name.par]] <- newlist
			}
		}
	}
	pars
}

## when synapse parameters are defined by list, it is not allowed to have parameter for non-existing synapse (where Wc.e is NULL)
check.nulls <- function(pars){
	## where Wc.e (Wc.i) is NULL, Tau.e, Ww.e and Ww.e2 (Ww.i, Tau.i) should also be NULL
	Wc.e <- pars$Wc.e
	Wc.i <- pars$Wc.i
	e.null <- unlist(lapply(Wc.e, is.null))
	pars$Ww.e[e.null] <- list(NULL)
	pars$Tau.e[e.null] <- list(NULL)
	if ('Ww.e2' %in% names(pars)) pars$Ww.e2[e.null] <- list(NULL)
	i.null <- unlist(lapply(Wc.i, is.null))
	pars$Ww.i[i.null] <- list(NULL)
	pars$Tau.i[i.null] <- list(NULL)
	pars
}


#############################################################
## a function that uncouples the parameters of a coupled model to allow different synapses
decouple.pars <- function(Jc, Wc.e, Wc.i, synpars){
	# Jc <- Jc.3Ns; Wc.e <- Wc.e.3Ns; Wc.i <- Wc.i.3Ns; synpars <- errs$pars.3LD$opt
	
	newpars <- list(Jc=Jc, Wc.e=Wc.e, Wc.i=Wc.i, v0=NULL, Jw=NULL, Ww.e=NULL, Ww.e2=NULL, Ww.i=NULL, Th=NULL, Tau.e=NULL, Tau.i=NULL, W.ahp=NULL, delay.t=NULL)
	newpars$v0 <- synpars$v0
	newpars$Jw <- synpars$Jw
	newpars$Th <- synpars$Th
	
	if ('delay.t' %in% names(synpars)) newpars$delay.t <- synpars$delay.t else newpars$delay.t <- NULL
	if ('W.ahp' %in% names(synpars)) newpars$W.ahp <- synpars$W.ahp else newpars$W.ahp <- NULL
	
	newpars$Ww.e <- decouple.single.par(Wc.e, synpars$Ww.e)
	newpars$Ww.i <- decouple.single.par(Wc.i, synpars$Ww.i)
	if ('Ww.e2' %in% names(synpars)) {
		if (!is.null(synpars$Ww.e2)){
			newpars$Ww.e2 <- decouple.single.par(Wc.e, synpars$Ww.e2)
		# newpars$Ww.e2 <- makeParListPositive(newpars$Ww.e2)
		}
	}

	# newpars$Ww.e <- makeParListPositive(newpars$Ww.e)
	# newpars$Ww.i <- makeParListPositive(newpars$Ww.i, negative=T)

	newpars$Tau.e <- decouple.single.par(Wc.e, synpars$Tau.e)
	if ('dTau.e' %in% names(synpars)) {
		if (!is.null(synpars$dTau.e)){
			newpars$dTau.e <- decouple.single.par(Wc.e, synpars$dTau.e)
		}
	}
	newpars$Tau.i <- decouple.single.par(Wc.i, synpars$Tau.i)
	newpars
}


decouple.single.par <- function(templateList, coupledList){
	newList <- list()
	L <- length(templateList)
	L2 <- length(coupledList)
	if (L!=L2) stop('template and coupled lists do not match!')
	for (i in 1:L){
		n <- length(templateList[[i]])
		m <- length(coupledList[[i]])
		if (m>1) {
			warning('coupled list contains more than 1 synapse per subunit, only the first will be used!') 
			coupledList[[i]] <- coupledList[[i]][1]
		}
		if (n > 0) newList[[i]] <- rep(coupledList[[i]], n) else newList[i] <- list(NULL)
	}
	newList
}

############################################################
## Wc.e, Wc.i: uncoupled connectivity parameters
## regp: mini parameters, as calculated by sim.hGLM(X, dt, pars, calc.minis=T)
## alpha.Ww, alpha.Tau: prior precisions
regpars.from.coupled <- function(Wc.e, Wc.i, minis, alpha.Ww, alpha.Tau){
	regpars <-  list(alpha.Ww=alpha.Ww, alpha.Tau=alpha.Tau)
	Ni <- length(unlist(minis$a.i))
	Ne <- length(unlist(minis$a.e))
	
	regpars$Ww.e1 <- unlist(decouple.regpar(Wc.e, minis$a.e))
	if ('a.e2' %in% names(minis)) regpars$Ww.e2 <- unlist(decouple.regpar(Wc.e, minis$a.e2))
			
	regpars$Ww.e1 <- makeParPositive(regpars$Ww.e1)
	regpars$Ww.e2 <- makeParPositive(regpars$Ww.e2)

	if (Ni > 0) {
		regpars$Ww.i <- unlist(decouple.regpar(Wc.i, minis$a.i))
		regpars$Ww.i <- makeParPositive(regpars$Ww.i, negative=T)
		regpars$logTau.i <- unlist(decouple.regpar(Wc.i, minis$logTau.i))
	} else {
		regpars$Ww.i <- list(NULL)
		regpars$logTau.i <- list(NULL)
	}

	regpars$logTau.e <- unlist(decouple.regpar(Wc.e, minis$logTau.e))
	regpars
}

decouple.regpar <- function(templateList, mini){
	newpar <- 0
	template.not.null <- which(!unlist(lapply(templateList, is.null)))
	templateList <- templateList[template.not.null]
	L <- length(templateList)
	L2 <- length(mini)
	if (L!=L2) stop('template and coupled lists do not match!')
	
	for (i in 1:L){
		n <- length(templateList[[i]])
		if (n > 0) newpar <- c(newpar, rep(mini[i], n))
	}
	newpar <- newpar[-1]
	newpar
}

makeParPositive <- function(par, negative=F){
	if (negative) {
		par.pos <- which(par > 0)
		par[par.pos] <- NA
		par[par.pos] <- max(par, na.rm=T)
	} else {
		par.neg <- which(par < 0)
		par[par.neg] <- NA
		par[par.neg] <- min(par, na.rm=T)
	}
	par
}

makeParListPositive <- function(par, negative=F){
	if (negative) {
		ppar <- unlist(par)
		ppar[ppar > 0] <- NA
		if (sum(is.na(ppar)) < length(ppar)) maxpar <- max(ppar, na.rm=T) else maxpar <- 0
		for ( i in 1:length(par)){
			if (length(par[[i]] > 0)){
				par[[i]][which(par[[i]] > 0)] <- maxpar
			}			
		}
	} else {
		ppar <- unlist(par)
		ppar[ppar < 0] <- NA
		if (sum(is.na(ppar)) < length(ppar)) minpar <- min(ppar, na.rm=T) else minpar <- 0
		for ( i in 1:length(par)){
			if (length(par[[i]] > 0)){
				par[[i]][which(par[[i]] < 0)] <- minpar
			}			
		}
	}
	par
}


###################################################################
## functions to prepare a parameter list form a vector version
## similar to reshapeParList - but these functions only handle a single parameter type, but for a complex hierarchy
## the first one applies a permutation in the vector - it uses the ordering from the template

reshapeVec2List <- function(templateList, synpars, Ne=0, sign=0){
	## from a parameter vector - synpars
	## make a list of parameters based on the connectivity data - templateList (Wc.e Wc.i)
	## used in reshape.pars.syn2()
	
	newList <- list()
	L <- length(templateList)
	for (i in 1:L){
		sublist <- templateList[[i]] - Ne
		if (length(sublist) > 0)	{
			newList[[i]] <- synpars[sublist] # here we select parameters according to the connectivity
			if (sign != 0) newList[[i]] <- sign * abs(newList[[i]])
		} else newList[i] <- list(NULL)
	}
	newList
}

# same as above, but without permutation - used in err.hGLM
reshapeVec2List.no.perm <- function(templateList, synpars, sign=0){
# same as above, but without permutation
	newList <- list()
	L <- length(templateList)
	K0 <- 1
	for (i in 1:L){
		K <- length(templateList[[i]])
		if (K > 0)	{
			newList[[i]] <- synpars[K0:(K0+K-1)] # the original order is restored
			if (sign != 0) newList[[i]] <- sign * abs(newList[[i]])
			K0 <- K0 + K
		} else newList[i] <- list(NULL)
	}
	if (K0 < length(synpars)) stop("Vec2List: not all parameters are used!")
	newList
}



#############################################################
## a function that reshapes the 1LsynD parameters into a hierarchical architecture
# Jc <- Jc.1N; Wc.e <- Wc.e.2N; Wc.i <- Wc.i.2N; synpars <- errs$pars.1Lsyn3DE$opt
# this is a formal transformation required before the init.hgLM function can match the hierarchy with the linear

reshape.pars.syn2 <- function(Jc, Wc.e, Wc.i, synpars){
	newpars <- list(Jc=Jc, Wc.e=Wc.e, Wc.i=Wc.i, v0=NULL, Jw=NULL, Ww.e=NULL, Ww.e2=NULL, Ww.i=NULL, Th=NULL, Tau.e=NULL, dTau.e=NULL, Tau.i=NULL)
	i.inh <- unique(unlist(Wc.i)); i.exc <- unique(unlist(Wc.e))
	N.i <- length(i.inh); N.e <- length(i.exc)
	M <- length(Wc.e) # number of subunits
	newpars$v0 <- synpars$v0
	newpars$Jw <- rep(log(4), M); if (!is.na(synpars$Jw)) newpars$Jw[1] <- synpars$Jw[1]
	
	newpars$Ww.e <- reshapeVec2List(Wc.e, synpars$Ww.e[[1]], sign=0)
	newpars$Ww.i <- reshapeVec2List(Wc.i, synpars$Ww.i[[1]], sign=0, Ne=N.e)
	if ('Ww.e2' %in% names(synpars)) {
		if (!is.null(synpars$Ww.e2)) {
			newpars$Ww.e2 <- reshapeVec2List(Wc.e, synpars$Ww.e2[[1]], sign=0)
		} else newpars['Ww.e2'] <- list(NULL)
	}
	newpars$Th <- rep(0, M); if (!is.na(synpars$Th)) newpars$Th[1] <- synpars$Th
	newpars$Tau.e <- reshapeVec2List(Wc.e, synpars$Tau.e[[1]])
	if ('dTau.e' %in% names(synpars)) {
		if (!is.null(synpars$dTau.e)) {
			newpars$dTau.e <- reshapeVec2List(Wc.e, synpars$dTau.e[[1]], sign=0)
		} else newpars['dTau.e'] <- list(NULL)
	}
	newpars$Tau.i <- reshapeVec2List(Wc.i, synpars$Tau.i[[1]], Ne=N.e)
	newpars
}

#############################################################
## functions to prepare the inputs for the hGLM
#############################################################

inputs.hGLM <- function(X, we, wi, branches.subunits, syn.branch=Inf, pool.subunits=T, xename=NULL){
	## preparing the inputs for a hGLM with a given subunit structure
	##
	## args:
	## X: original input spike train
	## w: mapping from the synapses to dendritic branches (2 x N matrix, branch, position within branch)
	## branches.subunits: mapping from branches to subunits
	## syn.branch: leaves only n synapses per branch
	## pool.subunits: TRUE - pools all synapses per subunit;
	## xename: file that contains information about the location of the synapses within branches
	##
	## returns: a list containing the following elements:
	## Xnew: the new input spike train
	## Wc.e, Wc.i: connectivity between synapses and subunits
	if (nrow(X) != (nrow(we) + nrow(wi))) stop('input spikes and synapses do not match!')
	
	## 1. pool synapses within braches
	if (syn.branch < Inf){
		if (!is.null(xename)){
			if (file.exists(xename)){
				load(xename)
				X <- st$X
				we <- st$we
				wi <- st$wi
			} else {
				st <- reshape.stim(X, we, wi, nseg=syn.branch)
				save(st, file=xename)
				X <- st$X
				we <- st$we
				wi <- st$wi
			}
		} else {
			st <- reshape.stim(X, we, wi, nseg=syn.branch)
			X <- st$X
			we <- st$we
			wi <- st$wi
		}
	}
	
	## 2. leave X as it is, and calculate Wc from we and wi and branches.subunits
	Wc.e <- list()
	for (i in 1:length(branches.subunits)){
		Wc.e[[i]] <-  which(we[,1] %in% branches.subunits[[i]])
		if (length(Wc.e[[i]]) == 0) Wc.e[[i]] <- NULL
	}
	
	Wc.i <- list()
	for (i in 1:length(branches.subunits)){
		Wc.i[[i]] <-  which(wi[,1] %in% branches.subunits[[i]]) + nrow(we)
		if (length(Wc.i[[i]]) == 0) Wc.i[[i]] <- NULL
	}
	
	# 3. pool the spikes that belong to the same subunit
	# when pool, make sure, that both the input (X) and the connectivity (Wc) are modified!
	if (pool.subunits) newXW <- pool.spikes(X, Wc.e, Wc.i) else newXW <- list(X=X, Wc.e=Wc.e, Wc.i=Wc.i)
	newXW
}

####################################
## pool spikes - one spike train per subunit
####################################

pool.spikes <- function(X, Wc.e, Wc.i){
	## pool spikes that belong to the same subunit
	
	if (length(Wc.e) != length(Wc.i)) stop("Wc.e and Wc.i has to be of the same length!")
	nsyn.e <- sum(!unlist(lapply(Wc.e, is.null))) ## we don't include branches with no innervation
	nsyn.i <- sum(!unlist(lapply(Wc.i, is.null)))
	L <- ncol(X)
	
	X.new <- matrix(0, nsyn.e + nsyn.i, L)
	k <- 1
	for (i in 1:length(Wc.e)){
		if (!is.null(Wc.e[[i]])){
			xx <- matrix(X[Wc.e[[i]],], ncol=L)
			X.new[k,] <- colSums(xx)
			Wc.e[[i]] <- k
			k <- k + 1
		}
	}
	for (i in 1:length(Wc.i)){
		if (!is.null(Wc.i[[i]])){
			xx <- matrix(X[Wc.i[[i]],], ncol=L)
			X.new[k,] <- colSums(xx)
			Wc.i[[i]] <- k
			k <- k + 1
		}
	}
	list(X=X.new, Wc.e=Wc.e, Wc.i=Wc.i)
}

####################################
## pool spikes within a branch
####################################

reshape.stim <- function(X, wee, wii, nseg=3){
	# pool inputs (X) and branch-maps (wee and wii) and leaves max nseg inputs per branches
	if (nrow(X) != (nrow(wee) + nrow(wii))) stop("stimulus and synapse number does not match")

	# separate excitatory and inhibitory input
	Xe <- X[1:nrow(wee),]
	Xi <- X[(nrow(wee)+1):nrow(X),]

	# sort branches to increasing order
	Xe <- Xe[sort(wee[,1], index.ret=T)$ix,]
	Xi <- Xi[sort(wii[,1], index.ret=T)$ix,]
	wee <- wee[sort(wee[,1], index.ret=T)$ix,]
	wii <- wii[sort(wii[,1], index.ret=T)$ix,]
		
	dends.e <- sort(unique(wee[,1])) # branches receiving E input
	dends.i <- sort(unique(wii[,1]))
	n.br.e <- length(dends.e) # number of branches receiving E input
	n.br.i <- length(dends.i) # ... I input
	
	i.x <- 0

	for (i in 1:length(dends.e)){
		rWX <- reshape.X(Xe, nseg, dends.e[i], wee)
		wee <- rWX$w
		Xe <- rWX$X
		i.x <- c(i.x, max(i.x) + rWX$i.x)
	}
	
	for (i in 1:length(dends.i)){
		rWX <- reshape.X(Xi, nseg, dends.i[i], wii)
		wii <- rWX$w
		Xi <- rWX$X
		i.x <- c(i.x, max(i.x) + rWX$i.x)
	}

	i.x <- i.x[-1]

	stim <- rbind(Xe, Xi)
	st <- list(X=stim, wee=wee, wii=wii, i.x=i.x)
	st
}

reshape.X <- function(X, nseg, d, w){
	L <- ncol(X)
	
	Xd <- matrix(X[w[,1]==d,], ncol=L)
	w.d <- matrix(w[w[,1]==d,], ncol=2)
	X.nd <- matrix(X[w[,1]!=d,], ncol=L)
	w.nd <- matrix(w[w[,1]!=d,], ncol=2)

	w.d <- w.d[,2]
	# ww.d <- floor(w.d*(nseg+1)) # - why nseg + 1? It just doesn't make sense!
	ww.d <- floor(w.d*nseg) # 
	
	groups <- sort(unique(ww.d))
	# pos <- seq(1/nseg, 1, by=1/nseg) - 1/(2*nseg)
	pos <- seq(0, 1, length=nseg)
	pos <- pos[groups+1]
	
	i.x <- NA	
	XX.d <- rep(NA, ncol(X))

	for (i in 1:length(groups)){
		XX.d <- rbind(XX.d, colSums(matrix(Xd[ww.d==groups[i],], ncol=ncol(Xd), byrow=F) ) )
		i.x <- c(i.x, rep(i, sum(ww.d==groups[i])))
	}
	i.x <- i.x[-1]
	
	XX.d <- XX.d[-1,]
	X <- rbind(X.nd, XX.d)
	w <- rbind(w.nd, cbind(rep(d, length(groups)), pos))
	list(X=X, w=w, i.x=i.x)
}

