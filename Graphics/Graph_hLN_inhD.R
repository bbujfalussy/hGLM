
# # source('../InhExcDouble/Sim_hLN_InhD.R', chdir = TRUE)
# library(igraph)

# graph.hLN.inhD <- function(X, dt, pars, vv=NULL, logpars=T, ret.resp=F, graphics=T, term1, term2){
# ## X: NxT binary input matrix of presynaptic spikes; N: # of input neurons; T: # of timesteps
# ## dt: the time resolution of X in miliseconds
# ## v0: baseline of somatic voltage
# ## Jc: M length Connectivity vector for dendritic branches; M: # of dendritic branches
# ## 	for each branch is mother is given. 1 is the root which has a value 0. 
# ## 	A 3 layer deep binary tree has 7 subunits is given by: [0, 1, 1, 2, 2, 3, 3]
# ## Jw: M length Weight vector for the coupling weights associated with the dendritic branches
# ## Wc: a list of M components. The component m is a vector indicating the neurons connected to branch m
# ## Ww: M length vector indicating the synaptic weights associated to the neurons connected to branch m
# ## Th: M length vector of the thresholds for the sigmoidal nonlinearities in the branches
# ## Tau: M length vector of the synaptic time constants in miliseconds
# ## dTau: Same as Tau, strictly positive (DExp synapses) or NULL (if alpha synapses are used)
# ## Delay.t: M length vector of the propagation times in miliseconds

	# N <- nrow(X)
	# L <- ncol(X)
	# Tmax <- L / dt
	# M <- length(pars$Jc)
	# if (!is.null(vv)) {
		# if (length(vv)!=L) stop("The training signal and the input must have the same length!")
	# }

	# if (!("dTau.e" %in% names(pars))) pars$dTau.e <- rep(NA, M)
	# if (!("dTau.i" %in% names(pars))) pars$dTau.i <- rep(NA, M)
	# if (!("delay.t" %in% names(pars))) {
		# if (logpars) pars$delay.t <- log(rgamma(M, 3, scale=1/3)) else pars$delay.t <- rgamma(M, 3, scale=1/3)
	# }

	# #######################################################################
	# resp <- with(pars, {		
		# # v0 <- pars$v0; Jc <- pars$Jc; Jw <- pars$Jw; Wc.e <- pars$Wc.e; Wc.i <- pars$Wc.i; Ww.e <- pars$Ww.e; Ww.e2 <- pars$Ww.e2; Ww.i <- pars$Ww.i; Th <- pars$Th; Tau.e <- pars$Tau.e; Tau.i <- pars$Tau.i; delay.t <- pars$delay.t

		# test.pars.inhD(v0, Jc, Jw, Wc.e, Wc.i, Ww.e, Ww.e2, Ww.i, Th, Tau.e, Tau.i, delay.t, N)	
		
		# ## we will mostly constrain Jw and Tau.e to be positive. The simplest way to achieve this is to define them by their log
		# ## in this case we have to transform them back to the noral form before processing
		# if (logpars) {
			# Jw <- exp(Jw)
			# Tau.e <- exp(Tau.e)	
			# Tau.i <- exp(Tau.i)	
			# Ww.e <- exp(Ww.e)	
			# Ww.e2 <- exp(Ww.e2)	
			# Ww.i <- (-1) * exp(Ww.i)	
			# delay.t  <- exp(delay.t)
		# }
		# Jw[is.na(Th)] <- 1 ## for linear subunits Jw is redundant with Ww.e - so we set it to 1
		
		# ##########################################
		# ## first, we calculate the synaptic input to each dendritic branch
		
		# Y <- matrix(0, M, L) # a matrix for the synaptic inputs to the subunits
		# for (m in 1:M){
			# if ((length(Wc.e[[m]]) + length(Wc.i[[m]])) > 0){
				# Y[m,] <- int.spikes(X, dt, Wc.e[[m]], Ww.e[m], Tau.e[m], NA, delay.t[m]) + int.spikes(X, dt, Wc.i[[m]], Ww.i[m], Tau.i[m], NA, delay.t[m]) + int.spikes(X, dt, Wc.e[[m]], Ww.e2[m], 10.4+2.8*Tau.e[m], NA, delay.t[m])
			# }			
		# }

		# ##########################################
		# ## next, we start from the leaves and apply the nonlinearities as well as add the inputs from the children
		# ## 1. Where are the leaves?
		# Jc.orig <- Jc

		# R <- matrix(0, M, L) # a matrix with the activation of the subunits
		# subunits.done <- 0 # leaves already processed
		# i.soma <- NULL # the index of the root subunit which is the soma
		# while (length(subunits.done) < (M+1)){ 
		# ## we repeat this until all the subunits ar eprcessed ...
			# i.leaves <- !seq(1,M) %in% Jc ## leaves are those subunits that do not receive input from other subunits
			# i.remain <- !seq(1,M) %in% subunits.done ## we only need the remaining leaves!
			# ii.leaves <- seq(1,M)[i.leaves & i.remain]
			# if (length(ii.leaves)==0) {
				# cat("No further leaves found, the graph may contain a loop among the following subunits: ", seq(1,M)[i.remain], " ")
				# stop(" We will quit.")
			# }
			# for (ii in ii.leaves){
			# # 2. apply the sigmoidal nonlinearity for every leaves
				# if (is.na(Th[ii])) {
					# if (ii !=1) warning(paste("Subunits are by definition nonlinear. All linear subunits should be merged into their parents. Subunit", ii, "seems to be linear."))
					# R[ii,] <- Y[ii,] 
				# } else R[ii,] <- sigm(Y[ii,], c=Th[ii])
	
			# # 3. ad the input from the child to the parent
				# ii.parent <- Jc[ii]
				# if (ii.parent>0){
					# Y[ii.parent,] <- Y[ii.parent,] + Jw[ii] * R[ii,]
				# } else {
					# if (ii != 1) warning(paste("The soma should be the first subunit. Now we found that subunit ", ii, "has no parent!"))
					# if (!is.null(i.soma)) stop(paste("This neuron has more than one some! Subunits", ii, "and", i.soma, "are both root subunits!"))
					# i.soma <- ii
				# }
				# Jc[ii] <- NA ## The subunit is already processed - does not give further input
				# subunits.done <- c(subunits.done, ii) ## register the subunits already processed				
			# }
		# }
		# v <- Jw[i.soma] * R[i.soma,] + v0

		# Jc <- Jc.orig
		# ncols <- length(term1)
		# nrows <- 2
		
		# names.dends <- paste("dend", 1:length(Jc)); names.dends[1] <- "soma"
		# ## plot the structure first
		# if ((length(Jc) > 1) & (graphics==T)){
			# # quartz("graph")
			# # el <- cbind(seq(2,length(Jc)), Jc[-1])
			# # gdend <- graph.edgelist(el)
			# # if (n.plots > 6) marg <- -0.5 else marg <- 0 
			# # plot(gdend, edge.arrow.size=1/3, margin=marg, vertex.label.cex=1, vertex.size=0.1)
			# # plot(gdend, edge.arrow.size=1/3, margin=marg, vertex.label.cex=1, vertex.size=0.1)
			# # quartz("nonlin", 12, 4)
			# par(mfrow=c(nrows, ncols)); par(mar=c(3,2,3,1))

		# }
		
		# all.inh <- rep(NA, length(Jc)) # we store the total inhibition to a subunit here
		# all.exc <- rep(NA, length(Jc)) # we store the total inhibition to a subunit here
		# E1perE2 <- rep(NA, length(Jc)) # we store the total inhibition to a subunit here
		# for (i in c(term1, term2)){
			# i.den <- i
			# ## the histogram of the input
			# if (diff(range(Y[i.den,])) > 0){
				# # 1. histogram of the inputs
				# hY <- hist(Y[i.den,], plot=F)
					# if (graphics){
						# plot(hY, ylab="", xlab="", col=grey(0.9), axes=F, main=""); axis(1, padj=-.8); mtext(names.dends[i.den], line=1.5, cex=0.7, font=2)
						# # cat("mean=", mean(Y[i.den,]), ", sd=", sd(Y[i.den,]), ", threshold=", Th[i.den], "output variance", var(Jw[i.den] * R[i.den,]), "\n")
					# }

				# # 2. plot of the sigmoid
				# max.hY <- max(hY$counts)
				# xx <- seq(min(Y[i.den,]), max(Y[i.den,]), length=100)
				# if (is.na(Th[i.den])) lines(c(xx[1], xx[100]), c(0, max.hY), col=3, lwd=Jws[j]) else {
					# yy <- sigm(xx, c=Th[i.den]); yy <- yy * max.hY
					# if (graphics) lines(xx, yy, col=3, lwd=Jws[j])
				# }
				# if ((length(Wc.e[[i.den]]) > 0) + (length(Wc.i[[i.den]]) > 0) ) {# we plot the kernel if there is synaptic inputs to the subunit
					# t <- seq(0, 100, by=0.5)
					# s <- matrix(0, 1, 201); s[1, 20] <- 1	
					# s2 <- matrix(0, 1, 601); s2[1, 20] <- 1	
					# if (length(Wc.e[[i.den]]) > 0) {
						# re <- as.vector(Ww.e[i.den] * int.spikes(s, dt=0.5, 1, 1, Tau.e[i.den], delay.t=delay.t[i.den]))
						# re2 <- as.vector(Ww.e2[i.den] * int.spikes(s2, dt=0.5, 1, 1, 10.4+2.8*Tau.e[i.den], delay.t=delay.t[i.den]))
					# } else re <- rep(0, length(t))
					# if (length(Wc.i[[i.den]]) > 0) {
						# ri <- as.vector(Ww.i[i.den] * int.spikes(s, dt=0.5, 1, 1, Tau.i[i.den], delay.t=delay.t[i.den]))
						# all.exc[i.den] <- sum(re)+sum(re2)
						# all.inh[i.den] <- sum(ri)/all.exc[i.den]
						# E1perE2[i.den] <- sum(re)/sum(re2)
					# } else ri <- rep(0, length(t))
					# if (graphics){
						# par(mfg=par()$mfg)
						# matplot(t, cbind(re, re2[1:201], ri), t="l", axes=F, xlab="", ylab="", main="", col=c(2,2,4), lwd=1, lty=1)
						# axis(3, col=2, line=0, padj=0.8)#, at=c(xx[1], labels=xx[100]), c(0, 200))
					# }
				# } else {
					# if (graphics){
						# plot(c(1,2,3,4,5), c(1,2,3,4,5), t="n", axes=F, main=names.dends[i.den]); axis(1)
					# }
				# }
			# }
		# }
		
		# # if (length(Jc) > 1){
			# # el <- cbind(seq(2,length(Jc)), Jc[-1])
			# # gdend <- graph.edgelist(el)
			# # if (n.plots > 6) marg <- -0.5 else marg <- 0 
			# # plot(gdend, edge.arrow.size=0.3, margin=marg, vertex.label.cex=1, vertex.size=n.plots*6)
		# # }
		
		# R <- normal.input(Y, Th, Tau.e)
		# R$pars <- cbind(R$pars, all.exc, all.inh, E1perE2)
		# R
	# })
	# if (ret.resp==T) resp else NULL
# }

# # #############################################################
# # # Testing this function
# # X <- matrix(rbinom(30*2000, 1, 0.005), 30, 2000)
# # dt <- 1
# # Jc <- c(0,1,2,2,3,3,5,5,6,6)
# # Wc.e <- list(NULL, c(1,2), c(3), c(4,5,6), c(7,8), NULL, c(9,10,11), c(12,13,14), c(15,16,17), c(18,19,20))
# # Wc.i <- list(NULL, 21, NULL, 22, 23, NULL, c(24, 25, 26), c(22, 27), c(23, 28, 29), c(21, 30))
# # initpars <- list(Jc=Jc, Wc.e=Wc.e, Wc.i=Wc.i)

# # pars <- init.hLN.inh(X, dt=1, initpars, lin.soma=T)

# # raster(X, dt=1)
# # v <- sim.hLN.inh(X, dt, pars)
# # lines(v+5, t="l", col=4)


# # graph.hLN.inh(X, dt, pars)
# # raster(X, dt=1)
# # lines(v+5, t="l", col=4)

# normal.input <- function(Y, th, tau){
	# n.subunits <- nrow(Y)
	# L <- ncol(Y)
	# Y.normal <- matrix(NA, n.subunits, L)
	# pars.sigm <- matrix(NA, n.subunits, 3, dimnames=list(NULL, c("slope", "threshold", "Tau.e")))#, "Tau.i", "W.e", "W.i", "Jj")))
	# for (subunit in 1:n.subunits){
		# yy <- Y[subunit,]
		# m.yy <- mean(yy)
		# s.yy <- sd(yy)
		# Y.normal[subunit,] <- (yy - m.yy) / s.yy

		# slope.norm <- s.yy 
		# th.norm <- (th[subunit] - m.yy) / s.yy 
		# pars.sigm[subunit,] <- c(slope.norm, th.norm, tau[subunit])
		
	# }
	# list(Y=Y.normal, pars=pars.sigm)
# }

