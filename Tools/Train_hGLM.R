## A script to train hLN model
# setwd('Tools')
source('./Err_hGLM.R')
source('./Sim_hGLM_sp.R')
source('../Utils/Extract_spiketimes.R')
source('../Utils/Spike_Basis.R')
source('../Utils/Estimate_threshold.R', chdir = TRUE)
source('../Utils/Spike_Reliability.R', chdir = TRUE)

# train <- NULL; test <- NULL; Tmax <- 5000; dt <- 1; seed <- 317; graphics.on <- 0; linear.soma <- F; double <- F; init <- T; preinit <- F; sim.spikes <- F; regpars <- NULL; maxit <- 100
# graphics.on <- 2

train.hGLM <- function(pars, train=NULL, test=NULL, dt=1, Tmax=5000, seed=317, graphics.on=0, linear.soma=F, double=F, alphasyn=T, init=T, preinit=F, sim.spikes=F, regpars=NULL, maxit=100){
	## a function to train hGLM models 
	## 
	## args:
	## pars: list of parameters defining the model
	## 		must contain  Jc, Wc.e, Wc.i
	## 		may contain parameters:
	##				Tau.e, Tau.i, delay.t: if omitted, these parameters will be initialised randomly
	##				Jw, Ww.e, Ww.i, Th, v0: these parameters are initialised randomly, but to ensure that all ' have aproximately equal effect on the output.
	## train: list with X and v, training input and training data
	## test: list with X and v, testing input and testing data
	## Tmax: if no data is provided, random data will be generated with length
	## 
	## returns:
	## pars: the list of the optimal parameter values
	## s.exp: the training and the test signal explained initially and after optimization
	## v: the predicted response
	## pars <- ipv1M; train=train.1N; test=test.1N; dt=1; Tmax=5000; seed=317; graphics.on=1; linear.soma=T; double=F; init=F; preinit=preinit; sim.spikes=F; regpars=NULL; maxit=100
	
	if (!alphasyn & double) stop('error: alphasyn is FALSE and double is TRUE - double exponential synapses are only implemented in single-kernel mode.')

	## 1. are there spikes in the training / test data?
	set.seed(seed)
	if (!is.null(train)){
		print('extracting spiketimes')
		spt <- extract.spiketimes(train$v, sampling.freq=1000/dt, limite=0, graphics=graphics.on)
		if (is.null(spt$st)) {
				sim.spikes <- F 
				cat('no spike have been found - spiking model will NOT be fitted \n')
			} else {
				sim.spikes <- T
				cat(length(spt$st), 'spikes have been found \n')
			}
	}	

	N <- max(c(unlist(pars$Wc.e), unlist(pars$Wc.i))) # number of synapses
	M <- length(pars$Wc.e)		# number of branches
	n.basis <- 8 # number of basis functions used to capture spike shape

	generate.train <- F
	 if (is.null(train)){
		L <- Tmax / dt; L.test <- L
	 	## generate training data	
		## 1) input - training data
		X <- matrix(rbinom(N*L, 1, 0.005*dt), N, L) # 5Hz training data
		X.test <- matrix(rbinom(N*L, 1, 0.005*dt), N, L) # 5Hz test data
		pars.sub <- init.hGLM(X, dt, pars, lin.soma=linear.soma, double=double, alphasyn=alphasyn, rand.seed=seed, rescale.soma=T, mean.v=-60, sd.v=5, nSD=1)

		if (sim.spikes){
			if (linear.soma) y.soma <- (pars.sub$v - pars.sub$v0)/ exp(pars.sub$Jw[1]) else y.soma <- isigm((pars.sub$v - pars.sub$v0)/exp(pars.sub$Jw[1]), c=pars.sub$Th[1])
			pars.spike <- default.spikes(dt=dt, n.basis=n.basis, Tmax=200, v=pars.sub$v, pars.sub=pars.sub, graphics=graphics.on, rate=5)
			pars.sub$W.ahp <- pars.spike$W.ahp
			resp <- sim.hGLM.sp(X, dt, pars.sub, pars.spike, rseed=seed+1, double=double, Nrep=1)
			train <- list(v=resp$v[1,] + resp$sp[1,], X=X)

			resp2 <- sim.hGLM.sp(X.test, dt, pars.sub, pars.spike, rseed=seed+2, double=double, Nrep=1)
			test <- list(v=resp2$v[1,] + resp2$sp[1,], X=X.test)
		} else {
			train <- list(v=pars.sub$v, X=X)

			resp2 <- sim.hGLM(X.test, dt, pars.sub, double=double)
			test <- list(v=resp2, X=X.test)			
		}
		
		if (graphics.on>1){
			# print('plotting the raster')
			raster(X)
			plot(train$v, t='l', main='train and test data')
			lines(test$v, t='l', col=2)
		}
		generate.train <- T
		print('Training and test signals are generated.')
		print(paste('The model has ', M, ' subunits and we generated ', L, ' training and ', L.test, ' test points.', sep=''))
	} else {
		L <- length(train$v);	L.test <- length(test$v)
		print('Training and test signals are received.')
		print(paste('The model has ', M, ' subunits and we will have ', L, ' training and ', L.test, ' test points.', sep=''))
	}
	
	################################################################
	## Prepare data for training
	################################################################
	## we need to set the v at the time of the spikes NA 
	## and rescale the subthreshold response to 0 mean and 1 variance
	## v - original
	## vsub - no APS

	v.orig <- train$v
	
	if (sim.spikes) {
		print('preparing subthreshold Vm for optimization')
		spt <- extract.spiketimes(train$v, sampling.freq=1000/dt, limite=-30, graphics=graphics.on)

		vsub <- as.vector(remove.spikes(train$v, spt, dt))
		LL <- sum(!is.na(vsub))

		## we need to convolve the spikes with the basis functions to get X.ahp
		basis.a <- gen.basis(Tmax=200, dt=dt, graphics=graphics.on, n.basis=n.basis, dur=spt$dur)	
		X.ahp <- spike.basis(spt, L, basis.a, dt=dt)
	} else {
		vsub <- train$v
		X.ahp <- NULL
		LL <- length(train$v)
	}
	
	m.v <- mean(vsub, na.rm=T)
	sd.v <- sd(vsub, na.rm=T)
	# v.scaled <- train$v - m.v; v.scaled <- v.scaled / sd.v
	# vsub.scaled <- vsub - m.v; vsub.scaled <- vsub.scaled / sd.v

	if (graphics.on > 1){
		plot(train$v, t='l')
		lines(vsub, col=2)

		# plot(v.scaled, t='l')
		# lines(vsub.scaled, col=2)
	}
	
	################################################################
	## initial parameters for learning
	##################################################################
	if (init){
		## random initialisation: 0 mean, unit variance, sparse AHP
		print('initialising all parameters... ')
		pars.init <- init.hGLM(train$X, dt, pars, X.ahp=X.ahp, lin.soma=linear.soma, rescale.soma=T, mean.v=m.v, sd.v=sd.v, double=double, alphasyn=alphasyn, rand.seed=seed+3, nSD=1)
		# if (graphics.on) {
			# matplot(cbind(pars.init$v, vsub), t='l', lty=1, xlab="time (bins)", ylab="response")
			# legend('topleft', bty="n", legend=c('init', 'reference'), lty=1, col=c(1,2))
		# }
		
		# if (sim.spikes){
			# pars$W.ahp <- init.ahp(train$v, dt, spt, basis.a, graphics=graphics.on)
			# if (graphics.on > 1) plot(as.vector(pars$W.ahp %*% basis.a), t='l')
		# }	
	} else {
		## if we don't initialize, then initial parameter guesses must be provided
		# print('scaling initial parameters... ')
		pars.init <- pars
	}
	
	if (sim.spikes) {
		if (!('W.ahp' %in% names(pars.init))) {
			pars.init$W.ahp <- rep(0, n.basis) # if W.ahp is missing
		} else {
			if (is.null(pars.init$W.ahp)) pars.init$W.ahp <- rep(0, n.basis) # if W.ahp is missing
		}
	}

	## 2. preinit: optimize the parameters of the after-spike potential - fast and easy

	if (preinit & sim.spikes){
		print('pre - initialising W.ahp ... ')
		args.parvec <- make.args.parvec(pars=pars.init, v=vsub, train$X, dt, X.ahp=X.ahp, optimvec=c('W.ahp'), regpars=regpars)
		parvec <- args.parvec$parvec; args.train <- args.parvec$args
		test.args(parvec, args.train)

		# gv.train.init <- grad.err.hGLM(parvec, args.train)
		# err.hGLM(parvec, args.train)
		p.opt <- optim(parvec, err.hGLM, grad.err.hGLM, args=args.train, method='BFGS', control=list(maxit=25, trace=10))
		presp <- err.hGLM(p.opt$par, args.train, ret.parlist=T)
		pars.init <- presp
		# print(pars.init)
	}

	args.parvec <- make.args.parvec(pars=pars.init, v=vsub, train$X, dt, X.ahp=X.ahp, regpars=regpars)
	parvec <- args.parvec$parvec; args.train <- args.parvec$args
	test.args(parvec, args.train)
	# err.hGLM(parvec, args.train)
	# gv.train.init <- grad.err.hGLM(parvec, args.train)

	if (graphics.on>0){
		plot(train$v, t='l', col=1, lwd=3)
		v.train.init <- err.hGLM(parvec, args.train, ret.resp=T)
		lines(v.train.init, col=3)
		legend('topright', legend=c('test signal', 'initial response', 'final response'), col=c(1, 3, 2), lty=1, lwd=c(3,1,1), bg=hsv(1,0,1,0.8))
	}

	print(err.hGLM(parvec, args.train))
	# grad.err.hGLM(parvec, args.train)

	##################################################################################
	print('Initial parameters tested, starting optimization ...')
	t0 <- proc.time()[1]
	p.opt <- optim(parvec, err.hGLM, grad.err.hGLM, args=args.train, method='BFGS', control=list(maxit=maxit, trace=10, abstol=1, reltol=1e-8))
	t1 <- proc.time()[1]
	t.opt <- ceiling(t1-t0)
	##################################################################################

	print('Optimization done, evaluating the results... ')
	print(paste('function: ', p.opt$counts[1], ' gradient: ', p.opt$counts[2]))
	print(paste('convergence: ', p.opt$convergence))
	
	##################################################################################
	#### Evaluate the training error and get the optimal response 
	parvec.opt <- p.opt$par
	pars.opt.sub <- err.hGLM(parvec.opt, args.train, ret.parlist=T)
	
	v.train.opt <- err.hGLM(parvec.opt, args.train, ret.resp=T)
	if (graphics.on>0){
		lines(v.train.opt, t='l', col=2)
		legend('topright', legend=c('test signal', 'initial response', 'final response'), col=c(1, 3, 2), lty=1, lwd=c(3,1,1), bg=hsv(1,0,1,0.8))
	}


	err.init.train <- err.hGLM(parvec, args.train)
	err.opt.train <- err.hGLM(parvec.opt, args.train)
	
	sig.exp <- matrix(c(init=1 - sqrt(err.init.train / LL) / sd.v, opt=1 - sqrt(err.opt.train / LL) / sd.v), 1)
	change.train.err <- 10^round(log((err.init.train/err.opt.train), 10))
	print(paste('We achieved ', change.train.err, '-fold reduction in the training rms error in ', t.opt,' seconds.', sep=''))
	print(paste('Signal explained in the training data,  initially: ', round(sig.exp[1], 4), ', after optimisation: ', round(sig.exp[2],4) , sep=''))
	
	######################################################################
	## Test error
	if (!is.null(test)){
		if (sim.spikes){
			## preparing the test data		
			print('extracting spiketimes')
			spt.test <- extract.spiketimes(test$v, sampling.freq=1000/dt, limite=-30, graphics=graphics.on)
	
			vsub.test <- as.vector(remove.spikes(test$v, spt.test, dt))
			LL.test <- sum(!is.na(vsub.test))
	
			## we need to convolve the spikes with the basis functions to get X.ahp
			X.ahp.test <- spike.basis(spt.test, L, basis.a, dt=dt)
		} else {
			vsub.test <- test$v
			X.ahp.test <- NULL
			LL.test <- length(test$v)
		}
		
		
		m.v.test <- mean(vsub.test, na.rm=T)
		sd.v.test <- sd(vsub.test, na.rm=T)
		
		args.test <- args.train
		args.test$inp <- test$X; args.test$X.ahp <- X.ahp.test; args.test$v <- vsub.test

		err.init.test <- err.hGLM(parvec, args.test)
		err.opt.test <- err.hGLM(parvec.opt, args.test)
		change.test.err <- 10^round(log((err.init.test/err.opt.test), 10))
		print(paste('We achieved ', change.test.err, '-fold reduction in the test rms error in ', t.opt,' seconds.', sep=''))

		sig.exp <- rbind(c(init=1 - sqrt(err.init.test / LL.test) / sd.v.test, opt=1 - sqrt(err.opt.test / LL.test) / sd.v.test), sig.exp)
		rownames(sig.exp) <- c('test', 'training')
		print(paste('Signal explained in the test data,  initially: ', round(sig.exp[1,1], 4), ', after optimisation: ', round(sig.exp[1,2],4) , sep=''))
	}

	
	if (sim.spikes & (sig.exp[1,2] > 0.25)){
		#######################################################################################
		## fit the spiking response
		#######################################################################################
		
		print("subthreshold membrane potential fitted, fitting the dynamical threshold ...")
	
		##################################
		## 3. gradient based maximum likelihood for the parameters of the spiking
		time <- seq(dt, length=L, by=dt)
		delay.spike <- spt$delay
		dur.spike <- spt$dur
		i.peaks <- spt$st
		i.thresholds <- i.peaks - delay.spike / dt
		i.peaks <- i.peaks[i.thresholds>0]
		i.thresholds <- i.thresholds[i.thresholds>0]
		fit.dynamic <- F
		
		## 3.1. To encourage short refractoryness, we use a slightly different set of basis functions for the threshold 
		
		n.basis.t <- min(2 + floor(length(spt$st)/20), 8)
		basis.t <- gen.basis(Tmax=200, dt=dt, graphics=graphics.on, n.basis=8, dur=spt$dur, first=T)	
		basis.t <- matrix(basis.t[1:n.basis.t,], n.basis.t)
		i.basis <- which(colSums(basis.t) > 1/1000)
		basis.t <- basis.t[,i.basis]
		L.basis <- ncol(basis.t)
		if (graphics.on > 2) matplot(t(basis.t), t="l", col=rainbow(n.basis), lty=1)
	
		# 3.1.1. convolve the spikes with the NEW basis functions
		sp.basis.th <- spike.basis(spt, L, basis.t, dt=dt)
		if (graphics.on > 2) matplot(t(sp.basis.th), t="l", col=rainbow(n.basis), lty=1)
	
		# 3.2.1. Evaluate the model to get the subthrehold responses GIVEN the observed spikes
		VV <- v.train.opt
	
		
		# 3.2. Extract the voltages before spikes (threshold crossings) and the voltage at silence
		Z <- rbind(VV, rep(-1, length(VV)), -sp.basis.th) # build matrix of the regressors
		Z.spike <- Z[,i.thresholds] # variables at THRESHOLD
		ZZ <- cut.spikes(Z, i.thresholds, dur.spike/dt)
	
		Z.static <- rbind(VV, rep(-1, length(VV)))
		Z.static.spike <- Z.static[,i.thresholds]
		ZZ.static <- cut.spikes(Z.static, i.thresholds, dur.spike/dt)
			
		# 3.3 Learn the parameters of the model - static threshold
		theta.0 <- c(1/5, mean(train$v) /5 - log(1/1000)) # beta=1/5, p=1/1000 5Hz baseline firing
		tol.theta <- 1e-4
		max.iter <- 10000
		apply(ZZ.static, 1, sd)
		theta.static <- compute.theta(ZZ.static, Z.static.spike, theta.0, tol.theta, max.iter)
		
		
		# 3.4 Learn the parameters of the dynamic threshold model
		# - initial value is from the static threshold
		theta.0 <- c(theta.static$theta, rep(0,n.basis.t))
		theta <- compute.theta(ZZ, Z.spike, theta.0, tol.theta, max.iter, eta=1/10)
		th <- theta$theta
		beta <- th[1]
		v.th <- (th[2] + log(dt/1000))/beta
		w.th <- th[3:(n.basis.t+2)]/beta
		th.opt <- as.vector(w.th %*% basis.t)
		# th.opt[1:(round(dur.spike / dt))] <- Inf # ablsolut refractoriness	
		if (graphics.on>2) plot(th.opt, t="l")
		pars.opt.spike = list(beta=beta, v.th=v.th, w.th=w.th)

		#############################################
		## testing spike timing accuracy
		#############################################
		if (!is.null(test)){
			print("model fitting DONE, testing the spike timing accuracy.")
		
			## first the reference response
			if (generate.train) {
				resp2 <- sim.hGLM.sp(test$X, dt, pars.sub, pars.spike, rseed=seed+2, double=double, Nrep=50)
				V.ref <- resp2$v + resp2$sp
			} else {
				V.ref <- matrix(test$v, 1)
			}
			
			## second the optimal prediction	
			resp.opt <- sim.hGLM.sp(test$X, dt, pars.opt.sub, pars.opt.spike, rseed=seed+103, double=double, Nrep=50)
			V.opt <- resp.opt$v + resp.opt$sp

			# resp.opt2 <- sim.hGLM.sp(test$X, dt, pars.opt.sub, pars.opt.spike, rseed=seed+203, double=double, Nrep=50)
			# V.opt2 <- resp.opt2$v + resp.opt2$sp
				
			sp.rel <- spike.reliability(V.ref, V.opt, dt, graphics=graphics.on, 0, 5000, 4)
			# sp.rel2 <- spike.reliability(V.opt, V.opt2, dt, graphics=graphics.on, 0, 5000, 4)
		
			if (graphics.on > 0) {
				spt <- extract.spiketimes(test$v, sampling.freq=1000/dt, limite=-30, graphics=F)
				V.opt <- resp.opt$v[1,] + as.vector(resp.opt$sp[1,])
				plot(test$v, t="l", lwd=3)
				lines(time, V.opt, col=2)
				
				resp.opt <- sim.hGLM.sp(test$X, dt, pars.opt.sub, pars.opt.spike, rseed=seed+3, ref.spikes=spt$st*dt - spt$del, double=double, Nrep=1)
				V.opt2 <- resp.opt$v + as.vector(resp.opt$sp)
		
				lines(V.opt2[1,], col=3)
				legend("topright", leg=c("true", "sub + spikes", "sub"), lty=1, col=c(1,2,3), bg=grey(1, alpha=0.75))
			}
			out <- list(s.exp=sig.exp, sp.rel=sp.rel, pars.sub=pars.opt.sub, pars.spike=pars.opt.spike, convergence=p.opt$convergence)
		} else {
			out <- list(s.exp=sig.exp, pars.sub=pars.opt.sub, pars.spike=pars.opt.spike, convergence=p.opt$convergence)			
		}
		
	} else {
		out <- list(s.exp=sig.exp, pars.sub=pars.opt.sub, convergence=p.opt$convergence)
	}
	out
}

# ############################################
# diff.pars.ahp <- function(p.target, p.opt, Wc.e, Wc.i){
	# ## For time constants the distance in the log(tau) is used. 
	# ## For all other parameters distance in the natural place is considered.
	# d.pars <- (p.opt[1] - p.target$v0)^2
	# pp.opt <- p.opt[-1]
	# M <- length(p.target$Jw)	
	# ## certain parameters are irrelevant (time constant and synaptic weight of sobunits not used.)
	# i.pars.e <- !unlist(lapply(Wc.e, is.null))
	# i.pars.i <- !unlist(lapply(Wc.i, is.null))
	# d.pars <- d.pars + sum(((exp(pp.opt[1:M]) - exp(p.target$Jw))^2)[i.pars.e])
	# pp.opt <- pp.opt[-(1:M)]
	# d.pars <- d.pars + sum(((pp.opt[1:M] - p.target$Ww.e)^2)[i.pars.e])
	# pp.opt <- pp.opt[-(1:M)]
	# d.pars <- d.pars + sum(((pp.opt[1:M] - p.target$Ww.i)^2)[i.pars.i])
	# pp.opt <- pp.opt[-(1:M)]
	# d.pars <- d.pars + sum((pp.opt[1:M] - p.target$Th)^2, na.rm=T)
	# pp.opt <- pp.opt[-(1:M)]
	# d.pars <- d.pars + sum(((pp.opt[1:M] - p.target$Tau.e)^2)[i.pars.e])
	# pp.opt <- pp.opt[-(1:M)]
	# d.pars <- d.pars + sum(((pp.opt[1:M] - p.target$Tau.i)^2)[i.pars.i])
	# d.pars
# }

# # # 1L
# pars <- list(Jc=0, Wc.e=list(seq(1,20)), Wc.i=list(seq(21,40)))
# parserr <- train.hLN.ahp(dt=0.2, pars=pars, seed=3, graphics.on=T, Tmax=20000, linear.soma=T)

# # # # ## 2N
# pars <- list(Jc=c(0,1,1), Wc.e=list(NULL, seq(1,10), seq(11,20)), Wc.i=list(NULL, seq(21,30), seq(31,40))); 
# parserr <- train.hLN.ahp(dt=0.2, pars=pars, seed=3, graphics.on=T, Tmax=5000)

# parserr <- train.hLN.ahp(pars=pars, seed=4, graphics.on=T, Tmax=5000, linear.soma=T)

