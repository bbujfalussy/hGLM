
#################################################
## filters the input spike train with an alpha kernel defined by the time constant tau: f(t) = t/tau exp(-t\tau)
## the amplitude of the function is constant, 1/e, while its integral is tau
## if delay.t > 0 it also delays the spikes with delay.t: f(t) = (t-delay.t)/tau exp(-(t-delay.t)\tau)
## the input (spikes) can be vector or matrix, not necessary binary
## if matrix, than each row is filetered independently using the same tau
## the output is the filtered spike train - same length as the input!
## it turns out that for sparse spikes and relatively short filters it is much faster - this is implemented in the first function.

alpha.filt.ampl <- function(tau, spikes, delay.t=0, dt=0.5, grad.tau=F, grad.delay=F, run.for=NULL){
	if (is.na(tau))	{
		# print("Tau is NA")
		return (-spikes)
	}
	if (is.nan(tau)) {
		# print("Tau is NAN")
		return (-spikes)
	}
	if (tau <= 1e-10) {	
		# print("Tau is too small"); 
		return (-spikes)
	}
	if (tau > 500) {
		# print("Tau is too large"); 
		return (-spikes)
	}
	
	
	if ((grad.tau) & (grad.delay)) stop("only one gradient please!")
	if (delay.t < 0) delay.t <- 0
	if (delay.t > tau) {
		warning("propagation delay is too long!")
		# cat('delay: ', delay.t, ', tau: ', tau, '\n')
		if (delay.t > 2*tau) delay.t <- 2*tau # this is not nice, but otherwise we can get missing values in the filter
	}
	if (is.vector(spikes)) spikes <- matrix(spikes, 1, length(spikes))
	N <- nrow(spikes); Lmax <- ncol(spikes)
	Tmax <- length(spikes) * dt; t <- seq(dt,Tmax,dt)
	if (grad.tau) tf <- seq(0, 15 * tau, by=dt) else tf <- seq(0, 10 * tau, by=dt); 
	L <- length(tf)
	filt <- (tf-delay.t)/tau * exp(-(tf-delay.t)/tau); filt[filt<0] <- 0
	if (grad.tau == T) filt <- filt * ((tf-delay.t)/tau^2 - 1/tau) # only if we compute the gradient of the error wrt tau!
	if (grad.delay == T) { # only if we compute the gradient of the error wrt delay.t!
		f0 <- as.numeric(filt>0)
		filt <- f0 * ((tf-delay.t)/tau^2 * exp(-(tf-delay.t)/tau) - 1/tau * exp(-(tf-delay.t)/tau))
	}
	if (is.null(run.for)) {
		# cat(sum((spikes^2)>0), length(spikes), "\n")
		if ((sum((spikes^2)>0) / length(spikes)) > 0.1) run.for <- F else run.for <- T
	}

	if (run.for){ # if the for cycle seems to be faster...
		f0.spikes <- matrix(0, N, L+Lmax)
		for (n in 1:N){
			ispn <- which(spikes[n,] != 0)
			if (length(ispn > 0)){
				for (isp in ispn){
					f0.spikes[n,isp:(isp+L-1)] <- f0.spikes[n,isp:(isp+L-1)] + spikes[n,isp] * filt
				}
			}
		}		
		f.spikes <- f0.spikes[,1:Lmax]
	} else { # if the filter seems faster ...
		f0.spikes <- cbind(matrix(0, N, L), spikes)
		f.spikes <- t(filter(t(f0.spikes), filt, sides=1))
		f.spikes <- f.spikes[,-(1:L)]		
	}
	f.spikes
}


# # # ####################################
# # # test the evaluation of the function

# spikes <- rep(0, 2000); spikes[c(120,550,900)] <- c(1,1.2, 1.5)
# plot(spikes, t="h")
# tau <- 10
# dt <- 0.5
# grad.tau <- F

# fs <- alpha.filt.ampl(tau, spikes, grad.tau=F)
# plot(fs, t="l")
# lines(spikes, t="h", col=3)

# spikes <- rep(0, 200); spikes[c(12,55,90)] <- c(1,1.2, 1.5)
# fs <- alpha.filt.ampl(tau, spikes, grad.tau=F, dt=5)
# lines(seq(10, 2000, by=10), fs, t="l", col=2)


# fs <- alpha.filt.ampl(tau, spikes, delay.t=10, grad.tau=F)
# lines(fs, t="l", col=2)
# abline(h=exp(-1), col=4, lty=2)

# fs <- alpha.filt.ampl(tau, spikes, grad.tau=T)
# plot(fs, t="l")
# lines(spikes, t="h", col=3)
# fs <- alpha.filt.ampl(tau, spikes, delay.t=10, grad.tau=T)
# lines(fs, t="l", col=2)

# fs <- alpha.filt.ampl(tau, spikes, grad.delay=T)
# plot(fs, t="l")
# lines(spikes, t="h", col=3)
# fs <- alpha.filt.ampl(tau, spikes, delay.t=10, grad.delay=T)
# lines(fs, t="l", col=2)

# #################################################
# ## find the correct parameters with optimisation
# pars <- c(tau=20, delay=2.34)
# spikes <- rep(0, 2000); spikes[c(120,550,900)] <- c(1,1.2, 1.5)
# dt <- 0.5

# resp <- alpha.filt.ampl(tau=pars[1], spikes, delay.t=pars[2], dt=dt)
# args <- list(dt=dt, resp=resp, spikes=spikes)

# err.par <- function(pars, args){
	# resp <- alpha.filt.ampl(tau=pars[1], args$spikes, delay.t=pars[2], dt=args$dt)
	# err <- sum((resp - args$resp)^2)
	# err
# }

# grad.par <- function(pars, args){
	# resp <- alpha.filt.ampl(tau=pars[1], args$spikes, delay.t=pars[2], dt=args$dt)
	# g.tau <- alpha.filt.ampl(tau=pars[1], args$spikes, delay.t=pars[2], dt=args$dt, grad.tau=T)
	# g.del <- alpha.filt.ampl(tau=pars[1], args$spikes, delay.t=pars[2], dt=args$dt, grad.delay=T)
	# grad.tau <- sum(2 * (resp - args$resp) * g.tau)
	# grad.delay <- sum(2 * (resp - args$resp) * g.del)
	# gr <- c(grad.tau, grad.delay)
	# gr
# }


# grad.par(pars, args)

# pars.init <- c(tau=5, delay=1)
# optim(pars.init, err.par, gr=grad.par, args=args, method="BFGS")


### double exponential approximation

# # L <- 1000; dt <- 0.5; t <- seq(dt, by=dt, length=L)
# spikes <- rep(0, L); spikes[50] <- exp(1)
# grad.tau <- F
# fs <- alpha.filt.ampl(tau, spikes, grad.tau=F)
# plot(t, fs, axes=F, t="n"); axis(1); axis(2, las=2)

# source('~/Programs/hGLMs/Rlib/DexpFilt.R', chdir = TRUE)

# for (tau in c(2, 4, 10, 20, 50)){
	
	# fs <- alpha.filt.ampl(tau, spikes, grad.tau=F)
	# lines(fs, t="l")
	
	# dtau <- tau / 1000
	# tau1 <- tau - dtau
	# tpeak <- tau1 * (tau1 + dtau) / dtau * log((tau1+dtau)/tau1)
	# ampl <- 1 / (exp(-tpeak / (tau1 + dtau)) - exp(-tpeak / (tau1)))
	
	# fs2e <- dexp.filt(tau1, dtau, spikes * ampl * exp(-1))
	# lines(fs2e, t="l", col=2, lty=2)
# }



