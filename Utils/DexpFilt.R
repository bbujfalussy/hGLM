#################################################
## filters the input spike train with a double exponential kernel defined by the two time constants tau1 and tau2:
## f(t) = exp(-t/tau2) - exp(-t/tau1)
## if delay.t > 0 it also delays the spikes with delay.t: f(t) = exp(-(t-delay.t)\tau2) - exp(-(t-delay.t)\tau1)
## the input (spikes) can be vector or matrix, not necessary binary
## if matrix, than each row is filetered independently using the same taus
## the output is a vector of the same length as the input
## the parametrisation of the function is a bit unusual: tau2 = tau1 + dtau
## in this case both tau1 and dtau has to be positive 


dexp.filt <- function(tau1, dtau, spikes,  delay.t=0, dt=0.5, grad.tau=F, grad.delay=F, run.for=NULL){
	if (is.na(tau1))	{
		# print("Tau is NA")
		return (-spikes)
	}
	if (is.nan(tau1)) {
		# print("Tau is NAN")
		return (-spikes)
	}
	if (tau1 <= 1e-10) {	
		# print("Tau is too small"); 
		return (-spikes)
	}
	if (tau1 > 500) {
		# print("Tau is too large"); 
		return (-spikes)
	}
	if (is.na(dtau))	{
		# print("Tau is NA")
		return (-spikes)
	}
	if (is.nan(dtau)) {
		# print("Tau is NAN")
		return (-spikes)
	}
	if (dtau <= 1e-10) {	
		# print("dTau is too small"); 
		return (-spikes)
	}
	if (dtau > 500) {
		# print("Tau is too large"); 
		return (-spikes)
	}

	tau2 <- tau1 + dtau

	if ((grad.tau>0) & (grad.delay)) stop("only one gradient please!")
	if (delay.t < 0) delay.t <- 0; if (delay.t > tau2) {
		warning("propagation delay is too long!")
		if (delay.t > 2*tau2) delay.t <- 2 * tau2
	}
	if (is.vector(spikes)) spikes <- matrix(spikes, 1, length(spikes))
	N <- nrow(spikes); Lmax <- ncol(spikes)
	Tmax <- length(spikes) * dt; t <- seq(dt,Tmax,dt)
	# if (grad.tau>0) tf <- seq(0, 12 * tau2, by=dt) else tf <- seq(0, 7 * tau2, by=dt); 
	tf <- seq(0, 10 * tau2, by=dt)
	L <- length(tf)
	filt <- exp(-(tf-delay.t)/tau2) -  exp(-(tf-delay.t)/tau1); filt[filt<0] <- 0
	f0 <- as.numeric(filt>0)
	# we compute the gradient of the error wrt tau1
	if (grad.tau == 1) filt <- f0 * ((tf-delay.t) / tau2^2 * exp(-(tf-delay.t)/tau2) - (tf-delay.t) / tau1^2 * exp(-(tf-delay.t)/tau1) )
	# we compute the gradient of the error wrt dtau
	if (grad.tau == 2) filt <- f0 * ((tf-delay.t) / tau2^2 * exp(-(tf-delay.t)/tau2)) # only if we compute the gradient of the error wrt dtau!
	if (grad.delay == T) { # only if we compute the gradient of the error wrt delay.t!
		filt <- f0 * (1/tau2 * exp(-(tf-delay.t)/tau2) - 1 / tau1 * exp(-(tf-delay.t)/tau1))
	}
	filt[is.nan(filt)] <- 0
	if (is.null(run.for)) {
		if ((sum((spikes)>0) / length(spikes)) > 0.1) run.for <- F else run.for <- T
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
	} else {
		f0.spikes <- cbind(matrix(0, N, L), spikes)
		f.spikes <- t(filter(t(f0.spikes), filt, sides=1))
		f.spikes <- f.spikes[,-(1:L)]		
	}
	f.spikes
}


# ####################################
# test the evaluation of the function
# spikes <- rep(0, 2000); spikes[c(120,550,900)] <- c(1,1.2, 1.5)
# plot(spikes, t="h")
# lines(spikes, t="h", col=3)
# tau1 <- 50; dtau=100
# fs <- dexp.filt(tau1, dtau, spikes, delay.t=10)
# lines(fs, t="l", col=2)

# plot(spikes, t="h")
# tau1 <- 2; dtau=70
# fs <- dexp.filt(tau1, dtau, spikes, delay.t=0.10)
# lines(fs, t="l", col=2)



# fs <- dexp.filt(tau1, dtau, spikes, grad.tau=1)
# plot(fs, t="l")
# lines(spikes, t="h", col=3)
# fs <- dexp.filt(tau1, dtau, spikes, delay.t=10, grad.tau=1)
# lines(fs, t="l", col=2)

# fs <- dexp.filt(tau1, dtau, spikes, grad.tau=2)
# plot(fs, t="l")
# lines(spikes, t="h", col=3)
# fs <- dexp.filt(tau1, dtau, spikes, delay.t=10, grad.tau=2)
# lines(fs, t="l", col=2)

# fs <- dexp.filt(tau1, dtau, spikes, grad.delay=T)
# plot(fs, t="l")
# lines(spikes, t="h", col=3)
# fs <- dexp.filt(tau1, dtau, spikes, delay.t=10, grad.delay=T)
# lines(fs, t="l", col=2)


#################################################
## find the correct parameters with optimisation
# pars <- c(tau1=2, dtau=20, delay=2.34)
# spikes <- rep(0, 2000); spikes[c(120,550,900)] <- c(1,1.2, 1.5)
# dt <- 0.5

# resp <- dexp.filt(tau1=pars[1], dtau=pars[2], spikes, delay.t=pars[3], dt=dt)
# args <- list(dt=dt, resp=resp, spikes=spikes)

# err.par <- function(pars, args){
	# resp <- dexp.filt(tau1=pars[1], dtau=pars[2], args$spikes, delay.t=pars[3], dt=args$dt)
	# err <- sum((resp - args$resp)^2)
	# err
# }

# grad.par <- function(pars, args){
	# resp <- dexp.filt(tau1=pars[1], dtau=pars[2], args$spikes, delay.t=pars[3], dt=args$dt)
	# g.tau1 <- dexp.filt(tau1=pars[1], dtau=pars[2], args$spikes, delay.t=pars[3], dt=args$dt, grad.tau=1)
	# g.dtau <- dexp.filt(tau1=pars[1], dtau=pars[2], args$spikes, delay.t=pars[3], dt=args$dt, grad.tau=2)
	# g.del <- dexp.filt(tau1=pars[1], dtau=pars[2], args$spikes, delay.t=pars[3], dt=args$dt, grad.delay=T)
	# grad.tau1 <- sum(2 * (resp - args$resp) * g.tau1)
	# grad.dtau <- sum(2 * (resp - args$resp) * g.dtau)
	# grad.delay <- sum(2 * (resp - args$resp) * g.del)
	# gr <- c(grad.tau1, grad.dtau, grad.delay)
	# gr
# }


# grad.par(pars, args)

# pars.init <- c(tau1=1, dtau=15, delay=1)
# optim(pars.init, err.par, gr=grad.par, args=args, method="BFGS")


