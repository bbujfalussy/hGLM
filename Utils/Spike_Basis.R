# t <- seq(0, 100, len=1000)
# a <- 3.75
# c <- 0.01
# M <- 10
# phis <- seq(-2,14, len=M)

basf <- function(t=seq(1, 1000, len=1000), a=3.75, c=0.01, phis=seq(3,22, len=10), monoton=F, first=T){
	# low-level function - Generates cosyne bump basis functions
	#
	# Args:
	# t: time (ms)
	# a, c, phis: parameters of the basis functions - don't change them!
	# monoton: only the monotonic first part is used
	# first: the first value of the first basis should be 0 (F) or 1 (T)?
	# 
	# Returns: 
	# matrix, with rows being the basis functions

	M <- length(phis)
	mat.f <- matrix(0, length(phis), length(t))
	for (i in 1:length(phis)){
		phi <- phis[i]
		x <- a * log(t + c)
		f <- 0.5 * cos(x-phi) + 0.5
		f[x>phi+pi] <- 0
		f[x<phi-pi] <- 0
		if ((first) & (i==1)) {
			i.max <- which.max(f)
			f[1:i.max] <- max(f)
		}
		if (monoton){
			if (i>1){
				i.max <- which.max(f)
				f[1:i.max] <- max(f)			
			}
		}
		mat.f[i,] <- f		
	}
	mat.f
}

# matplot(t(basf()), t="l", col=rainbow(10), lty=1)

###################################################################
gen.basis <- function(Tmax=200, dt=1, n.basis=10, dur=1, first=F, graphics=0){
	# a more user friendly way to generate basis functions
	L.basis <- Tmax/dt
	T.basis <- seq(dt, Tmax, by=dt)
	phi.min <- ceiling(3.75 * log(dur + 0.01)) 	# the minimum duration of the spike
	# n.basis <- ceiling((17-phi.min) / (pi/2))
	# phi.min <- 1
	basis.a <- basf(T.basis, phis=seq(phi.min, 17, len=n.basis), first=first)
	if (graphics>1) matplot(T.basis, t(basis.a), t="l", pch=21, col=rainbow(n.basis*1.3), lty=1, log="x")
	basis.a
}

spike.basis <- function(spt, L, basis.a, dt=1){
	# for loop to calculate the basis activations from spike times

	i.ths <- round(spt$st - spt$delay/dt)
	n.basis <- nrow(basis.a)
	L.basis <- ncol(basis.a)
	sp.basis <- matrix(0, n.basis, L + L.basis)	
	for (i.th in i.ths){
		inds <- seq(max((i.th+1),1), (i.th+L.basis))
		# print(length(inds))
		sp.basis[,inds] <- sp.basis[,inds] + basis.a[,(L.basis - length(inds)+1):L.basis]
	}
	sp.basis <- sp.basis[,1:L]
	sp.basis
}

default.spikes <- function(dt=1, n.basis=8, Tmax=200, v.soma, pars.sub, graphics=0, rand.pars=F, seed=12, rate=2){
	# set default spike parameters to show adaptation and refractoriness
	# rate: the desired firing rate - without adptation and threshold modification, assuming Gaussian Vm
	linear.soma <- F
	if (is.na(pars.sub$Th[1])) linear.soma <- T
	if (linear.soma) y.soma <- (v.soma - pars.sub$v0)/ exp(pars.sub$Jw[1]) else y.soma <- isigm((v.soma - pars.sub$v0)/exp(pars.sub$Jw[1]), c=pars.sub$Th[1])
	
	mu.v <- mean(v.soma)
	sd.v <- sd(v.soma)

	mu.y <- mean(y.soma)
	sd.y <- sd(y.soma)
	# m.2sd <- mu.v + 2.5*sd.v

	# if (dt == 3) l.mu <- log(2/dt) else l.mu <- log(3/dt)
	# if (dt==300) l.sd <- log(200/dt) else l.sd <- log(300/dt)
	# v.th <- (l.sd*mu.v - l.mu*m.2sd)/(l.sd-l.mu)
	# beta <- l.mu / (mu.v-v.th)
	beta <- 1 / sd.v
	# v.th <- -1 * log(rate * exp(-beta * mu.v - beta^2 * sd.v^2/2)) / beta # simplify it
	v.th <- mu.v + beta * sd.v^2/2 - log(rate) / beta
	
	w.th <- c(2*sd.v, 2*sd.v, seq(2*sd.v, 0, length=n.basis-4), (-1) * sd.v, 0) 
	W.ahp <- seq(-2*sd.y, 0, length=n.basis)
	# w.th <- rep(0, n.basis)
	# W.ahp <- rep(0, n.basis)


	if(graphics>1){		
		# par(mfcol=c(2,1))
		vv <- hist(v.soma, 100, freq=F)$mids
		lines(vv, dt/1000 * exp(beta*(vv-v.th)), col=2)
		 # cat(round(v.th, 3), " ", round(beta,3), "\n")
		
		basis.t <- gen.basis(200, dt, n.basis, first=T)
		basis.a <- gen.basis(200, dt, n.basis, first=F)
		L.basis <- ncol(basis.a)
		T.basis <- seq(dt, Tmax, by=dt)


		th <- as.vector(w.th %*% basis.t)
		adapt <- as.vector(W.ahp %*% basis.a)

		matplot(T.basis, cbind(th, adapt), t="l", pch=21, col=c(1,2), log="x")
	}
	pars.spike <- list(beta=beta, v.th=v.th, w.th=w.th, W.ahp=W.ahp)
	pars.spike
}


init.ahp <- function(v, dt, spt, basis.a, graphics=0){
	## a function to initialise the W.ahp parameters to match observed data
	## - uses a linear regression between mean spike shape and the basis activations
	
	st <- spt$st
    nbr.spikes <- length(st)
	size.spike.shape <- 30/dt;          # size of the spike shape (30 ms)
	
	temp.spikes <- st - spt$del/dt             # set the spiketimes 5 ms before the maximum of the AP
	spike.shape <- matrix(NA, size.spike.shape,nbr.spikes)     # i.e. just used for plotting
	k<-1
	for (i in 1:(nbr.spikes-1)){
		# only spikes not followed by another spike within the test period ...
	    if ((temp.spikes[i] + size.spike.shape <= temp.spikes[i+1]) & (temp.spikes[i]>0)){ 
		    	spike.shape[,k] <- v[temp.spikes[i]:(temp.spikes[i]+size.spike.shape-1)]
	        k <- k+1
	    }
	}
		
	m.spike.shape <- rowMeans(spike.shape, na.rm=T)             # Compute the mean spike shape
	m.spike.shape <- m.spike.shape - m.spike.shape[1]
	m.spike.shape <- m.spike.shape[-1]

	l.shape <- length(m.spike.shape)
	n.bas <- sum(apply(basis.a, 1, which.max) < l.shape)
	basis.aa <- t(basis.a[1:n.bas, 1:l.shape])
	
	n.cut <- spt$dur/dt - 1
	
	mm.shape <- m.spike.shape[-c(1:n.cut)]
	basis.aa <- basis.aa[-c(1:n.cut),]
	
	ww.ahp <- lm(mm.shape ~ 0 + basis.aa)$coeff
	w.ahp <- as.vector(c(ww.ahp, rep(0, nrow(basis.a) - n.bas)))
	if (graphics>1){
		plot(mm.shape, t="l", ylim=c(min(mm.shape), 0))
		matplot(basis.aa*(-1), t="l", col=rainbow(12), add=T, lty=1)
		lines(basis.aa %*% ww.ahp, col=2, lwd=3, lty=2)
	}	
	
	w.ahp
}
