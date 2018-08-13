
errsin <- function(pars, data){
	## calculates the error between the datapoints and 
	## a sinusoid function used to model the data defined by its parameters:
	## phase (phi), amplitude ("ampl) and minimum ("min"); (frequency = 1/pi)
	oris <- seq(0, by=90/4, length=16)
	p.up <- (sin((oris-pars["phi"])*2*pi*2/360) + 1) /2 * pars["ampl"] + pars["min"]
	P.up <- matrix(rep(p.up, ncol(data)), ncol=ncol(data), byrow=F)
#	P.up <- cbind(p.up, p.up, p.up, p.up, p.up, p.up)
	err <- sum((data - P.up)^2)
	err
}



er.p <- function(m.alpha, p.up, beta=10){
	## calculates the max.alpha parameters for the stimulus orientations
	## alpha = (sin(phi) +1) * max.alpha / 2
	### max.alpha depends on stimulus orientation
	### for each orientation, max.alpha is fitted to match the total up state probability
	t <- seq(0, 3042)
	alphas <- (sin(t/500 *2*pi+150) + 1) * m.alpha / 2	
	err <- (mean(alphas/(alphas+beta)) - p.up)^2
	err
}



gen.events.sin <- function(Tmax, m.alpha, beta, maxL=NULL){
	## T: time in ms
	## m.alpha: max of the down-to up transition rate in Hz
	## beta: up - down transition rate in Hz
	## maxL: maximum duration of an up state
	
	m.alpha <- m.alpha / 1000; beta <- beta / 1000 # conversion from Hz to 1/ms
	max.rate <- max(m.alpha, beta)
	t.events <- gen.Poisson.events(Tmax, max.rate)
	state <- 0 # start from down state
	st <- rep(0, Tmax)
	t.events.kept <- c(0, 0) # time (ms), event type (down - 0; up - 1)

	if (length(t.events > 0)){
		for (tt in t.events){
			## we keep the transition with some probability
			if (state==0) rr <- (sin(tt/500 *2*pi+150) + 1) * m.alpha / 2 else rr <- beta
			p.keep <- rr/max.rate
			if (runif(1) < p.keep){
				# print("transition kept")
				state <- 1-state
				st[round(tt):Tmax] <- state
				if (!is.null(maxL)){
					if (state == 0){
						t.last <- tail(t.events.kept,1)[2]
						Li <- tt-t.last
						if (Li > maxL){
							t1 <- runif(1, t.last + Li/5, t.last + 2*Li/5)
							t2 <- runif(1, t.last + 3*Li/5, t.last + 4*Li/5)
							t.events.kept <- rbind(t.events.kept, c(0, t1))
							t.events.kept <- rbind(t.events.kept, c(1, t2))
							st[round(t1):round(t2)] <- 0							
						}
					}
				}
				t.events.kept <- rbind(t.events.kept, c(state, tt))
			} else {
				# print ("rejected")
			}
		}
	}
	attributes(st) <- list(spt = t.events.kept)
	st
}

gen.spikes.states <- function(Tmax, N, mu, sd, tau, x, graphics=F){
	# Tmax, dt - ms
	# N: number of neurons
	# mu: vector of firing rates (Hz)
	# sd: sd of firing rates (Hz)
	# tau: time constant (ms)
	# x: states 1 - 2
	## Tmax <- 50000; N <- 100; mu <- c(5, 5); sd <- c(2, 2); tau <- 65; x <- gen.events.sin(Tmax, 2, 20, maxL = 150) + 1
	## N <- Ensyn[i.den]; mu <- Erate; sd <- Esd; tau <- 500; x <- st+1; graphics <- T
	dt.rates <- 1
	L <- Tmax / dt.rates; 

  	## the firing rates
	# Q <- sqrt(sd^2 * 2 / tau)
	rates <- rep(mu[x[1]], L)
	z <- rep(0, L) # we simulate an OU process with sd=1 and add it to the rates
	z0 <- 0
	Q <- sqrt(2/tau)
	
	for (i in 2:L){
		z[i] <- z[i-1] + (z0 - z[i-1]) * dt.rates / tau + Q * rnorm(1,0,1) * sqrt(dt.rates)
		r0 <- mu[x[i]]
		rates[i] <- r0 + sd[x[i]] * z[i]
	}
	
	if (graphics){
		plot(rates, t="l")
		lines(x*10-5, col=2)
	}
	
	## the spikes
  	rates[rates<0] <- 0
	t.sp <- gen.Poisson.train.rescale(rates, N)
	t.sp
}


gen.Poisson.events <- function(Tmax, rate){
	## homogeneous Poisson process
	## Tmax: time in ms
	## rate: event frequency in 1/ms
	t0 <- rexp(1, rate)
	if (t0 > Tmax) {
		sp <- NULL
	} else {
		sp <- t0
		tmax <- t0
		while(tmax < Tmax){
			t.next <- tmax + rexp(1, rate)
			if (t.next < Tmax) sp <- c(sp, t.next)
			tmax <- t.next
		}
	}
	return(sp)		
}

gen.Poisson.train <- function(rates, N){
	## inhomogeneous Poisson process
	Tmax <- length(rates) # Tmax in ms
	max.rate <- max(rates) # rates in Hz
	t.sp.kept <- c(0, 0) # cell, time (ms)
	
	for (cell in 1:N){	
		t.sp <- gen.Poisson.events(Tmax, max.rate/1000)
		if (length(t.sp > 0)){
			for (tt in t.sp){
				## we keep the transition with some probability
				rr <- rates[ceiling(tt)]
				p.keep <- rr/max.rate # both in Hz!
				if (runif(1) < p.keep){
					t.sp.kept <- rbind(t.sp.kept, c(cell-1, tt))
				}
			}
		}
	}
	t.sp.kept <- t.sp.kept[-1,]
	ii <- sort(t.sp.kept[,2], index.return=T)$ix
	t.sp <- t.sp.kept[ii,]
	t.sp
}


gen.Poisson.train.rescale <- function(rates, N){
	## inhomogeneous Poisson process using time rescaling theorem
	## idea from Fig 1B of Pillow: Time-rescaling methods for the estimation and assessment of non-Poisson neural encoding models
	## should be very fast!

	Tmax <- length(rates) # Tmax in ms
	mean.rate <- mean(rates) # rates in Hz
	t.spikes <- c(0, 0)	# cell, time (ms)
	
	for (cell in 1:N){	
		t.sp <- gen.Poisson.events(Tmax, mean.rate/1000) # convert from 1/s to 1/ms
		if (length(t.sp > 0)){
			T.sp <- t.sp * mean.rate
			
			Cr <- cumsum(rates)
			t <- seq(1, Tmax)
			t.rescaled <- approx(x=Cr, y=t, xout=T.sp)$y
			t.sp.cell <- cbind(rep(cell-1, length(t.rescaled)), t.rescaled)
			t.spikes <- rbind(t.spikes, t.sp.cell)
		}
	}
	t.spikes <- t.spikes[-1,]
	ii <- sort(t.spikes[,2], index.return=T)$ix
	t.sp <- t.spikes[ii,]
	t.sp
}


gen.Poisson.spikes <- function(r){
	## r is a vector of rates (Hz), spike counts are generated with the same resolution + jitter
	sp <- c(0,0)
	L <- length(r)
	for (i in 1:L){
		n.sp <- rpois(1, r[i]/1000)
		if (n.sp>0)	sp <- rbind(sp, cbind(  rep(0, n.sp), i + sort(runif(n.sp, 0, 1)) -1 ) )
	}
	sp <- sp[-1,]
}

up.state.duration <- function(resp, l.filt = 5, th=-35, graphics=F){
	N <- ncol(resp)
	dur.ups <- NA

	for (i in 1:N){
		filt <- dnorm(seq(-5*l.filt,5*l.filt), 0, l.filt); filt <- filt / sum(filt)
		rr <- filter(resp[,i], filt)
		LL <- length(rr)
		ups <- which((rr[-1]>th) & (rr[1:(LL-1)]<th))
		downs <- which((rr[-1]<th) & (rr[1:(LL-1)]>th))
		
		if (graphics){
			plot(resp[,i], t="l")
			lines(rr, col=2)
			abline(h=th)
			abline(v=ups, col=3)
			abline(v=downs, col=4)
		}

		for (j in 1:length(ups)){
			t.up <- ups[j]
			dd <- which(downs>t.up)
			if (length(dd) > 0){
				t.down <- downs[min(dd)]
				v.up <- rr[t.up:t.down]
				if (max(v.up) > (th+10)){
					dur.ups <- c(dur.ups, t.down-t.up)
				}
			}
		}
	}
	dur.ups <- dur.ups[!is.na(dur.ups)]
	dur.ups
}


add_gamma <- function(rate, tau=5, Tmin=20, Tmax=60, dt=1, rate_bg=1){
	# function to modulate the firing rate to have bursts of high activity
	# rate: firing rate as the function of time, Hz
	# dt: resolution of the firing rate data, in ms
	# tau: length of gamma bursts, in ms
	# Tmin: minimum period of gamma bursts, in ms
	# Tmax: maximum period of gamma bursts, in ms
	# rate_bg: background rate, Hz
	#
	# rate <- rates; tau <- 5; Tmin <- 20; Tmax <- 60; dt <- 1; rate_bg <- 5
	
	L <- length(rate)
	L_burst <- tau / dt
	rate_new <- rep(rate_bg, L)

	i_first <- 1
	LL <- round(runif(1, Tmin, Tmax) / dt)
	i_last <- min(i_first + LL - 1 , L)
	if ((L-i_last) < Tmin) {
		i_last <- L
		LL <- i_last - i_first + 1
	}
	

	while(i_last < L){
		sum_rate <- sum(rate[i_first:i_last])
		if (sum_rate < rate_bg * LL) {
			rate_new[i_first:i_last] <- rate[i_first:i_last]
		} else {
			ii_first <- i_last - L_burst + 1
			
			rate_new[ii_first:i_last] <- rate_new[ii_first:i_last] + (sum_rate - rate_bg * LL) / L_burst
		}

		i_first <- i_last + 1
		LL <- round(runif(1, Tmin, Tmax) / dt)
		i_last <- min(i_first + LL - 1 , L)
		if ((L-i_last) < Tmin) {
			i_last <- L
			LL <- i_last - i_first + 1
		}
	}
	
	sum_rate <- sum(rate[i_first:i_last])
	if (sum_rate < rate_bg * LL) {
		rate_new[i_first:i_last] <- rate[i_first:i_last]
	} else {
		ii_first <- i_last - L_burst + 1
		rate_new[ii_first:i_last] <- rate_new[ii_first:i_last] + (sum_rate - rate_bg * LL) / L_burst
	}
		
	graphics <- F
	if (graphics == T){
		plot(rate_new, t='l')
		lines(rates, col=2)
	}
	if (sum(rate_new) != 	sum(rate)) warning('number of spikes have changed, check the function add_gamma()')

	rate_new
	
}