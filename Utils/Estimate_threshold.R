## we take P(s=0) = exp(-r) and 
## P(s>0) = 1 - exp(-r) ~ r*exp(-r/2) - this is much better approximation than the usual P(s>0) ~ r
## we derived it by having a second order Taylor expansion for the spiking - the first order was REALLY biased!
## but I could not find a way to reproduce this derivation

# par(mfcol=c(1,2))
# r <- seq(1/1000, 1, length=1000) # rate
# plot(r, log(1-exp(-r)), t="l", ylim=c(-7, 0)) # true log P(s>0); the true P(s=0) = exp(-r) if spiking is Poisson
# lines(r, log(r), col=2) # first order
# lines(r, log(r) - r/2, col=3) # second order

# plot(r, 1-exp(-r), t="l", ylim=c(0,1)) # true log P(s>0); the true P(s=0) = exp(-r) if spiking is Poisson
# lines(r, r, col=2) # first order
# lines(r, r*exp(- r/2), col=3) # second order


# r <- seq(1/1000, 1, length=1000) # rate
# plot(r, r, t="l", log="", col=2) # true - Bernoulli P(s=1)
# lines(r, (1-exp(-r)), col=1) # first order - Poisson approximation for P(s=1)
# lines(r, exp(log(r) -r/2), col=3) # second order ?

compute.theta <- function(X, X.spike, theta.0, tol.theta, max.iter, eta=0){
	# columns: time steps
	# rows: regressors
	# eta: prior on the higher terms, not on theta[1:2]
	
	eta <- c(0, 0, rep(eta, length(theta.0)-2))
	if (is.vector(theta.0)) theta.0 <- matrix(theta.0, 1, length(theta.0))

	lrate.spike <- theta.0 %*% X.spike           # log rate at spike
	lrate <- theta.0 %*% X	                        # log rate at silence
	rate.spike <- exp(lrate.spike) 			        # rate at silence
	rate <- exp(lrate) 			                  		# rate at silence	
	
	L <- sum(lrate.spike) -sum(rate.spike)/2 - sum(rate) - (sqrt(eta) * theta.0) %*% t(sqrt(eta) * theta.0)   # likelihood - see Dayan and Abbott E. 1.37
	G <- X.spike %*% t(1-rate.spike/2) - X %*% t(rate)  - eta * t(theta.0)       # gradient
	temp <- matrix(NA, nrow(X), ncol(X))
	for (j in 1:nrow(X)) temp[j,] <- X[j,] * rate
	temp.spike <- matrix(NA, nrow(X.spike), ncol(X.spike))
	for (j in 1:nrow(X.spike)) temp.spike[j,] <- X.spike[j,] * rate.spike
	H <- (-1/2) * (temp.spike %*% t(X.spike)) - (temp %*% t(X)) - eta * diag(nrow(X))                             # Hessian
	 
	hist.theta <- matrix(NA, max.iter, length(theta.0)) # past values of the parameters
	hist.L <- rep(NA, max.iter)					                   # past values of the likelihood
	hist.theta[1,] <- theta.0
	hist.L[1] <- L
	theta <- theta.0
	 
	for (i in 2:max.iter){
	    # cat(i, "\n")
	    if((i>3) && (abs((L-hist.L[i-2])/hist.L[i-2]) < tol.theta)){
	        nbr.iter <- i
	        hist.theta[i,] <- theta
	        hist.L[i] <- L
	        	convergence <- "tolerance"
	        break    	
	    } else {
	        nbr.iter <- i

	        es <- tryCatch(solve(H), error=function(e) e)
	        if (!is.matrix(es)) {
				warning("the system is singlular - we stop fitting the threshold")
		        nbr.iter <- i
		        hist.theta[i,] <- theta
		        convergence <- "singular"
		        hist.L[i] <- L
        		break
	        }

	        theta <- theta - t(((solve(H)) %*% G))    # Fisher Scoring
			lrate.spike <- theta %*% X.spike           	# log rate at spike
			lrate <- theta %*% X	                        		# log rate at silence
			rate.spike <- exp(lrate.spike) 			        # rate at silence
			rate <- exp(lrate) 			                  		# rate at silence
			
			L <- sum(lrate.spike) -sum(rate.spike)/2 - sum(rate) - (sqrt(eta) * theta.0) %*% t(sqrt(eta) * theta.0)     		# likelihood - see Dayan and Abbott E. 1.37
			G <- X.spike %*% t(1-rate.spike/2) - X %*% t(rate) - eta * t(theta.0)       # gradient       # gradient
			temp <- matrix(NA, nrow(X), ncol(X))
			for (j in 1:nrow(X)) temp[j,] <- X[j,] * rate
			temp.spike <- matrix(NA, nrow(X.spike), ncol(X.spike))
			for (j in 1:nrow(X.spike)) temp.spike[j,] <- X.spike[j,] * rate.spike
			H <- (-1/2) * (temp.spike %*% t(X.spike)) - (temp %*% t(X)) - eta * diag(nrow(X)) 	                              # Hessian
	        hist.theta[i,] <- theta
	        hist.L[i] <- L
	 	}
 	}
 	 if (i == max.iter) convergence <- "maxit" 
 	hist.theta <- hist.theta[!is.na(hist.L),]
 	hist.L <- hist.L[!is.na(hist.L)]
 	out <- list(theta=theta, iter=nbr.iter, hist.L=hist.L, hist.theta=hist.theta, convergence=convergence)  
 	out
}

