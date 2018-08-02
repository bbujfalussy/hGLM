
################################################################

sigm <- function(x, a=1, b=1, c=0, d=0){
	# a / (1 + exp(-b(x-c)))
	#a - amplitude: limit in inf; b - slope; c - threshold
	# d: the derivative of the sigmoid function at x. d=0: the function; d=1: first derivative; d=2: second derivative
	# derivative with respect to the parameters b or c
	bc <- exp(-b*(x-c))
	if (is.numeric(d)){		
		if (d>2) return ("d must be smaller than 2!")
		if (d==0) {
			sx <- a / (1+bc)
			} else {
				if (d==1) {
					sx <- (a*b*bc)/(1+bc)^2
				} else {
					sx <- (a*b^2*bc*(bc-1))/(1+bc)^3
				}
				sx[(1+bc == bc)] <- 0
			}
		} else {
			if (d == "a"){
				sx <- 1 / (1+bc)
			}
			if (d == "b"){
				sx <- a * (x-c) * bc / (1+bc)^2	
			}
			if (d == "c"){
				sx <- -a * b * bc / (1+bc)^2	
			}
			if (!(d %in% c("a", "b", "c"))) return("invalid d value")
		}
	sx
}

# x <- seq(-10,10, le=201)
# plot (x, sigm(x, a=1, b=5, c=3), t="l")
# plot (x, sigm(x, a=1, b=5, c=3, d="b"), t="l")
# plot (x, sigm(x, a=1, b=5, c=3, d=2), t="l")
# plot (x, sigm(x, a=1, b=5, c=3, d="c"), t="l")
# plot (x, sigm(x, a=1, b=500, c=5, d=1), t="l")



#######################################################
## inverz of the sigmoid function
isigm <- function(y, a=1, b=1, c=0){
	# a / (1 + exp(-b(x-c)))
	#a - amplitude: limit in inf; b - slope; c - threshold
	if (min(y)<0) stop("y must be greater than 0")
	if (max(y)>a) stop("y must be smaller than a")
 	x <- c + 1/b * log(y/(a-y))
	x	
}

# a <- 2; th<- 3; sl<-1/10
# y <- sigm(seq(1,10), a, sl, th)
# isigm(y, a, sl, th)
# sigm(Inf)
# isigm(1)
# isigm(0)


# # ################################################################
# ## Calculates the error between a reference signal (args$signal) and the sigmoid function of the input (args$fspikes - filtered spikes). The nonlinearity has 4 parameters
sigm4.err <- function(pars, args){
	## calculates the error between the signal and a sigmoidal function of the estimated \sum_i u_i
	# pars: parameters of the nonlinear function (see: sigm()).
	# args is a list containing the signal, fspikes and err
	fspikes <- args$fspikes
	offset <- pars[1]; ampl <- pars[2]; slope <- pars[3]; threshold <- pars[4]
	nfspikes <- sigm(fspikes, ampl, slope, threshold) + offset
	if (args$err == T){
		signal <- args$signal
		if (is.vector(signal)){
			if (is.matrix(nfspikes)) nfspikes <- colMeans(nfspikes)
			err <- mean( (signal - nfspikes)^2 )
			out <- err / var(signal)
			} else {
				signal <- colMeans(signal)
				nfs <- colMeans(nfspikes)
				err <- mean( (signal - nfs)^2 )
				out <- err / var(signal)							
			}
	} else {
		out <- nfspikes
	}
	out
}
