
spike.reliability <- function(V.ref, V.opt, dt, graphics=F, tmin=4000, tmax=6000, dt.spike=4){

	if (is.matrix(V.ref)) {
		if (ncol(V.ref) < nrow(V.ref)) {
			V.ref <- t(V.ref)
			warning("if V.ref is matrix: columns are time points, rows are different traces")	
		}
	} else {
		V.ref <- matrix(V.ref, nrow=1)
	}
	N.ref <- nrow(V.ref)	

	if (is.matrix(V.opt)) {
		if (ncol(V.opt) < nrow(V.opt)) {
			V.opt <- t(V.opt)
			warning("if V.opt is matrix: columns are time points, rows are different traces")	
		}
	} else {
		V.opt <- matrix(V.opt, nrow=1)
	}
	N.opt <- nrow(V.opt)	

	
	Lmax <- ncol(V.ref)
	if (Lmax != ncol(V.opt)) stop('V.ref and V.opt have different duration')
	Tmax <- Lmax * dt # ms
	time <- seq(dt, Tmax, by=dt)

	###############################################
	## Test reliability of spiking
	## gathering spiketimes in the reference data
	sp.ref <- matrix(0, Lmax, N.ref)
	spt.ref <- list()

	for (i in 1:N.ref){
		spt <- extract.spiketimes(V.ref[i,], sampling.freq=1000/dt, limite=-30, graphics=F) 	# Hz 
		spt.ref[[i]] <- spt$st *dt - spt$delay # the threshold corssings
		spt.ref[[i]] <- spt.ref[[i]][spt.ref[[i]] > 0 ] # only at positive times..
		sp.ref[round(spt.ref[[i]]/dt),i] <- 1
	}
	
	## We can test spiking reliability ONLY if the stimuli are the same
	sp.opt <- matrix(0, Lmax, N.opt)
	spt.opt <- list()
	
	for (i in 1:N.opt){
		spt <- extract.spiketimes(V.opt[i,], sampling.freq=1000/dt, limite=-30, graphics=F) 	# Hz 
		spt.opt[[i]] <- spt$st *dt - spt$delay # the threshold corssings
		spt.opt[[i]] <- spt.opt[[i]][spt.opt[[i]] > 0 ] # only at positive times..
		sp.opt[round(spt.opt[[i]]/dt),i] <- 1
	}
	cat("\n")
	
	if (graphics > 0){
		plot(c(tmin, tmax), c(-N.opt, N.ref), axes=F, t="n", xlab="time (ms)", ylab="spikes")
		ii <- (spt.ref[[1]] > tmin) & (spt.ref[[1]] < tmax)
		if (N.ref > 1) {
			for (i in 1:N.ref) {
				ii <- (spt.ref[[i]] > tmin) & (spt.ref[[i]] < tmax)
				points(spt.ref[[i]][ii], rep(i, length(spt.ref[[i]][ii])), pch=21, col=1, bg=1, cex=0.3) 
			}
		} else abline(v=spt.ref[[1]][ii], col=1)
		for (i in 1:N.opt){
			ii <- (spt.opt[[i]] > tmin) & (spt.opt[[i]] < tmax)
			points(spt.opt[[i]][ii], rep(-i, length(spt.opt[[i]][ii])), pch=21, col=3, bg=3, cex=0.3)
		}
		axis(1, seq(tmin, tmax, by=500), seq(0, tmax- tmin, by=500))
	}


	msp.ref <- rowMeans(sp.ref)
	msp.opt <- rowMeans(sp.opt)
	
	# N.init <- 200/dt	# initial transient - 200ms, the cell is baheving in a different way...
	# msp.ref <- msp.ref[-(1:N.init)]
	# msp.opt <- msp.opt[-(1:N.init)]
	
	rel.sp <- 2 * inprod.th(msp.ref, msp.opt, dt.spike/dt) / (inprod.th(msp.ref, msp.ref, dt.spike/dt) + inprod.th(msp.opt, msp.opt, dt.spike/dt))
	print(paste("reliability of spiking (in ", dt.spike, "ms window): ", round(rel.sp,4)))
	
	rel.sp
}


###########################################################

inprod.th <- function(x, y, delta){
   # Compute inner product th with replacement according to
   # <X,Y> = int_{0}^{T} int_{-inf}^{inf} int_{-inf}^{inf} ...
   # K_{Delta}(s,s')X(t-s)Y(t-s')dsds'dt
   # with K_{Delta}(s,s') = the coincidence detector windows (i.e a rect)

   # Here we have:
   # K_{Delta}(s,s') = dirac(s)h_{r}(s';Delta)
   # h_{r}(s;Delta) = Heaviside(s-Delta)Heaviside(Delta-s)
   
   # X and Y are N-dimensionnal vectors
 
	filt <- rep(1, 2*delta + 1)
 
	x1 <- filter(c(rep(0, delta), x, rep(0, delta+1)), filt) 
	x1 <- x1[(delta+1):(length(x)+delta)]
	out <- sum(x1 * y)
 	out
}

# x <- rep(0, 100)
# x[c(34, 38, 56, 87)] <- 1
# inprod.th(x, x, 4)

# y <- rep(0, 100)
# y[c(36, 36, 56, 87)] <- 1
# inprod.th(x, y, 4)
