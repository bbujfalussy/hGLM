###########################################
## merge inputs according to the different subunits + combine them in time
merge.inputs.timebins <- function(X, Wc.e, Wc.i, Tpoi){
	N.edend <- unlist(lapply(Wc.e, length))
	i.edend <- which(N.edend>0)
	N.e <- sum(N.edend>0)
	
	N.idend <- unlist(lapply(Wc.i, length))
	i.idend <- which(N.idend>0)
	N.i <- sum(N.idend>0)

	Lmax <- ncol(X)
	Lpoi <- Lmax / Tpoi
	Npoi <- N.e + N.i
	
	Xpoi <- matrix(NA, Npoi, Lpoi)
	for (i.e in 1:length(i.edend)){
		xrange <- Wc.e[[i.edend[i.e]]]
		for (i.t in 1:Lpoi){
			trange <- ((i.t-1) * Tpoi + 1):(i.t * Tpoi)
			Xpoi[i.e,i.t] <- sum(X[xrange, trange])
		}
	}
	for (i.i in 1:length(i.idend)){
		xrange <- Wc.i[[i.idend[i.i]]]
		i.ii <- i.i + i.e
		for (i.t in 1:Lpoi){
			trange <- ((i.t-1) * Tpoi + 1):(i.t * Tpoi)
			Xpoi[i.ii,i.t] <- sum(X[xrange, trange])
		}
	}
	Xpoi
}

#######################################
## merge inputs according to the different subunits
merge.inputs <- function(X, Wc.e, Wc.i){
	N.edend <- unlist(lapply(Wc.e, length))
	i.edend <- which(N.edend>0)
	N.e <- sum(N.edend>0)
	
	N.idend <- unlist(lapply(Wc.i, length))
	i.idend <- which(N.idend>0)
	N.i <- sum(N.idend>0)

	Lmax <- ncol(X)
	Npoi <- N.e + N.i
	
	Xpoi <- matrix(NA, Npoi, Lmax)
	for (i.e in 1:N.e){
		xrange <- Wc.e[[i.edend[i.e]]]
		Xpoi[i.e,] <- colSums(X[xrange, ])
		# print(xrange)
	}
	for (i.i in 1:N.i){
		xrange <- Wc.i[[i.idend[i.i]]]
		i.ii <- i.i + i.e
		Xpoi[i.ii,] <- colSums(X[xrange,])
		# print(xrange)
	}
	Xpoi
}

