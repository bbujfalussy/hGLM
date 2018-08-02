## this function plots the mean of a set of samples and a shaded interval around them at mean ± nsd * sd.
polyplot <- function(t=NULL, samples, col.bg=grey(0.8), col=1, xlab=NULL, ylab=NULL, ylim=NULL, add=F, axes=T, tit=NULL, nsd=1, sem=F){


	if (is.null(t)) t <- seq(1, ncol(samples))
	
	if (sem) N <- ncol(samples) else  N <- 1	
	m.mat <- colMeans(samples, na.rm=T)
	sd.mat <- nsd * apply(samples, 2, sd, na.rm=T) / sqrt(N)

	tt <- c(t, rev(t))
	ss <- c(m.mat+sd.mat, rev(m.mat-sd.mat))
	if (is.null(ylim)) ylim <- range(ss)
	if (add==F)	plot(t, m.mat, ylim=ylim, t="n", axes=F, xlab=xlab, ylab=ylab, main=tit)
	polygon(tt,ss, col=col.bg, border=grey(1))
	if (axes) {axis(1); axis(2, las=2)}
	lines(t, m.mat, lwd=2, col=col)
	lines(t, m.mat+sd.mat, lwd=1, lty=2, col=col)
	lines(t, m.mat-sd.mat, lwd=1, lty=2, col=col)

}

## plots two functions and colors the area between them - similar to polygon, but a bit nicer result though less general
polyplot2 <- function(t, sig1, sig2, col1=1, col2=2, col.bg=grey(0.8), xlab="", ylab="", main="", axes=F, ylim=NULL){
# col1=1; col2=2; col.bg=grey(0.8); xlab=""; ylab=""; main=""; axes=F; ylim=NULL;
	if (length(t) > 5000) {
		warning("data size is too big - only the first 5000 points are shown")
		sig1 <- sig1[1:5000]
		sig2 <- sig2[1:5000]
		t <- t[1:5000]
	}

	if (is.null(ylim)) ylim <- range(c(sig1, sig2))
	tt <- c(t, rev(t))
	plot (tt, c(sig1, rev(sig2)), xlab=xlab, ylab=ylab, main=main, axes=F, t="n", ylim=ylim)
	polygon(tt, c(sig1, rev(sig2)), col=col.bg)
	lines(t, sig1, col=col1)
	lines(t, sig2, col=col2)
	if (axes) {axis(1); axis(2, las=2)}	
}

## similar to matplot but handles confidence intervals
## plots a matrix of means and confidence intervals
matplot.CI <- function(x=NULL, ms, sdu=NULL, sdl=sdu, cols=1, bgs=1, xlab=NULL, ylab=NULL, ylim=NULL, xlim=NULL, add=F, axes=T, tit=NULL, nsd=1, log="", pch=21, cex=1){
	
	library(gplots)
	if (is.matrix(ms)) n <- nrow(ms) else n <- 1
	if (is.matrix(ms)) m <- ncol(ms) else n <- length(ms)
	if (is.null(x)) x <- seq(1, m)
	if (is.null(sdu)) {
		sdu <- matrix(0, nrow(ms), ncol(ms))
		sdl=sdu
	}
	if (is.null(ylim)) ylim = range(c(ms+sdu, ms-sdl))
	if (is.null(xlim)) xlim = range(x)
	cols <- cols[seq(0,n-1) %% length(cols) + 1]
	bgs <- bgs[seq(0,n-1) %% length(bgs) + 1]
	if (is.vector(x)) x <- matrix(rep(x, n), n, byrow=T)
	
	if (n==1) {
		plotCI(x[1,], ms, sdu, sdl, pch=pch[1], pt.bg=bgs[1], col=cols[1], xlab=xlab, ylab=ylab, add=add, axes=F, log=log, gap=0, cex=cex, ylim=ylim, main=tit[1])
		lines(x,ms,t="o", col=cols[1], pch=pch, bg=bgs[1], cex=cex)
	} else {
		if (length(cex)<n) cex <- rep(cex[1], n)
		if (length(pch)<n) pch <- rep(pch[1], n)
		plotCI(x[1,], ms[1,], sdu[1,], sdl[1,], pch=pch[1], col=cols[1], xlab=xlab, ylab=ylab, add=add, axes=F, ylim=ylim, log=log, xlim=xlim, gap=0, cex=cex[1], main=tit[1])
		lines(x[1,],ms[1,], col=cols[1], t="o", pch=pch[1], bg=bgs[1], cex=cex[1])
		for(i in 2:n){
			plotCI(x[i,], ms[i,], sdu[i,], sdl[i,], pch=pch[i], pt.bg=bgs[i], col=cols[i], add=T, gap=0, cex=cex[i], main=tit[i])
			lines(x[i,],ms[i,], col=cols[i], t="o", pch=pch[i], bg=bgs[i], cex=cex[i])
		}
	}
	if (axes) {axis(1); axis(2, las=2)}
}

## similar to matplot but handles confidence intervals
## plots a matrix of means and confidence intervals + the original point in the background
points.CI <- function(samples, x=NULL, cols=1, bgs=1, xlab=NULL, ylab=NULL, axes=T, tit=NULL, nsd=1, log="", t.pt="o", pch=21, ...){
	
	library(gplots)
	ms <- colMeans(samples)
	sds <- apply(samples, 2, sd)

	if (is.null(x)) x <- seq(1, ncol(samples))
	matplot(x, t(samples), pch=23, col=grey(0.8), bg=grey(0.8), t=t.pt, lty=1, cex=0.7, xlab=xlab, ylab=ylab, axes=F, log=log, ...)
	plotCI(x, ms, sds, pt.bg=bgs, col=cols, add=T, gap=0, cex=2, pch=pch)
	if (axes) {axis(1); axis(2, las=2)}
}

# samples <- matrix(runif(120), ncol=6)
# points.CI(samples)

######################################################################
## from a binary matrix of spikes plots a rastergram

raster <- function(sp, dt=0.5, col=1,ymax=NULL, xmax=NULL, assemblies=F, x=NULL, u.r=NULL, add=F){
	# if no states or emissions than we can not mark assemblies
	# col <- 1; ymax <- NULL; xmax <- NULL; assemblies <- F; x <- NULL; u.r <- NULL; add <- F

	if (assemblies) {
		if (is.null(x) + is.null(u.r)){
			warning("missing state or emission matrix")
			assemblies <- F
		}
	}

	N <- nrow(sp)
	L <- ncol(sp)
	if (max(sp)>1) sp <- !!sp
	t <- seq(0, length=L, by=dt)
	if (is.null(xmax)) Tmax <- ceiling(max(t)/100) * 100 else Tmax <- xmax
	xlim=c(0, Tmax)
	if (is.null(ymax)) ylim=c(0, N+1) else ylim=c(0, ymax)
	if (!add) plot(t[!!sp[1,]], rep(1, sum(sp[1,])), t="n", ylim=ylim, xlim=xlim, axes=F, xlab="time (s)", ylab="cells")
	# mark assemblies if required
	if (assemblies) mark.assemblies(x, u.r, t)
	# show the spikes
	points(t[!!sp[1,]], rep(1, sum(sp[1,])), pch="|", col=col)
	if (!add) axis(1, seq(0, Tmax, length=6), seq(0, Tmax/1000, length=6)); axis(2, las=2)
	if (N>1) for (i in 2:N) 	points(t[!!sp[i,]], rep(i, sum(sp[i,])), pch="|")
}


######################################################################
## I want a function that draws nice background around each cell assembly activated
mark.assemblies <- function(x, u.r, t){
	n.cells <- ncol(u.r)
	n.st <- nrow(u.r)
	col.box <- hsv(seq(0,1,length=(n.st)), 0.15, 1)


	# we will use this to select the position of the boxes
	ur <- (u.r==max(u.r))

	## generate a matrix with the corners of the boxes
	Y.box <- matrix(c(0.5, 0.5, 1.5, 1.5), 1, 4)
	if (n.cells>1){
		for (i in 2:n.cells){
			Y.box <- rbind(Y.box, Y.box[1,]+i-1)
		}
	}
	
	## which is the quiescent state?
	st.0 <- which(rowSums(ur)==0)

	## the activations of the assemblies
	L <- length(x)
	i.up <- which((x[1:(L-1)]==st.0) & (x[-1]!=st.0))
	i.do <- which((x[1:(L-1)]!=st.0) & (x[-1]==st.0))
	if (x[1] != st.0) i.up <- c(1, i.up)
	if (x[L] != st.0) i.do <- c(i.do, L)

	if (length(i.do) < length(i.up)) stop("number of up and down states is different!")

	## draw a box around each assembly
	for (i in 1:length(i.up)){
		## which time?
		X.box <- c(i.up[i], i.do[i], i.do[i], i.up[i])
		## how many neurons?
		st <- x[i.up[i]+1]
		nn <- sum(ur[st,])
		x.box <- t[rep(X.box, nn)]
		
		## which many neurons?
		i.cells <- ur[st,]
		y.box <- as.vector(t(Y.box[i.cells,]))

		## plot the boxes!
		polygon(x.box, y.box, col=col.box[st], border=NA)	
	}	
}


