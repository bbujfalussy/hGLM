
##############################################
## A function to test the subthreshold simulation and initialisation of GLM models

source('../Sim_hGLM_sub.R', chdir=T)
source('../../Utils/Spike_Basis.R', chdir = TRUE)
library(viridis)

############################################################
# setting up the model architectures

L <- 2000
X <- matrix(rbinom(20*L, 1, 0.005), 20, L)
Jc <- c(0)
Wc.e <- list(seq(1,10))
Wc.i <- list(seq(11,20))
pars <- list(Jc=Jc, Wc.e=Wc.e, Wc.i=Wc.i)
 
Jc2 <- c(0,1,1)
Wc.e2 <- list(seq(1,5), NULL, seq(6,10))
Wc.i2 <- list(seq(11,15), seq(16,20), NULL)
pars2 <- list(Jc=Jc2, Wc.e=Wc.e2, Wc.i=Wc.i2)

Jc3 <- c(0,1,1, 2,2, 3,3)
Wc.e3 <- list(NULL, 1,2, 3:4, 5:6, 7:8, 9:10)
Wc.i3 <- list(seq(11,15), 16, 17, 18, NULL, 19, 20)
pars3 <- list(Jc=Jc3, Wc.e=Wc.e3, Wc.i=Wc.i3)


# 1. 1L - linear - output variance should be 1
{
	raster(X, dt=1)
	newpars.1L <- init.hGLM(X, 1, pars, rand.seed=1, lin.soma=T, rescale.soma=T, nSD=1, mean.v=10, sd.v=1, alphasyn=F)
	lines(newpars.1L$v, col=2); abline(h=5)
	
	# pars <- newpars.1L; vv <- NULL; logpars <- T; regpars <- NULL; double <- F; scale <- 2.8; verbose <- F; calc.minis <- F; X.ahp <- NULL
	
	v <- sim.hGLM(X, 1, newpars.1L, alphasyn=F)
	lines(v, col=3, lty=2)
	if (abs(sd(newpars.1L$v) - 1) > 0.01) warning("1L sd is not 1!")
	if (max((v- newpars.1L$v)^2) > 1e-10) warning("sim and init is not doing the same!")
}

# 1. 1L - linear - output variance should be 1, no inhibitory cells
{
	raster(X, dt=1)
	pars <- list(Jc=Jc, Wc.e=Wc.e, Wc.i=list(NULL))
	
	newpars.1L <- init.hGLM(X, 1, pars, rand.seed=1, lin.soma=T, rescale.soma=T, nSD=1, mean.v=10, sd.v=1)
	lines(newpars.1L$v, col=2); abline(h=5)
	
	# pars <- newpars.1L; vv <- NULL; logpars <- T; regpars <- NULL; double <- F; scale <- 2.8; verbose <- F; calc.minis <- F; X.ahp <- NULL
	
	v <- sim.hGLM(X, 1, newpars.1L)
	lines(v+5, col=3, lty=2)
	if (abs(sd(newpars.1L$v) - 1) > 0.01) warning("1L sd is not 1!")
	if (max((v- newpars.1L$v)^2) > 1e-10) warning("sim and init is not doing the same!")
}



# 2. 1N - nonlinear - output variance should be 1
{
	for (nn in c(1/10, 1, 10)){
		newpars.1N <- init.hGLM(X, 1, pars, rand.seed=2, lin.soma=F, nSD=nn, rescale.soma=T, mean.v=10, sd.v=2)
		lines(newpars.1N$v, col=3); abline(h=5)
		if (abs(sd(newpars.1N$v) - 2) > 0.01) warning("1N sd is not 1!")
	}
}

# 3. 2N - random initialization - try a couple of different
{
	newpars.2N <- init.hGLM(X, 1, pars2, rand.seed=3, lin.soma=F, rescale.soma=T)
	lines(newpars.2N$v+5, col=4, lty=1); abline(h=5)
	if (abs(sd(newpars.2N$v) - 1) > 0.01) warning("2N sd is not 1!")
}

# 4. 3N - random initialization - try a couple of different
{
	raster(X)
	newpars.3N <- init.hGLM(X, 1, pars3, rand.seed=5, lin.soma=F, nSD=1)
	lines(newpars.3N$v+5, col=2, lty=1); abline(h=5)
	if (abs(sd(newpars.3N$v) - 1) > 0.01) warning("3N sd is not 1!")

	v <- sim.hGLM(X, 1, newpars.3N)
	lines(v+5, col=3, lty=2)
	if (max((v- newpars.3N$v)^2) > 1e-10) warning("sim and init is not doing the same!")
}

##############################################
# 5. make a 3L identical to an 1L
##############################################
{
	raster(X)
	pars <- list(Jc=Jc, Wc.e=Wc.e, Wc.i=Wc.i)
	newpars.1L <- init.hGLM(X, 1, pars, rand.seed=5, lin.soma=T, nSD=1)
	lines(newpars.1L$v+5, col=1); abline(h=5)
	
	newpars.1L$v0 <- 5
	newpars.1L$Ww.e[[1]] <- newpars.1L$Ww.e[[1]] *2
	v <- sim.hGLM(X, 1, newpars.1L)
	lines(v, col=2, lwd=3)
	
	plot(v, v, pch=16)
	initp.3L <- reshape.pars.syn2(Jc=Jc3, Wc.e=Wc.e3, Wc.i=Wc.i3, newpars.1L) # reshape the 1L parameters to match the 3L
	k <- 1
	for (nn in c(0.5, 1, 2, 4, 8, 16)){
		newpars.3L <- init.hGLM(X, 1, initp.3L, rand.seed=23, lin.soma=T, rescale.soma=F, nSD = nn) # initialize the 3L parameters to match the 1L
		# lines(newpars.3L$v, col=viridis(6)[k], lty=1); abline(h=5)
		points(v, newpars.3L$v, pch=16, cex=0.5, col=viridis(6)[k])
		k <- k + 1
	}
}

##############
## 3N to 1N
{
	raster(X)
	newpars.1N <- init.hGLM(X, 1, pars, rand.seed=6, lin.soma=F, alphasyn=F, double=F)
	lines(newpars.1N$v+5, col=1); abline(h=5)
	newpars.1N$v0 <- -5
	newpars.1N$Jw <- newpars.1N$Jw + 1/2
	v <- sim.hGLM(X, 1, newpars.1N, alphasyn=F)
	lines(v, col=2, lwd=3)
	
	
	
	plot(v, v, pch=16)
	initp.3N <- reshape.pars.syn2(Jc=Jc3, Wc.e=Wc.e3, Wc.i=Wc.i3, newpars.1N) # reshape the 1L parameters to match the 3L
	k <- 1
	for (nn in c(0.5, 1, 2, 4, 8, 16)){
		newpars.3N <- init.hGLM(X, 1, initp.3N, rand.seed=23, lin.soma=F, rescale.soma=F, nSD = nn, alphasyn=F) # initialize the 3L parameters to match the 1L
		# lines(newpars.3N$v+5, col=viridis(6)[k], lty=1); abline(h=5)
		points(v, newpars.3N$v, pch=16, cex=0.5, col=viridis(6)[k])
		k <- k + 1
	}
}

##############################################
# the same tests with double - synapses
##############################################
{
	# 6. 1L - linear - output variance should be 1
	raster(X)
	newpars.1L <- init.hGLM(X, 1, pars, rand.seed=1, lin.soma=T, rescale.soma=T, nSD=1, double=T)
	lines(newpars.1L$v+5, col=2); abline(h=5)
	v <- sim.hGLM(X, 1, newpars.1L, double=T)
	lines(v+5, col=3, lty=2)
	if (abs(sd(newpars.1L$v) - 1) > 0.01) warning("1L sd is not 1!")
	if (max((v- newpars.1L$v)^2) > 1e-10) warning("sim and init is not doing the same!")	
}

# 7. 1N - nonlinear - output variance should be 1
{
	k <- 1
	for (nn in c(1/10, 1, 10)){
		newpars.1N <- init.hGLM(X, 1, pars, rand.seed=2, lin.soma=F, nSD=nn, rescale.soma=T, double=T)
		lines(newpars.1N$v+5, col=rainbow(3)[k]); abline(h=5)
		if (abs(sd(newpars.1N$v) - 1) > 0.01) warning("1N sd is not 1!")
		k <- k + 1
	}
}

# 8. 2N - random initialization - try a couple of different
{
	newpars.2N <- init.hGLM(X, 1, pars2, rand.seed=3, lin.soma=F, rescale.soma=T, double=T)
	lines(newpars.2N$v+5, col=4, lty=1); abline(h=5)
	if (abs(sd(newpars.2N$v) - 1) > 0.01) warning("2N sd is not 1!")
}

# 9. 3N - random initialization - try a couple of different
{
	raster(X)
	newpars.3N <- init.hGLM(X, 1, pars3, rand.seed=5, lin.soma=F, nSD=1, double=T)
	lines(newpars.3N$v+5, col=2, lty=1); abline(h=5)
	if (abs(sd(newpars.3N$v) - 1) > 0.01) warning("3N sd is not 1!")
	
	v <- sim.hGLM(X, 1, newpars.3N, double=T)
	lines(v+5, col=3, lty=2)
	if (max((v- newpars.3N$v)^2) > 1e-10) warning("sim and init is not doing the same!")
}

##############################################
# 10. make a 3N identical to an 1L
##############################################
{
	raster(X)
	newpars.1L <- init.hGLM(X, 1, pars, rand.seed=5, lin.soma=T, nSD=1, double=T)
	lines(newpars.1L$v+5, col=1); abline(h=5)
	
	newpars.1L$v0 <- 5
	newpars.1L$Ww.e[[1]] <- newpars.1L$Ww.e[[1]] *2
	v <- sim.hGLM(X, 1, newpars.1L, double=T)
	lines(v, col=2, lwd=3)
	
	plot(v, v, pch=16)
	initp.3L <- reshape.pars.syn2(Jc=Jc3, Wc.e=Wc.e3, Wc.i=Wc.i3, newpars.1L) # reshape the 1L parameters to match the 3L
	k <- 1
	for (nn in c(0.5, 1, 2, 4, 8, 16)){
		newpars.3L <- init.hGLM(X, 1, initp.3L, rand.seed=23, lin.soma=T, rescale.soma=F, nSD = nn, double=T) # initialize the 3L parameters to match the 1L
		# lines(newpars.3L$v, col=viridis(6)[k], lty=1); abline(h=5)
		points(v, newpars.3L$v, pch=16, cex=0.5, col=viridis(6)[k])
		k <- k + 1
	}
}
	
##############
## 3N to 1N
{
	raster(X)
	newpars.1N <- init.hGLM(X, 1, pars, rand.seed=6, lin.soma=F)
	lines(newpars.1N$v+5, col=1); abline(h=5)
	newpars.1N$v0 <- -5
	newpars.1N$Jw <- newpars.1N$Jw + 1/2
	v <- sim.hGLM(X, 1, newpars.1N)
	lines(v, col=2, lwd=3)
	
	plot(v, v, pch=16)
	initp.3N <- reshape.pars.syn2(Jc=Jc3, Wc.e=Wc.e3, Wc.i=Wc.i3, newpars.1N) # reshape the 1L parameters to match the 3L
	k <- 1
	for (nn in c(0.5, 1, 2, 4, 8, 16)){
		newpars.3N <- init.hGLM(X, 1, initp.3N, rand.seed=23, lin.soma=F, rescale.soma=F, nSD = nn) # initialize the 3L parameters to match the 1L
		# lines(newpars.3N$v+5, col=viridis(6)[k], lty=1); abline(h=5)
		points(v, newpars.3N$v, pch=16, cex=0.5, col=viridis(6)[k])
		k <- k + 1
	}
}

############################################
# Check the after-spike currents...

n.basis <- 8
basis.a <- gen.basis(Tmax=100, dt=1, graphics=T, n.basis=n.basis, dur=1)	
spt <- list(st=seq(100, 1900, by=100), delay=0)
X.ahp <- spike.basis(spt, L, basis.a, dt=1)
matplot(t(X.ahp), t="l", lty=1, col=rainbow(8))

# 6. 1L - linear - output variance should be 1
raster(X)
newpars.1L <- init.hGLM(X, 1, pars, rand.seed=1, lin.soma=T, rescale.soma=T, nSD=1, double=T)
lines(newpars.1L$v+5, col=2); abline(h=5)

newpars.1L <- init.hGLM(X, 1, pars, rand.seed=1, lin.soma=T, rescale.soma=T, nSD=1, double=T, X.ahp=X.ahp)
lines(newpars.1L$v+5, col=2); abline(h=5)



v <- sim.hGLM(X, 1, newpars.1L, double=T, X.ahp=X.ahp)
lines(v+5, col=3, lty=2)
if (abs(sd(newpars.1L$v) - 1) > 0.01) warning("1L sd is not 1!")
if (max((v- newpars.1L$v)^2) > 1e-10) warning("sim and init is not doing the same!")


############################################
# Check regularization - just check the error
############################################
## 1L
{
	refpars.1L <- init.hGLM(X, 1, pars, rand.seed=13, lin.soma=T, rescale.soma=T, nSD=1)
	
	X.minis <- matrix(0, 20, L)
	for (i in 1:20) X.minis[i, i*90] <- 1
	v.minis <- sim.hGLM(X.minis, 1, refpars.1L, calc.minis=T)
	raster(X.minis, dt=1)
	lines(v.minis$v+5, t="l")
	points(seq(90, by=90, length=10), v.minis$a.e+5 + v.minis$v[1], col=2, pch="_")
	points(seq(990, by=90, length=10), v.minis$a.i+5 + v.minis$v[1], col=4, pch="_")
	
	regpars <- list(Ww.e1=v.minis$a.e, Ww.i=v.minis$a.i, logTau.e=unlist(refpars.1L$Tau.e), logTau.i=unlist(refpars.1L$Tau.i), alpha.Ww=10, alpha.Tau=20)
	err <- sim.hGLM(X, 1, refpars.1L, vv=refpars.1L$v)
	if (err > 1e-10) warning('error is too large, sim is not doing the same as init')
	
	err2 <- sim.hGLM(X, 1, refpars.1L, vv=refpars.1L$v, regpars=regpars, verbose=T)
	if ((err2-err) > 1e-20) warning('regularisation error - should be 0, but is not')
	
	newpars.1L <- init.hGLM(X, 1, pars, rand.seed=5, lin.soma=T, rescale.soma=T, nSD=1)
	err <- sim.hGLM(X, 1, newpars.1L, vv=refpars.1L$v, verbose=T)
	err2 <- sim.hGLM(X, 1, newpars.1L, vv=refpars.1L$v, regpars=regpars, verbose=T)
	if ((err2-err) < 1) warning('regularisation error - should be > 0, but is 0')
}

############################################
# double - alpha regularisation
{
	refpars.1L <- init.hGLM(X, 1, pars, rand.seed=13, lin.soma=T, rescale.soma=T, nSD=1, double=T)
	
	X.minis <- matrix(0, 20, L)
	for (i in 1:20) X.minis[i, i*90] <- 1
	v.minis <- sim.hGLM(X.minis, 1, refpars.1L, calc.minis=T, double=T)
	raster(X.minis, dt=1)
	lines(v.minis$v+5, t="l")
	points(seq(90, by=90, length=10), v.minis$a.e2 + v.minis$a.e+5 + v.minis$v[1], col=2, pch="_")
	points(seq(990, by=90, length=10), v.minis$a.i+5 + v.minis$v[1], col=4, pch="_")
	
	regpars <- list(Ww.e1=v.minis$a.e, Ww.e2=v.minis$a.e2, Ww.i=v.minis$a.i, logTau.e=unlist(refpars.1L$Tau.e), logTau.i=unlist(refpars.1L$Tau.i), alpha.Ww=10, alpha.Tau=20)
	err <- sim.hGLM(X, 1, refpars.1L, vv=refpars.1L$v, double=T)
	if (err > 1e-10) warning('error is too large, sim is not doing the same as init')
	
	err2 <- sim.hGLM(X, 1, refpars.1L, vv=refpars.1L$v, regpars=regpars, verbose=T, double=T)
	if ((err2-err) > 1e-20) warning('regularisation error - should be 0, but is not')
	
	newpars.1L <- init.hGLM(X, 1, pars, rand.seed=5, lin.soma=T, rescale.soma=T, nSD=1, double=T)
	err <- sim.hGLM(X, 1, newpars.1L, vv=refpars.1L$v, verbose=T, double=T)
	err2 <- sim.hGLM(X, 1, newpars.1L, vv=refpars.1L$v, regpars=regpars, verbose=T, double=T)
	if ((err2-err) < 1) warning('regularisation error - should be > 0, but is 0')
	
	raster(X, dt=1)
	lines(newpars.1L$v+5, col=1); abline(h=5)
	lines(refpars.1L$v+5, col=2)
}

############################################
# double - alpha regularisation - see what happens with the regularisation if we match a 3N model with a 1N
{
	refpars.1N <- init.hGLM(X, 1, pars, rand.seed=13, lin.soma=F, rescale.soma=T, nSD=1, double=T)
	raster(X, dt=1)
	lines(refpars.1N$v+5, col=2)
	
	
	X.minis <- matrix(0, 20, 2*L)
	for (i in 1:20) X.minis[i, i*180] <- 1
	v.minis <- sim.hGLM(X.minis, 1, refpars.1N, calc.minis=T, double=T)
	raster(X.minis, dt=1)
	lines(v.minis$v+5, t="l")
	points(seq(180, by=180, length=10), v.minis$a.e2 + v.minis$a.e+5 + v.minis$v[1], col=2, pch="_")
	points(seq(1980, by=180, length=10), v.minis$a.i+5 + v.minis$v[1], col=4, pch="_")
	
	regpars <- list(Ww.e1=v.minis$a.e, Ww.e2=v.minis$a.e2, Ww.i=v.minis$a.i, logTau.e=unlist(refpars.1N$Tau.e), logTau.i=unlist(refpars.1N$Tau.i), alpha.Ww=10, alpha.Tau=20)
	err <- sim.hGLM(X, 1, refpars.1N, vv=refpars.1N$v, double=T)
	if (err > 1e-10) warning('error is too large, sim is not doing the same as init')
	
	err2 <- sim.hGLM(X, 1, refpars.1N, vv=refpars.1N$v, regpars=regpars, verbose=T, double=T)
	if ((err2-err) > 1e-20) warning('regularisation error - should be 0, but is not')
	
	newpars.3N <- init.hGLM(X, 1, pars3, rand.seed=5, lin.soma=F, rescale.soma=T, nSD=1, double=T)
	err <- sim.hGLM(X, 1, newpars.3N, vv=refpars.1N$v, verbose=T, double=T)
	err2 <- sim.hGLM(X, 1, newpars.3N, vv=refpars.1N$v, regpars=regpars, verbose=T, double=T)
	if ((err2-err) < 1) warning('regularisation error - should be > 0, but is 0')
	
	raster(X, dt=1)
	lines(newpars.3N$v+5, col=1); abline(h=5)
	lines(refpars.1N$v+5, col=2)
	
	initp.3N <- reshape.pars.syn2(Jc=Jc3, Wc.e=Wc.e3, Wc.i=Wc.i3, refpars.1N) # reshape the 1L parameters to match the 3L
	newpars.3N <- init.hGLM(X, 1, initp.3N, rand.seed=23, lin.soma=F, rescale.soma=F, nSD = 8, double=T) # initialize the 3L parameters to match the 1L
	err <- sim.hGLM(X, 1, newpars.3N, vv=refpars.1N$v, verbose=T, double=T)
	if (err > 1) warning('regularisation error, should be near zero!')
	
	err2 <- sim.hGLM(X, 1, newpars.3N, vv=refpars.1N$v, regpars=regpars, verbose=T, double=T)
	if ((err2-err) > 1) warning('regularisation error - should be < 1, but it is larger')
	
	
	raster(X, dt=1)
	lines(newpars.3N$v+5, col=1); abline(h=5)
	lines(refpars.1N$v+5, col=2, lty=2)
}
