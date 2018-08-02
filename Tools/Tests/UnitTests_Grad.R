
##############################################
## A function to test the subthreshold gradients of GLM models
## First, we test the grad.hGLM function 
## Second the err.hGLM() function - it should do the same as the sim.hGLM()
## Third, the grad.err.hGLM function - it should do the same as the grad.hGLM() 
## Fourth, we test the derivatives using the pakage numderiv

###############################################
# - test the w.ahp
# - test regularisation
###############################################

source('../Grad_hGLM.R', chdir=T)
source('../../Utils/Spike_Basis.R', chdir = TRUE)
# source('./UnitTests_Init.R', chdir = TRUE)

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
raster(X)
par1 <- init.hGLM(X, 1, pars, rand.seed=1, lin.soma=T, rescale.soma=T, nSD=1)
par2 <- init.hGLM(X, 1, pars, rand.seed=2, lin.soma=T, rescale.soma=T, nSD=1)
lines(par1$v+5, col=2); abline(h=5)
lines(par2$v+5, col=3)


# dt <- 1; pars <- par1; vv <- par1$v; logpars <- T; regpars <- NULL; double <- F; scale <- 2.8; verbose <- F; calc.minis <- F; X.ahp <- NULL
g12 <- unlist(grad.hGLM(X, 1, par1, par2$v))
g11 <- unlist(grad.hGLM(X, 1, par1, par1$v))

plot(g12)
points(g11, col=grey(0.75), pch=16)
points(unlist(grad.hGLM(X, 1, par2, par2$v)), col=rgb(1, 0.6, 0.6))
points(unlist(grad.hGLM(X, 1, par2, par1$v)), col=2)

if (max(abs(g11)) > 1e-10) warning('gradient should be zero if response is perfect! - 1L')

# 2. 1N - nonlinear - output variance should be 1
raster(X)
par1 <- init.hGLM(X, 1, pars, rand.seed=1, lin.soma=F, rescale.soma=T, nSD=1/2)
par2 <- init.hGLM(X, 1, pars, rand.seed=2, lin.soma=F, rescale.soma=T, nSD=1/2)
lines(par1$v+5, col=2); abline(h=5)
lines(par2$v+5, col=3)


g12 <- unlist(grad.hGLM(X, 1, par1, par2$v))
g11 <- unlist(grad.hGLM(X, 1, par1, par1$v))

plot(g12)
points(g11, col=grey(0.75), pch=16)
points(unlist(grad.hGLM(X, 1, par2, par2$v)), col=rgb(1, 0.6, 0.6))
points(unlist(grad.hGLM(X, 1, par2, par1$v)), col=2)

if (max(abs(g11)) > 1e-10) warning('gradient should be zero if response is perfect! - 1L')


# 3. 3N - random initialization - try a couple of different
raster(X)
par1 <- init.hGLM(X, 1, pars3, rand.seed=1, lin.soma=F, rescale.soma=T, nSD=1/2)
par2 <- init.hGLM(X, 1, pars3, rand.seed=2, lin.soma=F, rescale.soma=T, nSD=1/2)
lines(par1$v+5, col=2); abline(h=5)
lines(par2$v+5, col=3)

g12 <- unlist(grad.hGLM(X, 1, par1, par2$v))
g11 <- unlist(grad.hGLM(X, 1, par1, par1$v))

plot(g12)
points(g11, col=grey(0.75), pch=16)
points(unlist(grad.hGLM(X, 1, par2, par2$v)), col=rgb(1, 0.6, 0.6))
points(unlist(grad.hGLM(X, 1, par2, par1$v)), col=2)

if (max(abs(g11)) > 1e-10) warning('gradient should be zero if response is perfect! - 1L')


###########################################################################
###########################################################################

