
##############################################
## A function to test the subthreshold simulation and initialisation of GLM models

source('../Err_hGLM.R', chdir=T)
source('../../Utils/Spike_Basis.R', chdir = TRUE)
library('numDeriv')
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

graphics <- T
############################################
## we prepare a target or reference output - refpars1L - and measure the error and the gradients by comparing the predistions with the reference
alphasyn <- F
refpars <- init.hGLM(X, 1, pars, rand.seed=17, lin.soma=T, rescale.soma=T, nSD=1, alphasyn=alphasyn)
v.ref <- refpars$v; X.ahp <- NULL; regpars <- NULL; double <- F
optimvec <- NULL

############################################
# 1. 1L
############################################
## newpars is a random initialisation
newpars <- init.hGLM(X, 1, pars, rand.seed=1, lin.soma=T, rescale.soma=T, nSD=1, alphasyn=alphasyn)
source('./test_numgrad.R')

optimvec=c('Tau.i', 'Th')
source('./test_numgrad.R')
optimvec <- NULL

############################################
# 2. 1N
############################################
newpars <- init.hGLM(X, 1, pars, rand.seed=2, lin.soma=F, rescale.soma=T, nSD=1, alphasyn=alphasyn)
source('./test_numgrad.R')

############################################
# 3. 2N
############################################
newpars <- init.hGLM(X, 1, pars2, rand.seed=3, lin.soma=F, rescale.soma=T, nSD=1, alphasyn=alphasyn)
source('./test_numgrad.R')

############################################
# 4. 3L
############################################
newpars <- init.hGLM(X, 1, pars3, rand.seed=3, lin.soma=T, rescale.soma=T, nSD=1, alphasyn=alphasyn)
source('./test_numgrad.R')

############################################
# 5. 3N
############################################
newpars <- init.hGLM(X, 1, pars3, rand.seed=3, lin.soma=F, rescale.soma=T, nSD=1, alphasyn=alphasyn)
source('./test_numgrad.R')

############################################
# 6. double synapses
############################################
double <- T
newpars <- init.hGLM(X, 1, pars3, rand.seed=3, lin.soma=F, rescale.soma=T, nSD=1, double=double, alphasyn=alphasyn)
source('./test_numgrad.R')


############################################
# 7. test the regularisation...
############################################
# 1L, simple synapses
double <- F
X.ahp <- NULL

# regularisation parameters
X.minis <- matrix(0, 20, L)
for (i in 1:20) X.minis[i, i*90] <- 1
v.minis <- sim.hGLM(X.minis, 1, refpars, calc.minis=T)
regpars <- list(Ww.e1=v.minis$a.e, Ww.i=v.minis$a.i, logTau.e=unlist(refpars$Tau.e), logTau.i=unlist(refpars$Tau.i), alpha.Ww=10, alpha.Tau=20)

newpars <- init.hGLM(X, 1, pars, rand.seed=5, lin.soma=F, rescale.soma=T, nSD=10, double=double, X.ahp=X.ahp)
source('./test_numgrad.R')

# 2N, simple synapses
newpars <- init.hGLM(X, 1, pars2, rand.seed=3, lin.soma=F, rescale.soma=T, nSD=1, double=double, X.ahp=X.ahp)
source('./test_numgrad.R')

# 3N, simple synapses
newpars <- init.hGLM(X, 1, pars3, rand.seed=3, lin.soma=F, rescale.soma=T, nSD=10, double=double, X.ahp=X.ahp)
source('./test_numgrad.R')

# 3N, double synapses, nonlinear soma
double <- T
refpars <- init.hGLM(X, 1, pars, rand.seed=17, lin.soma=F, rescale.soma=T, nSD=1, double=double)
v.ref <- refpars$v

X.minis <- matrix(0, 20, L)
for (i in 1:20) X.minis[i, i*90] <- 1
v.minis <- sim.hGLM(X.minis, 1, refpars, calc.minis=T, double=double)
regpars <- list(Ww.e1=v.minis$a.e, Ww.e2=v.minis$a.e2, Ww.i=v.minis$a.i, logTau.e=v.minis$logTau.e, logTau.i=v.minis$logTau.i, alpha.Ww=10, alpha.Tau=20)

newpars <- init.hGLM(X, 1, pars, rand.seed=3, lin.soma=F, rescale.soma=T, nSD=1, double=double, X.ahp=X.ahp)
source('./test_numgrad.R')

newpars <- init.hGLM(X, 1, pars3, rand.seed=3, lin.soma=F, rescale.soma=T, nSD=4, double=double, X.ahp=X.ahp)
source('./test_numgrad.R')


############################################
# 8. ahp
############################################
regpars <- NULL
n.basis <- 8
basis.a <- gen.basis(Tmax=100, dt=1, graphics=T, n.basis=n.basis, dur=1)	
spt <- list(st=seq(100, 1900, by=100), delay=0)
X.ahp <- spike.basis(spt, L, basis.a, dt=1)
matplot(t(X.ahp), t="l", lty=1, col=rainbow(8))

raster(X)
double <- F
refpars <- init.hGLM(X, 1, pars, rand.seed=17, lin.soma=F, rescale.soma=T, nSD=1, double=double, X.ahp=X.ahp)
v.ref <- refpars$v
lines(v.ref + 55)

newpars <- init.hGLM(X, 1, pars3, rand.seed=3, lin.soma=F, rescale.soma=T, nSD=1, double=double)
lines(newpars$v+55, col=2); abline(h=5)

newpars <- init.hGLM(X, 1, pars3, rand.seed=3, lin.soma=F, rescale.soma=T, nSD=1, double=double, X.ahp=X.ahp)
lines(newpars$v+55, col=4)

source('./test_numgrad.R')

optimvec <- c('Tau.i', 'W.ahp')
source('./test_numgrad.R')

optimvec <- NULL

############################################
# 9. test the regularisation and ahp, simple synapses
############################################

X.ahp <- spike.basis(spt, L, basis.a, dt=1)
double <- F
refpars <- init.hGLM(X, 1, pars, rand.seed=17, lin.soma=F, rescale.soma=T, nSD=1, double=double, X.ahp=X.ahp)
v.ref <- refpars$v
v.minis <- sim.hGLM(X.minis, 1, refpars, calc.minis=T)
regpars <- list(Ww.e1=v.minis$a.e, Ww.i=v.minis$a.i, logTau.e=unlist(refpars$Tau.e), logTau.i=unlist(refpars$Tau.i), alpha.Ww=10, alpha.Tau=20)

newpars <- init.hGLM(X, 1, pars3, rand.seed=3, lin.soma=F, rescale.soma=T, nSD=1, double=double, X.ahp=X.ahp)
source('./test_numgrad.R')

