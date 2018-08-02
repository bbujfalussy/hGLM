
##############################################
## A function to test the subthreshold simulation and initialisation of GLM models

source('../Sim_hGLM_sub.R', chdir=T)
source('../../Utils/Spike_Basis.R', chdir = TRUE)

############################################################
# setting up the model architectures

Ne <- 15
Ni <- 4
N <- Ne + Ni
L <- 200 * N

X <- matrix(0, N, L)
for ( i in 1:Ne) X[i,(i-1)*200 + 10] <- 1
X <- cbind(X, matrix(rbinom(N*L, 1, 0.005), N, L))

Jc <- c(0)
Wc.e <- list(seq(1,Ne))
Wc.i <- list(seq(Ne+1,N))
pars <- list(Jc=Jc, Wc.e=Wc.e, Wc.i=Wc.i)
  
Jc1b <- c(0)
Wc.e1b <- list(sample(1:Ne, Ne))
Wc.i1b <- list((Ne+1):N)
pars1b <- list(Jc=Jc1b, Wc.e=Wc.e1b, Wc.i=Wc.i1b)
 

Jc2 <- c(0,1,1)
Wc.e2 <- list(seq(1,Ne/3), seq(Ne/3+1, 2*Ne/3), seq(2*Ne/3+1,Ne))
Wc.i2 <- list(seq(Ne+1,Ne+Ni/2), seq(Ne+Ni/2+1,N), NULL)
pars2 <- list(Jc=Jc2, Wc.e=Wc.e2, Wc.i=Wc.i2)

Jc2b <- c(0,1,1)
Nee <- sample(1:Ne, Ne)
Wc.e2b <- list(Nee[1:(Ne/3)], Nee[(Ne/3+1): (2*Ne/3)], Nee[(2*Ne/3+1):Ne])
Wc.i2b <- list(seq(Ne+1,Ne+Ni/2), seq(Ne+Ni/2+1,N), NULL)
pars2b <- list(Jc=Jc2b, Wc.e=Wc.e2b, Wc.i=Wc.i2b)



# 1L - random initialization - connectivity (Wc.e) increasing - here the differences are caused by numerical errors
raster(X, dt=1)
newpars.1L <- init.hGLM(X, 1, pars, rand.seed=2, lin.soma=T, rescale.soma=T, nSD=1, sd.v=2, mean.v=10)
lines(newpars.1L$v, col=4, lty=1); abline(h=5)

epsp <- matrix(newpars.1L$v[1:(Ne*200)], ncol=Ne)
matplot(epsp[8:20,], t="l", col=rainbow(Ne), lt=1)

a.epsp <- apply(epsp, 2, max) - newpars.1L$v[1]
mini <- sim.hGLM(X, 1, newpars.1L, calc.minis=T, double=F)
plot(mini$a.e, a.epsp); abline(0, 1)
plot(exp(newpars.1L$Tau.e[[1]]), mini$a.e - a.epsp) 


# 1L - random initialization - connectivity (Wc.e) NOT increasing - still OK, just need to reorder
raster(X, dt=1)
newpars.1L <- init.hGLM(X, 1, pars1b, rand.seed=2, lin.soma=T, rescale.soma=T, nSD=1, sd.v=2, mean.v=10)
lines(newpars.1L$v, col=4, lty=1); abline(h=5)

epsp <- matrix(newpars.1L$v[1:(Ne*200)], ncol=Ne)
a.epsp <- apply(epsp, 2, max) - newpars.1L$v[1]
mini <- sim.hGLM(X, 1, newpars.1L, calc.minis=T, double=F)
plot(mini$a.e, a.epsp); abline(0, 1)
plot(mini$a.e, a.epsp[unlist(newpars.1L$Wc.e)]); abline(0, 1)
mini$a.e - a.epsp


# 2N - random initialization - connectivity (Wc.e) increasing
raster(X, dt=1)
newpars.2N <- init.hGLM(X, 1, pars2, rand.seed=2, lin.soma=F, rescale.soma=T, nSD=10)
lines(newpars.2N$v+5, col=4, lty=1); abline(h=5)

epsp <- matrix(newpars.2N$v[1:(Ne*200)], ncol=Ne)
a.epsp <- apply(epsp, 2, max) - newpars.2N$v[1]
mini <- sim.hGLM(X, 1, newpars.2N, calc.minis=T, double=F)
plot(mini$a.e, a.epsp)


# 2N - random initialization - connectivity (Wc.e) NOT increasing - still crap!
raster(X, dt=1)
newpars.2N <- init.hGLM(X, 1, pars2b, rand.seed=3, lin.soma=F, rescale.soma=T, nSD=10)
lines(newpars.2N$v+5, col=4, lty=1); abline(h=5)

epsp <- matrix(newpars.2N$v[1:(Ne*200)], ncol=Ne)
a.epsp <- apply(epsp, 2, max) - newpars.2N$v[1]
mini <- sim.hGLM(X, 1, newpars.2N, calc.minis=T, double=F)
plot(mini$a.e, a.epsp)
plot(mini$a.e, a.epsp[unlist(newpars.2N$Wc.e)], col=rep(1:3, each=5), pch=16); abline(0, 1)

