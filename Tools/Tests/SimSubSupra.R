##############################################
## A function to test the subthreshold simulation and initialisation of GLM models


source('../Sim_hGLM_sub.R', chdir=T)
source('../Sim_hGLM_sp.R', chdir = TRUE)
source('../../Utils/Spike_Basis.R', chdir = TRUE)

############################################################
# Testing this function
# init.hGLM(X, dt, pars, logpars=T, double=F, scale=2.8, lin.soma=F, rand.seed=12, nSD=4, X.ahp=NULL)

L <- 20000; dt <- 1
X <- matrix(rbinom(20*L, 1, 0.005), 20, L)

Jc3 <- c(0,1,1, 2,2, 3,3)
Wc.e3 <- list(NULL, 1,2, 3:4, 5:6, 7:8, 9:10)
Wc.i3 <- list(seq(11,15), 16, 17, 18, NULL, 19, 20)
pars3 <- list(Jc=Jc3, Wc.e=Wc.e3, Wc.i=Wc.i3)
linear.soma <- T

# 9. 3N - random initialization - try a couple of different
raster(X)
pars.sub <- init.hGLM(X, 1, pars3, rand.seed=5, lin.soma=linear.soma, nSD=1, double=T, rescale.soma=T)
lines(pars.sub$v+5, col=2, lty=1); abline(h=5)

##############################################################
## we need the following thing for simulating spikes + Vm after spikes
## 1. basis functions for Vm and Th
## 2. wights for the basis
## 3. simulate random spikes and Vm given the input
## 4. simulate random spikes given the Vm
## 5. simulate Vm given input + output spike times

##############################################################
## W.ahp and W.th must contain basis parameters as attributes - N, Tmax

basis.a <- gen.basis(Tmax=100, dt=1, n.basis=8, first=F, graphics=T)
Jw <- exp(pars.sub$Jw[1])
if (linear.soma) y.soma <- (pars.sub$v - pars.sub$v0)/ Jw[1] else y.soma <- isigm((pars.sub$v - pars.sub$v0)/Jw[1], c=pars.sub$Th[1])

pars.spike <- default.spikes(dt=dt, n.basis=8, Tmax=100, v.soma=pars.sub$v, pars.sub=pars.sub, graphics=T, rate=5)
pars.sub$W.ahp <- pars.spike$W.ahp
# dt <- 1; pars <- pars.sub; pars.spike <- pars.spike; rseed <- 317; ref.spikes <- NULL; logpars <- T; double <- T; scale <- 2.8; verbose=F
resp <- sim.hGLM.sp(X, 1, pars.sub, pars.spike, rseed=317, double=T, Nrep=10)


matplot(t(resp$v + resp$sp)[1:1000,], t="l", lty=1, col=1)
# plot(resp$v[1,] + resp$sp[1,], t="l")
lines(pars.sub$v, t="l", col=2)
# lines(y.soma, t="l", col=2)
lines(resp$	v[1,] + resp$sp[1,], t="l", col=3)
sum(resp$sp[1,]>0) / 20

plot(which(resp$sp[1,]>0), rep(1, sum(resp$sp[1,]>0),), pch="|", ylim=c(0, 11), xlim=c(0, L), axes=F, xlab="time (ms)", ylab="spikes"); axis(1); axis(2, las=2)
for (j in 2:10) points(which(resp$sp[j,]>0), rep(j, sum(resp$sp[j,]>0),), pch="|")
