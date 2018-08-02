##############################################
## A function to test the initialisation, simulation and training of hGLM models

source('./Tools/Train_hGLM.R', chdir=TRUE)

############################################################
# creating data

L <- 5000 # length of input data
dt <- 1 # time resolution
N <- 20 # number of pre neurons
X <- matrix(rbinom(N*L, 1, 0.005), N, L) # input spikes
X.test <- matrix(rbinom(N*L, 1, 0.005), N, L) # input spikes

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

linear.soma <- T

# plot the input spikes
raster(X, dt=1)
# initialising the model subthreshold part of the model
pars.sub <- init.hGLM(X, dt, pars3, rand.seed=5, lin.soma=linear.soma, nSD=1, double=T, rescale.soma=T, mean.v=-65, sd.v=2)
lines(pars.sub$v + 70, col=2, lty=1); abline(h=5)


## initialising the spikes
pars.spike <- default.spikes(dt=dt, n.basis=8, Tmax=200, v.soma=pars.sub$v, pars.sub=pars.sub, graphics=0, rate=5)
pars.sub$W.ahp <- pars.spike$W.ahp

# simulating the spiking response
resp <- sim.hGLM.sp(X, dt, pars.sub, pars.spike, rseed=1, double=T, Nrep=10)
lines(resp$v[1,] + resp$sp[1,]+70, t="l", lty=1, col=3)

## response + spikes - 10 repetitions
matplot(t(resp$v + resp$sp)+70, t="l", col=rainbow(10, end=0.75), lty=1, add=F)

## only the spike response
plot(which(resp$sp[1,]>0), rep(1, sum(resp$sp[1,]>0),), pch="|", ylim=c(0, 11), xlim=c(0, L), axes=F, xlab="time (ms)", ylab="spikes"); axis(1); axis(2, las=2)
for (j in 2:10) points(which(resp$sp[j,]>0), rep(j, sum(resp$sp[j,]>0),), pch="|")

##########################################
## without train and test arguments the train function can be used for validation - it generates training and test data for itself
opt.pars1 <- train.hGLM(pars=pars, graphics.on=2)
opt.pars2 <- train.hGLM(pars=pars, graphics.on=2, linear.soma=T)
opt.pars3 <- train.hGLM(pars=pars2, graphics.on=2, linear.soma=T)

# preinit often helps to start optimization if spikes are present
opt.pars4 <- train.hGLM(pars=pars3, graphics.on=2, linear.soma=T, sim.spikes=T, preinit=T)

##########################################
## with train and test arguments they are used as reference data
##########################################

## use a 3N model for training - subthreshold ...
pars.sub <- init.hGLM(X, dt, pars3, rand.seed=5, lin.soma=linear.soma, nSD=1, double=T, rescale.soma=T, mean.v=-65, sd.v=2)
v <- pars.sub$v
v.test <- sim.hGLM(X.test, dt=dt, pars=pars.sub, double=T)

train <- list(X=X, v=v)
test <- list(X=X.test, v=v.test)

## different hGLM architectures for explaining the data
opt.pars5 <- train.hGLM(pars=pars, train=train, test=test, dt=1, graphics.on=2, linear.soma=T, double=F)
opt.pars6 <- train.hGLM(pars=pars, train=train, test=test, dt=1, graphics.on=2, linear.soma=T, double=T)
opt.pars7 <- train.hGLM(pars=pars2, train=train, test=test, dt=1, graphics.on=2, linear.soma=T, double=T)
opt.pars8 <- train.hGLM(pars=pars3, train=train, test=test, dt=1, graphics.on=2, linear.soma=T, double=T, maxit=50)
opt.pars8b <- train.hGLM(pars=pars3, train=train, test=test, dt=1, graphics.on=2, linear.soma=T, double=T, maxit=500)

##########################################
## ... same with spikes
pars.spike <- default.spikes(dt=dt, n.basis=8, Tmax=200, v.soma=pars.sub$v, pars.sub=pars.sub, graphics=0, rate=5)
pars.sub$W.ahp <- pars.spike$W.ahp
resp <- sim.hGLM.sp(X, dt, pars.sub, pars.spike, rseed=1, double=T, Nrep=1)
v <- as.vector(resp$v + resp$sp)

resp.test <- sim.hGLM.sp(X.test, dt, pars.sub, pars.spike, rseed=1, double=T, Nrep=1)
v.test <- as.vector(resp.test$v + resp.test$sp)
matplot(cbind(v, v.test), t="l", ylim=c(-75, -40), lty=1)

train <- list(X=X, v=v)
test <- list(X=X.test, v=v.test)

## different hGLM architectures for explaining the data
opt.pars9 <- train.hGLM(pars=pars, train=train, test=test, dt=1, graphics.on=2, linear.soma=T, double=F, preinit=T)
opt.pars10 <- train.hGLM(pars=pars, train=train, test=test, dt=1, graphics.on=2, linear.soma=T, double=T, preinit=T)
opt.pars11 <- train.hGLM(pars=pars2, train=train, test=test, dt=1, graphics.on=2, linear.soma=T, double=T, preinit=T)
opt.pars12 <- train.hGLM(pars=pars3, train=train, test=test, dt=1, graphics.on=2, linear.soma=T, double=T, preinit=T)

## spike prediction accuracy can not be high if the reference cell is stochastic and we don't have multiple runs
## see opt.pars4$sp.rel 


