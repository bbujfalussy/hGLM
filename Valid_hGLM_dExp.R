##############################################
## A function to test the initialisation, simulation and training of hGLM models

source('./Tools/Train_hGLM.R', chdir=TRUE)

############################################################
# creating data

L <- 300 # length of input data
dt <- 1 # time resolution
N <- 1 # number of pre neurons
X <- matrix(0, N, L) # input spikes
X[1, 50] <- 1

Jc <- c(0)
Wc.e <- list(1)
Wc.i <- list(NULL)
pars <- list(Jc=Jc, Wc.e=Wc.e, Wc.i=Wc.i)
 linear.soma <- T

# plot the input spikes
raster(X, dt=1)
# initialising the model subthreshold part of the model
pars.sub <- init.hGLM(X, dt, pars, rand.seed=6, lin.soma=linear.soma, nSD=1, double=F, rescale.soma=T, mean.v=-65, sd.v=0.2)
lines(pars.sub$v + 65.5, col=2, lty=1); abline(h=5)

train <- list(X=X, v=pars.sub$v)

## different hGLM architectures for explaining the data
opt.pars <- train.hGLM(pars=pars, train=train, dt=1, graphics.on=2, linear.soma=T, double=F)
opt.pars2 <- train.hGLM(pars=pars, train=train, dt=1, graphics.on=2, linear.soma=T, double=F, alphasyn=F)
