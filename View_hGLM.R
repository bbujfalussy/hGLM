## A script to plot the response of various hLN models

setwd('~/Programs/Code/hGLM/')
source('./Tools/Train_hGLM.R', chdir=TRUE)
source('~/Programs/hGLMs/Rlib/Graphics/scalebar.R')

graphics <- F
refit <- F

group <- 'allDend'
stages <- c('Passive', 'Adend', 'Active')
stage <-stages[1]
alpha.Ww <- 10 # 50 
alpha.Tau <- 20 # 160

Tmax <- 48000
Tinit <- 0	

dt <- 1 # ms
Lmax <- (Tmax-Tinit) / dt
time <- seq(dt, (Tmax-Tinit), by=dt)

rep <- 7
source('./Data_IVL/read_data.R')
plot(resp.test, t="l")
cat(length(extract.spiketimes(resp.test, 1000, T, limite=-20)$st))

########################################
## initial parameters 
########################################

outdir <- paste('Data_IVL/Fits', stage, '/', sep="")

errfile <- paste(outdir, 'error_alphaW', alpha.Ww, 'T', alpha.Tau, '_n13_r2_r', rep, '.RData', sep='')
load(errfile)

## pooling spikes for the different models - 1 synaps / subunit
xefile.train <- paste('Data_IVL/xe_n13_r2_r', rep, '.RData', sep='')
xefile.test <- paste('Data_IVL/xe_n13_r2_r', 11-rep, '.RData', sep='')
source('./Data_IVL/modelStructures.R') # this can be somewhat slow...

########################################
## checking previously learned parameters 
########################################
## first check a couple of relevant models
## coupled models
if (stage != 'Active'){
	model <- errs$pars.1LDE
	err <- sim.hGLM(test.1N$X, dt=dt, pars=model$pars.sub, vv=test.1N$v, double=T)
	cat(1 - sqrt(err / length(test.1N$v) / var(test.1N$v)), model$s.exp[1,2])
	
	model <- errs$pars.2NDE
	err <- sim.hGLM(test.2N$X, dt=dt, pars=model$pars.sub, vv=test.2N$v, double=T)
	cat(1 - sqrt(err / length(test.2N$v) / var(test.2N$v)), model$s.exp[1,2])
	
	model <- errs$pars.3NDE
	err <- sim.hGLM(test.3N$X, dt=dt, pars=model$pars.sub, vv=test.3N$v, double=T)
	cat(1 - sqrt(err / length(test.3N$v) / var(test.3N$v)), model$s.exp[1,2])
	
	
	## ultimately, we are interested in the following models: 
	## uncoupled models
	
	# model <- errs$pars.1Lsyn3DE
	# err <- sim.hGLM(test.syn$X, dt=dt, pars=model$pars.sub, vv=test.syn$v, double=T)
	# cat(1 - sqrt(err / length(test.syn$v) / var(test.syn$v)), model$s.exp[1,2])
	
	model <- errs$pars.3Nsyn3DE
	# err <- sim.hGLM(test.syn$X, dt=dt, pars=model$pars.sub, vv=test.syn$v, double=T)
	# cat(1 - sqrt(err / length(test.syn$v) / var(test.syn$v)), model$s.exp[1,2])	
	resp <- sim.hGLM(test.syn$X, dt=dt, pars=model$pars.sub, double=T)
	
	t1 <- 12000; tdur <- 10000; t2 <- t1 + tdur


	outfile <- paste('predresp_', stage, '_rep', rep, '_t', t1, '.pdf', sep='')
	pdf(file=outfile, w=6, h=2, useDingbats=F)
	par(mar=c(1,1,1,1))
	matplot(cbind(test.syn$v[t1:t2], resp[t1:t2]), t="l", lty=1, axes=F, xlab="", ylab="")
	scalebar2(500, 5, '0.5 s', '5 mV', pos='topleft')
	dev.off()


}

#######################################################
## spiking models
if (stage == 'Active'){
	model <- errs$pars.1LDE
	resp <- sim.hGLM.sp(test.1N$X, dt=dt, pars=model$pars.sub, pars.spike=model$pars.spike, double=T, Nrep=10, verbose=2)
	matplot(t(resp$v + resp$sp), t="l", col=grey(0.5), lty=1)	
	lines(test.1N$v, col=2)
	
	model <- errs$pars.2NDE
	resp <- sim.hGLM.sp(test.2N$X, dt=dt, pars=model$pars.sub, pars.spike=model$pars.spike, double=T, Nrep=10, verbose=1)
	matplot(t(resp$v + resp$sp), t="l", col=grey(0.5), lty=1)	
	lines(test.1N$v, col=2)
	
	model <- errs$pars.3NDE
	resp <- sim.hGLM.sp(test.3N$X, dt=dt, pars=model$pars.sub, pars.spike=model$pars.spike, double=T, Nrep=10, verbose=2)
	matplot(t(resp$v + resp$sp), t="l", col=grey(0.5), lty=1)	
	lines(test.1N$v, col=2)
	
	## ultimately, we are interested in the following models: 
	## uncoupled models
	
	# model <- errs$pars.1Lsyn3DE
	# resp <- sim.hGLM.sp(test.syn$X, dt=dt, pars=model$pars.sub, pars.spike=model$pars.spike, double=T, Nrep=10, verbose=2)
	# matplot(t(resp$v + resp$sp), t="l", col=grey(0.5), lty=1)	
	# lines(test.1N$v, col=2)
	
	model <- errs$pars.3Nsyn3DE
	resp <- sim.hGLM.sp(test.syn$X, dt=dt, pars=model$pars.sub, pars.spike=model$pars.spike, double=T, Nrep=10, verbose=2)


	t1 <- 12000; tdur <- 10000; t2 <- t1 + tdur
	matplot(t(resp$v[,t1:t2] + resp$sp[,t1:t2]), t="l", col=rainbow(10), lty=1)	
	lines(test.1N$v[t1:t2], col=1)
	
	outfile <- paste('predresp_', stage, '_rep', rep, '_t', t1, '.pdf', sep='')
	pdf(file=outfile, w=6, h=2, useDingbats=F)
	par(mar=c(1,1,1,1))
	k <- 4
	v <- resp$v[k,] + resp$sp[k,]
	matplot(cbind(test.syn$v[t1:t2], v[t1:t2]), t="l", lty=1, axes=F, xlab="", ylab="", lwd=c(2, 1))
	scalebar2(500, 20, '0.5 s', '20 mV', pos='topleft')	
	dev.off()

}