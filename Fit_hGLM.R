## A script to train hLN model on compartmental data
## PAR1 <- 4; PAR2 <- 1


setwd('/data/bbu/hGLM/')
source('./Tools/Train_hGLM.R', chdir=TRUE)


graphics <- 0
refit <- F

group <- 'allDend'
stages <- c('Passive', 'Adend', 'Active', 'AdendRand')
stage <-stages[PAR1]
alpha.Ww <- 10 # 50 
alpha.Tau <- 20 # 160

preinit <- F
graphics <- 0

Tmax <- 48000
Tinit <- 0 #40000	

dt <- 1 # ms
Lmax <- (Tmax-Tinit) / dt
time <- seq(dt, (Tmax-Tinit), by=dt)

rep <- PAR2
source('./Data_IVL/read_data.R')


########################################
## initial parameters 
########################################

outdir <- paste('Data_IVL/Fits', stage, '/', sep="")
initdir <- paste('Data_IVL/FitsAdend/', sep="")

outfile <- paste(outdir, 'output_', stage, '_alphaW', alpha.Ww, 'T', alpha.Tau, '_n13_r2_r', rep, '.txt', sep='')
write('initiation done, evaluating the results of previous optimization ... ', file=outfile)

errfile <- paste(outdir, 'error_alphaW', alpha.Ww, 'T', alpha.Tau, '_n13_r2_r', rep, '.RData', sep='')
initfile <- paste(initdir, 'error_alphaW', alpha.Ww, 'T', alpha.Tau, '_n13_r2_r', rep, '.RData', sep='')

if (file.exists(errfile)){
	cat('errorfile loaded ... \n')
	load(errfile)
	names(errs)
	for (i in 1:length(errs)) {
		if ('s.exp' %in% names(errs[[i]]))	write(paste(names(errs)[i], '\t', round(errs[[i]]$s.exp[2,1], 2), '\t', round(errs[[i]]$s.exp[2,2], 4)), file=outfile, append=T)
	}
} else {
	cat('file loaded for initialization... \n')
	load(initfile)
	names(errs)
	for (i in 1:length(errs)) {
		if ('s.exp' %in% names(errs[[i]]))	write(paste(names(errs)[i], '\t', round(errs[[i]]$s.exp[2,1], 2), '\t', round(errs[[i]]$s.exp[2,2], 4)), file=outfile, append=T)
	}
	if (stage == 'Active')  preinit <- T
}

write('starting optimization ... \n ', file=outfile, append=T)
print(errfile)

## pooling spikes for the different models - 1 synaps / subunit
if (stage == 'AdendRand') {
	xefile.train <- paste('Data_IVL/xeRand_n13_r2_r', rep, '.RData', sep='')
	xefile.test <- paste('Data_IVL/xeRand_n13_r2_r', 11-rep, '.RData', sep='')
} else {
	xefile.train <- paste('Data_IVL/xe_n13_r2_r', rep, '.RData', sep='')
	xefile.test <- paste('Data_IVL/xe_n13_r2_r', 11-rep, '.RData', sep='')
}
source('./Data_IVL/modelStructures.R') # this can be somewhat slow...

########################################
## training using parameters learned for the stage Adend
########################################
## model architectures - for 1LDE - 3NDE models

# initpars1 <- list(Jc=Jc.1N, Wc.e=Wc.e.1N, Wc.i=Wc.i.1N)
# initpars2 <- list(Jc=Jc.2N, Wc.e=Wc.e.2N, Wc.i=Wc.i.2N)
# initpars3 <- list(Jc=Jc.3N, Wc.e=Wc.e.3N, Wc.i=Wc.i.3N)

# ipv1 <- init.hGLM(train.1N$X, dt=dt, pars=initpars1, double=T, mean.v=median(train.1N$v), sd.v=sd(train.1N$v))
# ipv2 <- init.hGLM(train.2N$X, dt=dt, pars=initpars2, double=T, mean.v=median(train.1N$v), sd.v=sd(train.1N$v))
# ipv3 <- init.hGLM(train.3N$X, dt=dt, pars=initpars3, double=T, mean.v=median(train.1N$v), sd.v=sd(train.1N$v))

### $opt should be changed to $pars.sub everywhere!
# ppars <- reshapeParList(ipv1, errs$pars.1LDE$opt)
ppars <- errs$pars.1LDE$pars.sub
errs$pars.1LDE <- train.hGLM(ppars, train=train.1N, test=test.1N, double=T, linear.soma=T, init=F, preinit=preinit, graphics.on=graphics)
write(paste("training signal explained: ", errs$pars.1LDE$s.exp[2,2], "test signal explained: ", errs$pars.1LDE$s.exp[1,2], "\n"), file=outfile, append=T)

## 1NDE
# ppars <- reshapeParList(ipv1, errs$pars.1NDE$opt)
ppars <- errs$pars.1NDE$pars.sub
errs$pars.1NDE <- train.hGLM(ppars, train=train.1N, test=test.1N, double=T, linear.soma=F, init=F, preinit=preinit, graphics.on=graphics)
write(paste("training signal explained: ", errs$pars.1NDE$s.exp[2,2], "test signal explained: ", errs$pars.1NDE$s.exp[1,2], "\n"), file=outfile, append=T)

## 2NDE
# ppars <- reshapeParList(ipv2, errs$pars.2NDE$opt)
# ppars <- check.nulls(ppars)
ppars <- errs$pars.2NDE$pars.sub
errs$pars.2NDE <- train.hGLM(ppars, train=train.2N, test=test.2N, double=T, linear.soma=F, init=F, preinit=preinit, graphics.on=graphics)
save(errs, file=errfile)
write(paste("training signal explained: ", errs$pars.2NDE$s.exp[2,2], "test signal explained: ", errs$pars.2NDE$s.exp[1,2], "\n"), file=outfile, append=T)

## 3NDE
# ppars <- reshapeParList(ipv3, errs$pars.3NDE$opt)
# ppars <- check.nulls(ppars)
ppars <- errs$pars.3NDE$pars.sub
errs$pars.3NDE <- train.hGLM(ppars, train=train.3N, test=test.3N, double=T, linear.soma=F, init=F, preinit=preinit, graphics.on=graphics)
save(errs, file=errfile)
write(paste("training signal explained: ", errs$pars.3NDE$s.exp[2,2], "test signal explained: ", errs$pars.3NDE$s.exp[1,2], "\n"), file=outfile, append=T)

#################################################################
## uncoupled models
#################################################################

init <- 'coupled'

#################################################################
## 1Lsyn3DE
## initialising from the 1LDE model:
if (init == "coupled"){
	# sim.hGLM(train.1N$X, dt, errs$pars.1LDE$pars.sub, train.1N$v, double=T)
	initpar <- errs$pars.1LDE$pars.sub
	minis <- sim.hGLM(X=train.1N$X[,1:2000], dt=dt, pars=initpar, calc.minis=T, double=T)
	regpars <- regpars.from.coupled(Wc.e.1Nsyn, Wc.i.1Nsyn, minis, alpha.Ww, alpha.Tau)
	initpars.1L <- decouple.pars(Jc=Jc.1N, Wc.e=Wc.e.1Nsyn, Wc.i=Wc.i.1Nsyn, initpar) 
	write("1Lsyn3DE model initialised from coupled", file=outfile, append=T)
	# sim.hGLM(train.syn$X, dt, initpars.1L, train.syn$v, double=T)
}

## initialising from the Adend simulations
if (init == 'Adend'){
	if ('regpars' %in% names(errs$pars.1Lsyn3DE)) {
		regpars <- errs$pars.1Lsyn3DE$regpars 
	} else {
		initpar <- errs$pars.1LDE$pars.sub
		minis <- sim.hGLM(X=train.1N$X[,1:2000], dt=dt, pars=initpar, calc.minis=T, double=T)
		regpars <- regpars.from.coupled(Wc.e.1Nsyn, Wc.i.1Nsyn, minis, alpha.Ww, alpha.Tau)
	}
	initpars.1L <- errs$pars.1Lsyn3DE$pars.sub
}

## optimization
runfit <- T
k <- 1

while(runfit){
	pars.1Lsyn3DE <- train.hGLM(initpars.1L, train=train.syn, test=test.syn, dt=dt, graphics.on=F,  linear.soma=T, init=F, preinit=preinit, regpars=regpars, maxit=10)
	pars.1Lsyn3DE$regpars <- regpars
	
	errs$pars.1Lsyn3DE <- pars.1Lsyn3DE
	save(errs, file=errfile)

	initpars.1L <- errs$pars.1Lsyn3DE$pars.sub
	if (pars.1Lsyn3DE$convergence==0) runfit <- F
	diff.sig <- diff(pars.1Lsyn3DE$s.exp[,1]) / pars.1Lsyn3DE$s.exp[2,1]
	if (diff.sig < 1/10000) runfit <- F
	
	write(paste("k=", k, " , diff=", round(diff.sig, 5), " , training sig. exp.=", round(pars.1Lsyn3DE$s.exp[2,1],5), sep=""), file=outfile, append=T)
	write(paste("current signal explained: ", pars.1Lsyn3DE$s.exp[2,2], "\n"), file=outfile, append=T)
	k <- k + 1
}
write("1Lsyn3DE model finished \n", file=outfile, append=T)

#################################################################
## 1Nsyn3DE
# regpar from the 1NDE model:
if (init == "coupled"){
	initpar <- errs$pars.1NDE$pars.sub
	minis <- sim.hGLM(X=train.1N$X[,1:2000], dt=dt, pars=initpar, calc.minis=T, double=T)
	regpars <- regpars.from.coupled(Wc.e.1Nsyn, Wc.i.1Nsyn, minis, alpha.Ww, alpha.Tau)
	initpars.1N <- decouple.pars(Jc=Jc.1N, Wc.e=Wc.e.1Nsyn, Wc.i=Wc.i.1Nsyn, initpar) 
	write("1Nsyn3DE model initialised from coupled", file=outfile, append=T)
}

## initialising from the Adend simulations
if (init == 'Adend'){
	if ('regpars' %in% names(errs$pars.1Nsyn3DE)) {
		regpars <- errs$pars.1Nsyn3DE$regpars 
	} else {
		initpar <- errs$pars.1NDE$pars.sub
		minis <- sim.hGLM(X=train.1N$X[,1:2000], dt=dt, pars=initpar, calc.minis=T, double=T)
		regpars <- regpars.from.coupled(Wc.e.1Nsyn, Wc.i.1Nsyn, minis, alpha.Ww, alpha.Tau)
	}
	initpars.1N <- errs$pars.1Nsyn3DE$pars.sub
}

## optimization
runfit <- T
k <- 1

while(runfit){
	pars.1Nsyn3DE <- train.hGLM(initpars.1N, train=train.syn, test=test.syn, dt=dt, graphics.on=F,  linear.soma=F, init=F, preinit=preinit, regpars=regpars, maxit=10)
	pars.1Nsyn3DE$regpars <- regpars
	
	errs$pars.1Nsyn3DE <- pars.1Nsyn3DE
	save(errs, file=errfile)

	initpars.1N <- errs$pars.1Nsyn3DE$pars.sub
	if (pars.1Nsyn3DE$convergence==0) runfit <- F
	diff.sig <- diff(pars.1Nsyn3DE$s.exp[,1]) / pars.1Nsyn3DE$s.exp[2,1]
	if (diff.sig < 1/10000) runfit <- F
	
	write(paste("k=", k, " , diff=", round(diff.sig, 5), " , training sig. exp.=", round(pars.1Nsyn3DE$s.exp[2,1],5), sep=""), file=outfile, append=T)
	write(paste("current signal explained: ", pars.1Nsyn3DE$s.exp[2,2], "\n"), file=outfile, append=T)
	k <- k + 1
}
write("1Nsyn3DE model finished \n", file=outfile, append=T)


#################################################################
## 2Nsyn3DE
## initialization: 
## 1. based on Adend - use its own regularisation
## 2. based on 1Nsyn - use the same regularisation
## 3. based on 2Lsyn - use the same regularisation
## 4. based on 2NDE - regularise to the 2NDE
## Here we implement only 1 and 4. 

# regpar from the 2NDE model:
if (init == "coupled"){
	initpar <- errs$pars.2NDE$pars.sub
	minis <- sim.hGLM(X=train.2N$X[,1:2000], dt=dt, pars=initpar, calc.minis=T, double=T)
	regpars <- regpars.from.coupled(Wc.e.2Nsyn, Wc.i.2Nsyn, minis, alpha.Ww, alpha.Tau)
	initpars.2N <- decouple.pars(Jc=Jc.2N, Wc.e=Wc.e.2Nsyn, Wc.i=Wc.i.2Nsyn, initpar) 
	write("2Nsyn3DE model initialised from coupled", file=outfile, append=T)
}

## initialising from the Adend simulations
if (init == 'Adend'){
	if ('regpars' %in% names(errs$pars.2Nsyn3DE)) {
		regpars <- errs$pars.2Nsyn3DE$regpars 
	} else {
		initpar <- errs$pars.2NDE$pars.sub
		minis <- sim.hGLM(X=train.2N$X[,1:2000], dt=dt, pars=initpar, calc.minis=T, double=T)
		regpars <- regpars.from.coupled(Wc.e.2Nsyn, Wc.i.2Nsyn, minis, alpha.Ww, alpha.Tau)
	}
	initpars.2N <- errs$pars.2Nsyn3DE$pars.sub
}

## optimization
runfit <- T
k <- 1

while(runfit){
	pars.2Nsyn3DE <- train.hGLM(initpars.2N, train=train.syn, test=test.syn, dt=dt, graphics.on=F,  linear.soma=F, init=F, preinit=preinit, regpars=regpars, maxit=10)
	pars.2Nsyn3DE$regpars <- regpars
	
	errs$pars.2Nsyn3DE <- pars.2Nsyn3DE
	save(errs, file=errfile)

	initpars.2N <- errs$pars.2Nsyn3DE$pars.sub
	if (pars.2Nsyn3DE$convergence==0) runfit <- F
	diff.sig <- diff(pars.2Nsyn3DE$s.exp[,1]) / pars.2Nsyn3DE$s.exp[2,1]
	if (diff.sig < 1/10000) runfit <- F
	
	write(paste("k=", k, " , diff=", round(diff.sig, 5), " , training sig. exp.=", round(pars.2Nsyn3DE$s.exp[2,1],5), sep=""), file=outfile, append=T)
	write(paste("current signal explained: ", pars.2Nsyn3DE$s.exp[2,2], "\n"), file=outfile, append=T)
	k <- k + 1
}
write("2Nsyn3DE model finished \n", file=outfile, append=T)




################################################################
# 3Nsyn3DE
# initialization: 
# regpar from the 3NDE model:
if (init == "coupled"){
	sim.hGLM(train.3N$X, dt, errs$pars.3NDE$pars.sub, train.3N$v, double=T)
	initpar <- errs$pars.3NDE$pars.sub
	minis <- sim.hGLM(X=train.3N$X[,1:2000], dt=dt, pars=initpar, calc.minis=T, double=T)
	regpars <- regpars.from.coupled(Wc.e.3Nsyn, Wc.i.3Nsyn, minis, alpha.Ww, alpha.Tau)

	initpars.3N <- decouple.pars(Jc=Jc.3N, Wc.e=Wc.e.3Nsyn, Wc.i=Wc.i.3Nsyn, initpar) 
	write("3Nsyn3DE model initialised from coupled", file=outfile, append=T)
	sim.hGLM(train.syn$X, dt, initpars.3N, train.syn$v, double=T)
}

# initialising from the Adend simulations
if (init == 'Adend'){
	if ('regpars' %in% names(errs$pars.3Nsyn3DE)) {
		regpars <- errs$pars.3Nsyn3DE$regpars 
	} else {
		initpar <- errs$pars.3NDE$pars.sub
		minis <- sim.hGLM(X=train.3N$X[,1:2000], dt=dt, pars=initpar, calc.minis=T, double=T)
		regpars <- regpars.from.coupled(Wc.e.3Nsyn, Wc.i.3Nsyn, minis, alpha.Ww, alpha.Tau)
	}
	initpars.3N <- errs$pars.3Nsyn3DE$pars.sub
}

# optimization
runfit <- T
k <- 1

train.syn$X <- train.syn$X[,1:5000]
train.syn$v <- train.syn$v[1:5000]

while(runfit){
	pars.3Nsyn3DE <- train.hGLM(initpars.3N, train=train.syn, test=test.syn, dt=dt, graphics.on=2,  linear.soma=F, init=F, preinit=preinit, regpars=regpars, maxit=10)
	pars.3Nsyn3DE$regpars <- regpars
	
	errs$pars.3Nsyn3DE <- pars.3Nsyn3DE
	save(errs, file=errfile)

	initpars.3N <- errs$pars.3Nsyn3DE$pars.sub
	if (pars.3Nsyn3DE$convergence==0) runfit <- F
	diff.sig <- diff(pars.3Nsyn3DE$s.exp[,1]) / pars.3Nsyn3DE$s.exp[2,1]
	if (diff.sig < 1/10000) runfit <- F
	
	write(paste("k=", k, " , diff=", round(diff.sig, 5), " , training sig. exp.=", round(pars.3Nsyn3DE$s.exp[2,1],5), sep=""), file=outfile, append=T)
	write(paste("current signal explained: ", pars.3Nsyn3DE$s.exp[2,2], "\n"), file=outfile, append=T)
	k <- k + 1
}
write("3Nsyn3DE model finished \n", file=outfile, append=T)

