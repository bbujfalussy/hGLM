# Generate synaptic inputs with mixing different orientations
# we need 10 repetitions, each 48s long
# 16 orientations, each shown for 3s, following each other

# m.alpha is the probability of going into the up state depending on the stimulus orientation (and phase...)
#m.alphas <- c(12.40,  7.40,  3.76,  2.73,  4.47,  8.66, 13.53, 15.39, 12.40,  7.40,  3.76, 2.73,  4.47,  8.66, 13.53, 15.39)# - a sinusoid fit
# m.alphas <- c(12.64,  6.60,  2.65,  0.49, 12.64, 12.64, 15.87, 22.42, 22.42, 10.72,  1.00,  0.00,  2.08,  2.08,  4.51,  8.96) # - median of data

m.alphas <- c(10.54,  4.96,  1.07,  0.5,  1.82,  6.33, 11.84, 14.01, 10.54,  4.96,  1.07,  0.5,  1.82,  6.33, 11.84, 14.01)
beta <- 20
### generate data ...
### we generate using thinning - a method based on rejection sampling
### we first generate a Poisson train with rate_max = 10 Hz.
### the we remove all transitions with the probability of the ratio of the true rate and rate_max

source('./Funts_dendr_orisel.R', chdir = TRUE)
Tmax <- 3000 # set to 3000; 125
set.seed(317)
ori.dends <- floor(rnorm(16, 0, 1.5) %% 16)
print(ori.dends)

Ensyn <- c(48,58, 52,34,45,39,44, 68,50,62,30, 60,39)
Insyn <- c(11,11, 9,6,8,5,8, 12,11,13,6, 11,8); N.inh <- sum(Insyn)

Erate <- c(5, 20); Esd <- c(2.5, 10)
Irate <- c(20, 20)
N.soma <- 420

rep <- 1; ori <- 14; i.den <- 1

saveinput <- F
showinput <- T

for (rep in 1:2){

	# rate.inh <- rep(0, Tmax)	# inhibitory rate

	for (ori in 1:16){ # 16 stimulus orientation
		rate.inh <- rep(0, Tmax)	# inhibitory rate
		set.seed(100*rep + ori)
		print(ori)
		for (i.den in 1:13){ # preferred orientation of the 13 dendrites
			st <- gen.events.sin(Tmax, m.alphas[(ori-ori.dends[i.den]) %% 16 + 1], beta, maxL=150)
			# while(sum(st) < 1) st <- gen.events.sin(Tmax, m.alphas[(ori-ori.dends[i.den]) %% 16 + 1], beta, maxL=150)
			rate.inh <- rate.inh + Insyn[i.den] * (st * (Irate[2]-Irate[1]) + Irate[1])
			# plot(st, col=rainbow(16)[ori], t="l", lty=1, axes=F, xlab="", ylab="", ylim=c(-0.1,1.1))
			spt.Eden <- gen.spikes.states(Tmax=Tmax, N=Ensyn[i.den], mu=Erate, sd=Esd, 500, st+1)
			if (i.den>1) {
				spt.Eden[,1] <- spt.Eden[,1] + sum(Ensyn[1:(i.den-1)])
				spt.E <- rbind(spt.E, spt.Eden)
				ii <- sort(spt.E[,2], index.return=T)$ix
				spt.E <- spt.E[ii,]
			} else spt.E <- spt.Eden
			cat(i.den, " ")
		}
		rate.inh <- rate.inh / N.inh
		sp.i.d <- gen.Poisson.train.rescale(rate.inh, N.inh)
		sp.i.d[,1] <- sp.i.d[,1] + 1
		sp.i.soma <- gen.Poisson.spikes(rate.inh*N.soma)
		
		sp.i <- rbind(sp.i.d, sp.i.soma)
		ii <- sort(sp.i[,2], index.return=T)$ix
		spt.I <- sp.i[ii,]
		
		if (ori==1) {
			Espikes <- spt.E 
			Ispikes <- spt.I
		} else {
			spt.E[,2] <- spt.E[,2] + (ori-1) * Tmax
			spt.I[,2] <- spt.I[,2] + (ori-1) * Tmax
			Espikes <- rbind(Espikes, spt.E)
			Ispikes <- rbind(Ispikes, spt.I)
		}
		cat("orientation", ori, "finished. \n")
	}
				
 	if (showinput){ 	
		plot(Espikes[,2], Espikes[,1], pch=".", ylim=c(0, 750))
	 	points(Ispikes[,2], Ispikes[,1] + 629, pch=".", col=3)
		# plot(spt.I[,2], spt.I[,1], pch=".")
	}

	if (saveinput){ 	
	    fname <- paste("Espikes_d", Tmax*16, '_r', i.rate, '_rep', rep, '_Ne', sum(Ensyn), '_e', Erate[1], '_E', Erate[2], '.dat', sep="")
		write.table(Espikes, file=fname, row.names=FALSE, col.names=FALSE)
	
	    fname <- paste("Ispikes_d", Tmax*16, '_r', i.rate, '_rep', rep, '_Ni', sum(Insyn)+1, '_i', Irate[1], '_I', Irate[2], '.dat', sep="")
		write.table(Ispikes, file=fname, row.names=FALSE, col.names=FALSE)
	}
}