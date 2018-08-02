
extract.spiketimes <- function(voltage, sampling.freq, graphics=0, limite=0, dvdt.th=5){

   # Given a voltage trace 'voltage', detect the spiketimes with the upward
   # zero crossing cirterion and set spiketimes to the maximum value
   # following the zero-crossing in a 1 ms windows
   # limite: threshold for spike detection

   # Output: list with 
   # st: the index of the spike peaks
   # delay: the delay between threshold crossing and spike peak
   # dur: of the spikes from the threshold crossing to the first point below the threshold


    #################################################
    ## extract spikes maxima
    
    dt <- 1e3/sampling.freq;	# ms
    voltage.prime <- c(0,diff(voltage))  		# voltage derivative
    Dt <- round(2/dt)
    limit.t.last <- floor(2/dt)
    t.last <- limit.t.last + 1
    spiketimes.t <- 0
    
    for (i in 1:(length(voltage)-Dt)){         # loop over the voltage traces
        if (voltage[i] >= limite && voltage.prime[i] > 0 && t.last >= limit.t.last){
            t.last <- 0
            temp.ind <- which.max(voltage[i:(i+Dt)])
            spiketimes.t <- c(spiketimes.t, i+temp.ind)
			# plot(voltage[(i+temp.ind-5):(i+temp.ind+15)], t="l")
			# abline(v=5)
			# readline(i)
       } else t.last <- t.last + 1
    }
    
    spiketimes.t <- spiketimes.t-1
    st <- spiketimes.t[-1]

    #################################################
    ## extract parameters of the spikes
	if (length(st) > 0){
	    nbr.spikes <- length(st)
		size.spike.shape <- 30/dt;          # size of the spike shape (30 ms) - the peak is at 5 + dt!
		
		temp.spikes <- st - 5/dt             # set the spiketimes 5 ms before the maximum of the AP
		tt.spikes <- c(temp.spikes, length(voltage))
		spike.shape <- matrix(NA, size.spike.shape,nbr.spikes)     # i.e. just used for plotting
		k<-1
		for (i in 1:nbr.spikes){
			# only spikes not followed by another spike within the test period ...
		    if ((tt.spikes[i] + size.spike.shape <= tt.spikes[i+1]) & (tt.spikes[i]>0)){ 
			    	spike.shape[,k] <- voltage[tt.spikes[i]:(tt.spikes[i]+size.spike.shape-1)]
		        k <- k+1
		    }
		}
			
		m.spike.shape <- rowMeans(spike.shape, na.rm=T)             # Compute the mean spike shape
		i.threshold <- min(which((diff(m.spike.shape)/dt) > dvdt.th)) # dV/dt > 5 mV/ms - start of the spike
		t.threshold <- i.threshold*dt # dV/dt > 5 mV/ms - start of the spike
		# end.spike <- (max(which((diff(m.spike.shape)/dt) < -5))+1)*dt # dV/dt > -5 mV/ms - end of the spike
		V.th <- m.spike.shape[t.threshold/dt]
			
		i.end.spike <- min(which(m.spike.shape[i.threshold:size.spike.shape] < V.th)) + i.threshold - 1
		end.spike <- i.end.spike * dt # voltage is below threshold
		
		delay.spike <- 5 + dt - t.threshold # ms - this is the threshold relative to the peak (positive)!
		dur.spike <- end.spike - t.threshold # ms - this is the first normal point relative to the threshold
		n.shapes <- sum(!is.na(spike.shape[1,]))
		
		if (graphics>1)	{
			t.spikes <- seq(dt, length=size.spike.shape, by=dt)
			matplot(t.spikes, spike.shape, t="p", lty=1, pch=23, cex=0.3, axes=F, col=rainbow(n.shapes)); axis(1); axis(2, las=2)
			lines(t.spikes, m.spike.shape, col=1, lwd=3)
			lines(t.spikes[(t.threshold/dt): (end.spike/dt)], m.spike.shape[(t.threshold/dt): (end.spike/dt)], col=2, lwd=3)
			abline(v=c(t.threshold, end.spike), col=c(1,2))
			text(t.threshold, max(m.spike.shape), "threshold", pos=1, col=1)
			text(end.spike, min(V.th), "end of spike", pos=1, col=2)
		}
	
	    out  <- list(st=st, delay=delay.spike, dur=dur.spike)
	} else out  <- list(st=NULL, delay=NA, dur=NA)
	out
}


remove.spikes <- function(X, spt, dt){
	i.ths <- spt$st - spt$delay/dt
	n.cut <- spt$dur/dt
	## i.th: indices of the first points to cut
	## n.cut: number of points to cut
	
	if (n.cut<1) stop("n.cut must be greater than 0")
	if (is.vector(X)) X <- matrix(X, 1, length(X))
	nx <- ncol(X)
	for (i.th in i.ths){
		# spt=10; delay=-2; dur=5; [10-2]-[10-2+5]; cut 9-10-11-12
		X[,(max(i.th,1)):(min(i.th + n.cut -1, nx))] <- NA
	}
	X
}



cut.spikes <- function(X, st, n.cut){
	## st: indices of the first points to cut
	## n.cut: number of points to cut
	if (n.cut<1) stop("n.cut must be greater than 0")
	if (is.vector(X)) X <- matrix(X, 1, length(X))
	nx <- ncol(X)
	for (spt in st){
		# spt=10; delay=-2; dur=5; [10-2]-[10-2+5]; cut 9-10-11-12
		X[,(max(spt,1)):(min((spt + n.cut[1] -1), nx))] <- NA
	}
	X <- X[,!is.na(X[1,])]
	X
}

# X <- c(seq(1,10), 100, seq(12, 20))
# st <- 10
# cut.spikes(X, st, 0)
# cut.spikes(X, st, 1)
# cut.spikes(X, st, 2)

