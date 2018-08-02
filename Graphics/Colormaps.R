

BlueRed <- function(n){
	if (n == 1) cols <- rgb(1,0,0)
	if (n == 2) cols <- c(rgb(1,0,0), rgb(0,0,1))
	if (n > 2) {
		n2 <- ceiling(n/2); n1 <- n-n2
		reds <- c(rep(1,n1-1), seq(1,0, length=n2+1))
		reds[1:(n1+1)] <- seq(1, reds[n1+1], length=n1+1)
		blues <- rev(reds)
		greens <- apply(rbind(reds, blues), 2, min)
		cols <- rgb(reds, greens, blues)		
	}		
	cols
}

BlueGreen <- function(n){
	if (n == 1) cols <- rgb(1,0,0)
	if (n == 2) cols <- c(rgb(1,0,0), rgb(0,0,1))
	if (n > 2) {
		n2 <- ceiling(n/2); n1 <- n-n2
		blues <- c(rep(1,n1-1), seq(1,0, length=n2+1))
		greens <- rev(blues)
		reds <- apply(rbind(greens, blues), 2, min)
		cols <- rgb(reds, greens, blues)		
	}		
	cols
}

# cols4 <- c(rgb(205/255, 221/255, 255/255), rgb(69/255, 131/255, 254/255), rgb(120/255, 218/255, 116/255), rgb(249/255, 78/255, 7/255))

# cols <- BlueRed(2)
# plot(seq(1,100), rnorm(100), pch=21, bg=cols)
# cols <- BlueRed(11)
# plot(seq(1,100), rnorm(100), pch=21, bg=cols)


# cols <- BlueRed(100)
# plot(seq(1,100), rnorm(100), pch=21, bg=cols)

# cols <- BlueRed(101)
# plot(seq(1,101), rnorm(101), pch=21, bg=cols)
# cols <- BlueGreen(101)
# plot(seq(1,101), rnorm(101), pch=21, bg=cols)

blueRed2 <- function(x){
	# x has to be between 0 and 1
	if (min(x) < 0){ 
		x <- x - min(x)
		x <- x / max(x)
	}
	if (max(x) > 1){ 
		x <- x - min(x)
		x <- x / max(x)
	}
	
	cols <- rep(NA, length(x))
	for (i in 1:length(x)){
		if (x[i]<0.5) cols[i] <- rgb(1, 2*x[i], 2*x[i]) else cols[i] <- rgb(2-2*x[i],2-2*x[i],1)
	}
	cols
}



BlackRed <- function(n){
	if (n == 1) cols <- rgb(1,0,0)
	if (n == 2) cols <- c(rgb(1,0,0), rgb(0,0,1))
	if (n > 2) {
		n2 <- ceiling(n/2); n1 <- n-n2
		reds <- c(seq(0,1, length=n2), rep(1,n1))
		blues <- c(rep(0,n2-1), seq(0, 0.8, length=n1+1))
		greens <- blues
		cols <- rgb(reds, greens, blues)		
	}		
	cols
}


BlackGreen <- function(n){
	if (n == 1) cols <- rgb(1,0,0)
	if (n == 2) cols <- c(rgb(1,0,0), rgb(0,0,1))
	if (n > 2) {
		n2 <- ceiling(n/2); n1 <- n-n2
		greens <- c(seq(0,1, length=n2), rep(1,n1))
		blues <- c(rep(0,n2-1), seq(0, 0.8, length=n1+1))
		reds <- blues
		cols <- rgb(reds, greens, blues)		
	}		
	cols
}

BlackBlue <- function(n){
	if (n == 1) cols <- rgb(1,0,0)
	if (n == 2) cols <- c(rgb(1,0,0), rgb(0,0,1))
	if (n > 2) {
		n2 <- ceiling(n/2); n1 <- n-n2
		blues <- c(seq(0,1, length=n2), rep(1,n1))
		reds <- c(rep(0,n2-1), seq(0, 0.8, length=n1+1))
		greens <- reds
		cols <- rgb(reds, greens, blues)		
	}		
	cols
}



