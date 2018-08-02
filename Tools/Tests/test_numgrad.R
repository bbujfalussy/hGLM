
## collect arguments into a list, and parameters into a vector
args.parvec <- make.args.parvec(pars=newpars, v=newpars$v, X, 1, X.ahp=X.ahp, optimvec=optimvec, regpars=regpars)
parvec.new <- args.parvec$parvec; args.newpars.newv <- args.parvec$args
args.newpars.refv <- args.newpars.newv; args.newpars.refv$v <- refpars$v

args.parvec <- make.args.parvec(pars=refpars, v=refpars$v, X, 1, X.ahp=X.ahp, optimvec=optimvec, regpars=regpars)
parvec.ref <- args.parvec$parvec; args.refpars.refv <- args.parvec$args



v.new <- sim.hGLM(X, 1, newpars, X.ahp=X.ahp, double=double, alphasyn=alphasyn)
vv.new <- err.hGLM(parvec.new, args.newpars.newv, ret.resp=T)
pp <- err.hGLM(parvec.new, args.newpars.newv, ret.parlist=T)
vvv.new <- sim.hGLM(X, 1, pp, X.ahp=X.ahp, double=double, alphasyn=alphasyn)

vv.ref <- err.hGLM(parvec.ref, args.refpars.refv, ret.resp=T)

## 1. plot: the responses
if (graphics){
	par(mfcol=c(2,2))
	# par(mfcol=c(1,1))
	raster(X)
	title("responses")
	lines(refpars$v+5, col=1, lwd=3); abline(h=5)
	lines(vv.ref+5, col=6, lty=1)
	lines(v.new+5, col=2, lwd=3)
	lines(vvv.new+5, col=5, lty=1)
	lines(vv.new+5, col=3, lty=1)
	legend('topright', legend=c('ref', 'ref - err', 'new - sim', 'new - err'), col=c(1,6,2, 3), lty=c(1,1,1,1), lwd=c(3,1,3,1), bg=grey(1))
}

cat('testing error ... \n')

if (max((v.new- vv.new)^2) > 1e-10) warning("sim and err is not doing the same!")
if (max((v.new- vvv.new)^2) > 1e-10) warning("sim and sim(err) is not doing the same!")
if (max((refpars$v - vv.ref)^2) > 1e-10) warning("reference voltage has changed")

err.nn <- sim.hGLM(X, 1, newpars, vv=newpars$v, double=double, X.ahp=X.ahp, regpars=regpars, alphasyn=alphasyn)
err.nr <- sim.hGLM(X, 1, newpars, vv=refpars$v, double=double, X.ahp=X.ahp, regpars=regpars, alphasyn=alphasyn)

errerr.nn <- err.hGLM(parvec.new, args.newpars.newv)
errerr.nr <- err.hGLM(parvec.new, args.newpars.refv)

if (abs(err.nn-errerr.nn) > 1e-10) warning("sim and err is not doing the same - ref!")
if (abs(err.nr-errerr.nr) > 1e-10) warning("sim and err is not doing the same - new!")


cat('testing the gradient ... \n')

gv.newpars.newv <- grad.err.hGLM(parvec.new, args.newpars.newv)
gv.newpars.refv <- grad.err.hGLM(parvec.new, args.newpars.refv)
gv.refpars.refv <- grad.err.hGLM(parvec.ref, args.refpars.refv)
# plot(gv.newpars.refv, gv.newpars.newv, xlim=range(gv.newpars.refv), ylim=range(gv.newpars.refv)); abline(h=0, col=2)
if (max(abs(gv.refpars.refv)) > 1e-5) warning("gradient is too large in the reference")
if (is.null(regpars)) {
	if (max(abs(gv.newpars.newv)) > 1e-5) warning("gradient is too large in the optimum - without regularisation")
}


gv.newpars.newv.num <- grad(err.hGLM, parvec.new, method="simple", args=args.newpars.newv, method.args=list(eps=1e-8))
gv.newpars.refv.num <- grad(err.hGLM, parvec.new, method="simple", args=args.newpars.refv, method.args=list(eps=1e-8))
gv.refpars.refv.num <- grad(err.hGLM, parvec.ref, method="simple", args=args.refpars.refv, method.args=list(eps=1e-8))

## look at the plots, and judge whether it works...
diff.grads.newpars.newv <- abs(gv.newpars.newv - gv.newpars.newv.num) / (abs(gv.newpars.newv)+1)
if(max(diff.grads.newpars.newv) > 1/10) {
	warning('gradient error - numerical and analytic gradients differs, when the match between target and prediction is perfect \n large gradients can be caused by the regularisation term in this case \n differences between the analytic and numeric are expected since the analytic gradients are calculated by \n assuming linear subunits (first order approximations) \n gradients caused by regularisation are therefore especially sensitive to the subunit nonlinearity parameters (Th and Jw) \n  gradient estimation must be much better when the target and prediction do not match \n')
}

diff.grads.newpars.refv <- abs(gv.newpars.refv - gv.newpars.refv.num) / (abs(gv.newpars.refv)+1)
if(max(diff.grads.newpars.refv) > 1/5) warning('gradient error - numerical and analytic gradients differ for the new pars ref v') 

diff.grads.refpars.refv <- abs(gv.refpars.refv - gv.refpars.refv.num) / (abs(gv.refpars.refv)+1)
if(max(diff.grads.refpars.refv) > 1/5) warning('gradient error - numerical and analytic gradients differ for the ref pars ref v') 

if (graphics){
	plot(gv.newpars.refv, gv.newpars.refv.num, main="gradients", xlab="analytic", ylab="numeric"); abline(0, 1, col=2)
	points(gv.newpars.newv, gv.newpars.newv.num, col=rainbow(length(parvec.new), end=0.7), pch=16, cex=0.7) # note, that numderiv is not exact - see also axis scales!
	points(gv.refpars.refv, gv.refpars.refv.num, col=3, pch=3) # note, that numderiv is not exact - see also axis scales!
	legend('topleft', legend=c('new pars ref v', 'new pars new v', 'ref pars ref v'), col=c(1,2,3), pch=c(1,16, 3), pt.cex=c(1,0.6, 1), bty="n")
	
	plot(gv.newpars.newv, gv.newpars.newv.num, col=rainbow(length(parvec.new), end=0.7), pch=16, cex=1, xlab="analytic", ylab="numeric", main="grads: new pars new v") # note, that numderiv is not exact - see also axis scales!

	matplot(cbind(diff.grads.newpars.refv, diff.grads.newpars.newv), t="p", pch=c(16,16), cex=c(1, 0.5), col=c(1,2), xlab="parameter", ylab="relative error in gradient", axes=F, main="relative error")
	points(diff.grads.refpars.refv, cex=1, pch=3, col=3)
	points(diff.grads.newpars.newv, cex=0.7, pch=16, col=rainbow(length(parvec.new), end=0.7))
	axis(1, 1:length(parvec.new), names(parvec.new), las=2, cex.axis=0.5); axis(2, las=2)
}
