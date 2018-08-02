

scalebar <- function(x1, x2, y1, y2, xscale=NULL, yscale=NULL){
	lines(c(x1, x1, x2), c(y2, y1, y1))
	if (!is.null(xscale)) text((x1 + x2)/2, y1, xscale, pos=1)
	if (!is.null(yscale)) text(x1, (y1 + y2)/2, yscale, pos=4)
	# text(125, -69, "2 mV", pos=4)
}
