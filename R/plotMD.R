plotMD.DataList <- function(object, column=1, xlab=NULL, ylab=NULL, main=NULL, status=object$genes$Status, zero.weights=FALSE, prior.count=3, values=names(table(status)), p.value=0.05, ...) {
#	Mean-difference plot with color coding for controls
#	Gordon Smyth
#	Created 24 June 2015. Last modified 24 June 2015.

	if (is.null(object$et_results)) {
		nlib <- ncol(object)
		if (nlib < 2L) {
			stop("Need at least two columns")
		}

		#Convert column to integer if not already
		j <- 1:nlib
		names(j) <- colnames(object)
		column <- j[column[1]]

		logCPM <- cpm(object, log=TRUE, prior.count=prior.count)
		AveOfOthers <- rowMeans(logCPM[,-column,drop=FALSE],na.rm=TRUE)
		Diff <- logCPM[,column]-AveOfOthers
		Mean <- (logCPM[,column]+AveOfOthers)/2

		if (!zero.weights && !is.null(object$weights)) {
			w <- as.matrix(object$weights)[,column]
			Diff[ is.na(w) | (w <= 0) ] <- NA
		}
		xlab <- "Average log CPM (this sample and others)"
		ylab <- "log-ratio (this sample vs others)"
		main <- colnames(object)[column]

		plotWithHighlights(x=Mean,y=Diff,xlab=xlab,ylab=ylab,main=main,status=status,...)
	}

	if (!(is.null(object$et_results))) {

		status <- decideTestsData(object, adjust.method=adjust.method, p.value=p.value)
		status <- c("Down", "non-DE", "Up")[status+2L]
		values <- c("Up","Down")
		col <- c("red","blue")

		xlab <- "Average log CPM"
		ylab <- "log-fold-change"
		main <- paste(object$comparison[2],"vs",object$comparison[1])

		plotWithHighlights(x=object$et_results$logCPM,y=object$et_results$logFC,xlab=xlab,ylab=ylab,main=main,status=status,values=values,col=col,...)
	}
}
