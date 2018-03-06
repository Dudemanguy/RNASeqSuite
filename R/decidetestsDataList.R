#  decideTestsData.R

decideTests.DataList <- function(object,adjust.method="BH",p.value=0.05,lfc=0,...) {

	decideTestsData(object=object,adjust.method=adjust.method,p.value=p.value,lfc=lfc)
}

decideTestsData <- function(object,adjust.method="BH",p.value=0.05,lfc=0) {
#	Accept or reject hypothesis tests across genes and contrasts
#	edgeR team. Original author was Davis McCarthy.
#	Created 15 August 2010. Last modified 14 Dec 2017.

#	Check object class
	if (!(is(object,"DataList"))) {
		stop("Need DataList object")
	}

#	Apply multiple testing
	p <- object$et_results$FDR
	isDE <- as.integer(p < p.value)

#	Extract logFC
	logFC <- object$et_results$logFC

#	Check for F-test with multiple logFC columns
	FTest <- is.null(logFC)

#	With single contrast, apply directionality and lfc threshold
	isDE[isDE & logFC<0] <- -1L
	SmallFC <- (abs(logFC) < lfc)
	isDE[SmallFC] <- 0L

#	Assemble TestResults object
	isDE <- matrix(isDE, ncol=1)
	row.names(isDE) <- row.names(object)
	colnames(isDE) <- paste(object$comparison,collapse="+")

#	Record possible values
	if (FTest) {
		attr(isDE,"levels") <- c(0L,1L)
		attr(isDE,"labels") <- c("NotSig","Sig")
	} 
	else {
		attr(isDE,"levels") <- c(-1L,0L,1L)
		attr(isDE,"labels") <- c("Down","NotSig","Up")
	}		

	new("TestResults", isDE)
}
