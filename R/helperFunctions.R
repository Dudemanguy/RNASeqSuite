#fetch additional annotation from org R packages and insert into DGEList
#TODO: remove hardcoded options for species and annotations

annotationFetch <- function(count, species, symbol=FALSE) {
	selection <- paste("org.", species, ".eg.db", sep="")
	suppressPackageStartupMessages(require(selection, character.only=TRUE))
	if (selection == 'org.Mm.eg.db') {
		idfound <- rownames(count) %in% mappedRkeys(org.Mm.egREFSEQ)
		count <- count[idfound,]
		egREFSEQ <- toTable(org.Mm.egREFSEQ)
		if (symbol == TRUE) {
			egSYMBOL <- toTable(org.Mm.egSYMBOL)
		}
	}
	if (selection == 'org.Hs.eg.db') {
		idfound <- rownames(count) %in% mappedRkeys(org.Hs.egREFSEQ)
		count <- count[idfound,]
		egREFSEQ <- toTable(org.Hs.egREFSEQ)
		if (symbol == TRUE) {
			egSYMBOL <- toTable(org.Hs.egSYMBOL)
		}
	}
	if (species == 'org.Rn.eg.db') {
		idfound <- rownames(count) %in% mappedRkeys(org.Rn.egREFSEQ)
		count <- count[idfound,]
		egREFSEQ <- toTable(org.Rn.egREFSEQ)
		if (symbol == TRUE) {
			egSYMBOL <- toTable(org.Rn.egSYMBOL)
		}
	}
	m <- match(rownames(count), egREFSEQ$accession)
	count$genes$EntrezGene <- egREFSEQ$gene_id[m]
	if (symbol == TRUE) {
		m <- match(count$genes$EntrezGene, egSYMBOL$gene_id)
		count$genes$Symbol <- egSYMBOL$symbol[m]
	}
	d <- duplicated(count$genes$EntrezGene)
	count <- count[!d,]
	count
}

validDGEList <- function(y) {
#	Check for standard components of DGEList object
#	Gordon Smyth
#	20 Nov 2013

	if (is.null(y$counts)) {
		stop("No count matrix")
	}
	y$counts <- as.matrix(y$counts)
	nlib <- ncol(y$counts)
	if (is.null(y$samples$group)) {
		y$samples$group <- gl(1,nlib)
	}
	if (is.null(y$samples$lib.size)) {
		y$samples$lib.size <- colSums(y$counts)
	}
	if (is.null(y$samples$norm.factors)) {
		y$samples$norm.factors <- rep.int(1,nlib)
	}
	y
}

getOffset <- function(y) {
#	Extract offset vector or matrix from data object and optional arguments.
#	By default, offset is constructed from the lib.size and norm.factors
#	but offset supplied explicitly takes precedence

#	Gordon Smyth
#	26 Jan 2011. Last modified 11 Jan 2012.

	offset <- y$offset
	lib.size <- y$samples$lib.size
	norm.factors <- y$samples$norm.factors
	
	if (!is.null(offset)) {
		return(offset)
	} else {		
		if (!is.null(norm.factors)) {
			lib.size <- lib.size*norm.factors
		}
		return(log(lib.size))
	}
}

condLogLikDerSize <- function(y, r, der=1L) {
#	Derivatives of the conditional log-likelihood function (given the row sum)
#	with respect to r=1/dispersion
#	for a single group of replicate libraries, all of the same total size

#	Vector interpreted as matrix of one row, i.e., one gene
	if (is.vector(y)) {
		y <- matrix(y,nrow=1)
	} else {
		y <- as.matrix(y)
	}

	n <- ncol(y)
	m <- rowMeans(y)

	switch(der+1L,
		rowSums(lgamma(y+r)) + lgamma(n*r) - lgamma(n*(m+r)) - n*lgamma(r),
		rowSums(digamma(y+r)) + n*digamma(n*r) - n*digamma(n*(m+r)) - n*digamma(r),
		rowSums(trigamma(y+r)) + n^2*trigamma(n*r) - n^2*trigamma(n*(m+r)) - n*trigamma(r)
	)
}

commonCondLogLikDerDelta <- function(y, delta, der=0) {
# Calculates the common conditional log-likelihood (i.e. summed over all tags) - necessary so that optimize can be applied in estimateCommonDisp
# Davis McCarthy, July 2009

	l0 <- 0
	for (i in 1:length(y)) {
		l0 <- condLogLikDerDelta(y[[i]],delta,der=der)+l0
	}
	sum(l0)
}

condLogLikDerDelta <- function(y,delta,der=1L) {
# Derivatives of log-likelihood function wrt to delta
# r=1/dispersion and delta=1/(1+r)=dispersion/(1+dispersion)
# der is order of derivative required (0th deriv is the function)
# Written by Mark Robinson, edited by Davis McCarthy, February 2009

#	Vector interpreted as matrix of one row, i.e., one gene
	if (is.vector(y)) {
		y <- matrix(y,nrow=1)
	} else {
		y <- as.matrix(y)
	}
	if ( !(length(delta)==1 | length(delta)==nrow(y)) ) {
		stop("delta must be of length 1 or nrow(y)")
	}

	r <- (1/delta)-1
	switch(der+1L,
		condLogLikDerSize(y,r,der=0L),
		condLogLikDerSize(y,r,der=1L)*(-delta^(-2)),
		condLogLikDerSize(y,r,der=1L)*2*(delta^(-3))+condLogLikDerSize(y,r,der=2)*(delta^(-4))
	)
}

getDispersion <- function(y) {
#	Get most complex dispersion values from DGEList object
#	Gordon Smyth
#	Created 12 Dec 2011.  Last modified 3 Oct 2012.

	if (!is.null(y$tagwise.dispersion)) {
		dispersion <- y$tagwise.dispersion
		attr(dispersion,"type") <- "tagwise"
	} else {
		if (!is.null(y$trended.dispersion)) {
			dispersion <- y$trended.dispersion
			attr(dispersion,"type") <- "trended"
		} else {
			if (!is.null(y$common.dispersion)) {
				dispersion <- y$common.dispersion
				attr(dispersion,"type") <- "common"
			} else {
				dispersion <- NULL
			}
		}
	}
	dispersion
}

locfitByCol <- function(y, x=NULL, weights=1, span=0.5, degree=0) {
#	Gordon Smyth
#	20 Aug 2012.  Last modified 15 June 2016.

	y <- as.matrix(y)
	ntags <- nrow(y)
	weights <- rep_len(weights,ntags)
	if (is.null(x)) {
		x <- 1:ntags
	}
	if (span*ntags<2 || ntags<=1) {
		return(y)
	}
	for (j in 1:ncol(y)) {
		y[,j] <- fitted(locfit(y[,j]~x,weights=weights, alpha=span,deg=degree))
	}
	y
}

maximizeInterpolant <- function(x, y) {
# maximizeInterpolant: written by Aaron Lun
#
# This function takes an ordered set of spline points and a likelihood matrix where each row 
# corresponds to a tag and each column corresponds to a spline point. It then calculates the 
# position at which the maximum interpolated likelihood occurs for each by solving the derivative
# of the spline function.

    if (is.vector(y)) {
        y <- rbind(y)
        warning("coverting vector of likelihoods to matrix format for interpolation")
    }
    if (length(x)!=ncol(y)) { 
        stop("number of columns must equal number of spline points")
    } else if (is.unsorted(x) || anyDuplicated(x)) {
        stop("spline points must be unique and sorted")
    }

#	Performing some type checking.
	if (!is.double(x)) {
		storage.mode(x) <- "double"
	}
	if (!is.double(y)) {
		storage.mode(y) <- "double"
	}
    out <- .Call(.cxx_maximize_interpolant, x, y)
    return(out)
}
