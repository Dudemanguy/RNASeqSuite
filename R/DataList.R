DataList <- function(counts=matrix(0,0,0), lib.size=colSums(counts), norm.factors=rep(1,ncol(counts)), samples=NULL, group=NULL, genes=NULL, remove.zeros=FALSE) {

	#Check counts
	counts <- as.matrix(counts)
	nlib <- ncol(counts)
	ntags <- nrow(counts)
	if (nlib>0L && is.null(colnames(counts))) {
		colnames(counts) <- paste0("Sample",1L:nlib)
	}
	if (ntags>0L && is.null(rownames(counts))) {
		rownames(counts) <- 1L:ntags
	}

	#Check lib.size
	if (nlib != length(lib.size)) {
		stop("length of 'lib.size' must equal number of columns in 'counts'")
	}
	if (any(lib.size==0L)) {
		warning("library size of zero detected")
	}

	#Check norm.factors
	if (nlib != length(norm.factors)) {
		stop("Length of 'norm.factors' must equal number of columns in 'counts'")
	}

	#Check samples
	if (!is.null(samples)) {
		samples <- as.data.frame(samples)
		if (nlib != nrow(samples)) {
			stop("Number of rows in 'samples' must equal number of columns in 'counts'")
		}
	}
	
	#Check group
	if (is.null(group)) {
		if (!is.null(samples))
			group <- samples[,1]
		else
			group <- rep(1,ncol(counts))
	}
	group <- .dropEmptyLevels(group)
	if (nlib != length(group)) {
		stop("Length of 'group' must equal number of columns in 'counts'")
	}

	#Make data frame of sample information
	sam <- data.frame(group=group, lib.size=lib.size, norm.factors=norm.factors)
	if (!is.null(samples)) {
		sam <- data.frame(sam, samples)
	}
	samples <- sam
	if (anyDuplicated(colnames(counts))) {
		message("Repeated column names found in count matrix")
		row.names(samples) <- 1L:nlib
	} 
	else {
		row.names(samples) <- colnames(counts)
	}

	#Make object
	x <- new("DataList", list(counts=counts,samples=samples))

	#Add data frame of gene information
	if (!is.null(genes)) {
		genes <- as.data.frame(genes, stringsAsFactors=FALSE)
		if (nrow(genes) != ntags) {
			stop("Counts and genes have different numbers of rows")
		}
		x$genes <- genes
	}

	x
}

.dropEmptyLevels <- function(x) {
	if(is.factor(x)) {
		i <- which(tabulate(as.integer(x))>0L)
		if (length(i) < nlevels(x)) {
			x <- factor(x, levels=levels(x)[i])
		}
		return(x)
	}
		else {
			return(factor(x))
	}
}
