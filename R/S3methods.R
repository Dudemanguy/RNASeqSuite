as.matrix.DataList <- function(x, ...) {
	as.matrix(x$counts)
}

dim.DataList <- function(x) {
	if (is.null(x$counts)) {
		c(0,0) 
	}
	else {
		dim(as.matrix(x$counts))
	}
}

dimnames.DataList <- function(x) {
	dimnames(x$counts)
}

assign("dimnames<.DataList", function(x, value) {
	dimnames(x$counts) <- value
	if (!is.null(x$samples)) {
		row.names(x$samples) <- value[[2]]
	}
	if (!is.null(x$genes)) {
		row.names(x$genes) <- value[[1]]
	}
	x
})
