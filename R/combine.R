cbind.DataList <- function(..., deparse.level=1) {
#  Combine samples from DataList objects with same genelists
#  Gordon Smyth
#  20 June 2017

	objects <- list(...)
	nobjects <- length(objects)
	out <- objects[[1]]
	out$design <- NULL
	if (nobjects > 1) {
		for (i in 2:nobjects) {
			out$counts <- cbind(out$counts,objects[[i]]$counts)
			out$weights <- cbind(out$weights,objects[[i]]$weights)
			out$offset <- cbind(out$offset,objects[[i]]$offset)
			out$samples <- rbind(out$samples,objects[[i]]$samples)
		}
	}
	out
}

rbind.DataList <- function(..., deparse.level=1) {
#  Combine genes from DataList objects with same samples
#  Gordon Smyth
#  20 June 2017

	objects <- list(...)
	nobjects <- length(objects)
	out <- objects[[1]]
	if (nobjects > 1) {
		for (i in 2:nobjects) {
			out$counts <- rbind(out$counts,objects[[i]]$counts)
			out$weights <- rbind(out$weights,objects[[i]]$weights)
			out$offset <- rbind(out$offset,objects[[i]]$offset)
			out$genes <- rbind(out$genes,objects[[i]]$genes)
		}
	}
	out
}
