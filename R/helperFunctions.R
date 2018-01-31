#general functions to check argument input types and ensure they are valid
#TODO: Make valid messages cleaner

.argumentValid <- function (x, y) {
	for (i in seq_along(x)) {
		if (is(x[[i]], "data.frame") & (!(y[[i]] == "data.frame"))) {
			stop(paste(deparse(substitute(x[[i]])), "is not a data frame."))
		}
		if (is(x[[i]], "numeric") & (!(y[[i]] == "numeric"))) {
			stop(paste(deparse(substitute(x[[i]])), "is not numeric."))
		}
		if (is(x[[i]], "character") & (!(y[[i]] == "character"))) {
			stop(paste(deparse(substitute(x[[i]])), "is not a character vector."))
		}
		if (is(x[[i]], "logical") & (!(y[[i]] == "logical"))) {
			stop(paste(deparse(substitute(x[[i]])), "is not a boolean value."))
		}
		if (is(x[[i]], "list") & (!(y[[i]] == "list"))) {
			stop(paste(deparse(substitute(x[[i]])), "is not a list."))
		}
	}
}

#general function to match strings across dataframes

.stringMatch <- function (data, index, strings) {
	for (i in 1:length(strings)) {
		if (!(strings[[i]] %in% data[,index])) {
			return(FALSE)
		}
	}
	return(TRUE)
}

#helper function for normalizing vectors

.norm_vector <- function(vec, norm_type) {
	normalized <- vec/norm(data.matrix(vec), type=norm_type)
	return(normalized)
}

#abstract selection operations for apply

.select <- function(x, y) {
	output <- y[x,]
	return(output)
}

#grouplist <- function(frame, groupselect) {
#	check <- list(frame, groupselect)
#	ref <- list("data.frame","list")
#	.argumentValid(check, ref)
#	if (!(.stringMatch(frame, 2, groupselect))) {
#		stop(paste("Error, some entries in", deparse(substitute(frame)), 
#					"do not match the arguments."))
#	}
#	else {
#		selection_grep <- lapply(groupselect, '==', frame[,2])
#		selection <- data.frame()
#		for (i in 1:length(selection_grep)) {
#			selection <- rbind(selection, frame[selection_grep[[i]],])
#		}
#		getgroup <- selection[,2]
#		getgroup <- factor(getgroup, levels=unique(getgroup))
#		groupselect[["factors"]] <- getgroup
#		return(groupselect)
#	}
#}
