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

.select <- function(sub, whole) {
	output <- whole[sub,]
	return(output)
}

#split by group and seperate dataframe into lists

.ctSplit <- function (data, frame, group) {
	selection_grep <- sapply(levels(group), '==', frame[,2])
	selection_list <- apply(selection_grep, 2, .select, frame)
	selection <- do.call("rbind", selection_list)
	columns <- as.vector(selection[,1])
	matframe <- subset(data, select = eval(parse(text=list(columns))))
	subframe <- data.frame()
	sublist <- list()
	j <- 1
	k <- 0
	for (i in 1:(nrow(selection))) {
		if (i != nrow(selection)) {
			if (selection[i,2] == selection[i+1,2]) {
				j <- j+1
				k <- k+1
			}
			if (!(selection[i,2] == selection[i+1,2])) {
				subframe <- matframe[,(j-k):j]
				sublist[[i]] <- subframe
				print(j)
				print(k)
				j <- j+1
				k <- 0
			}
		}
		if (i == nrow(selection)) {
			print(j)
			print(k)
			subframe <- matframe[,(j-k):j]
			sublist[[i]] <- subframe
		}
	}
	return(sublist)
}
