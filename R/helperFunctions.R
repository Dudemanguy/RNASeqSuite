#general functions to check argument input types and ensure they are valid

.argumentValid <- function (classnames, refnames) {
	classframe <- data.frame(names(classnames[]))
	classlist <- list()
	for (i in seq_along(classnames)) {
		classlist[[i]] <- class(classnames[[i]])
	}
	classvector <- unlist(classlist)
	classframe["class"] <- classvector
	classframe["ref"] <- refnames
	for (i in nrow(classframe)) {
		if (classframe[i,2] != classframe[i,3]) {
			stop(paste(classframe[i,2], " is an invalid class for ", classframe[i,1], ". It must be a ", classframe[i,3], sep=""), call.=FALSE)
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

#abstract operation to remove rows with rowsums equal to zeros

.removeZeros <- function(df) {
	df_nozero <- df[rowSums(df) != 0,]
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
				j <- j+1
				k <- 0
			}
		}
		if (i == nrow(selection)) {
			subframe <- matframe[,(j-k):j]
			sublist[[i]] <- subframe
		}
	}
	sublist[sapply(sublist, is.null)] <- NULL
	return(sublist)
}
