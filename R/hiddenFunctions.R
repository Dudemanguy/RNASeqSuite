#general function to check argument input types and ensure they are valid
.argumentValid <- function(classnames, refnames) {
	classframe <- data.frame(names(classnames[]))
	classlist <- list()
	for (i in seq_along(classnames)) {
		classlist[[i]] <- class(classnames[[i]])
	}
	classvector <- unlist(classlist)
	classframe["class"] <- classvector
	classframe["ref"] <- refnames
	for (i in nrow(classframe)) {
		if (classframe[i, 2] != classframe[i, 3]) {
			stop(paste(classframe[i, 2], " is an invalid class for ", classframe[i, 1], ". It must be a ", classframe[i, 3], sep=""), call.=FALSE)
		}
	}
}


#split by group and return list
.ctSplit <- function(data, group)
UseMethod(".ctSplit")

.ctSplit.DGEList <- function(data) {
	group <- data$samples$group
	ctSplit(data$counts, group=group)
}

.ctSplit.default <- function(data, group) {
	if (class(group) == "list") {
		selection <- group$frame
	}
	if (class(group) == "factor") {
		selection <- data.frame(group=group)
		rownames(selection) <- colnames(data)
	}
	columns <- rownames(selection)
	matframe <- subset(data, select=eval(parse(text=list(columns))))
	subframe <- data.frame()
	sublist <- list()
	j <- 1
	k <- 0
	for (i in 1:(nrow(selection))) {
		if (i != nrow(selection)) {
			if (selection[i, 1] == selection[i+1, 1]) {
				j <- j+1
				k <- k+1
			}
			if (!(selection[i, 1] == selection[i+1, 1])) {
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


#helper function for normalizing vectors
.norm_vector <- function(vec, norm_type) {
	normalized <- vec/norm(data.matrix(vec), type=norm_type)
	normalized
}


#abstract operation to remove rows with rowsums equal to zeros
.removeZeros <- function(df) {
	df_nozero <- df[rowSums(df) != 0,]
}


#abstract selection operations for apply
.select <- function(sub, whole) {
	output <- whole[sub,, drop=FALSE]
	output
}


#general function to match strings across dataframes
.stringMatch <- function(data, index, strings) {
	for (i in 1:length(strings)) {
		if (!(strings[[i]] %in% data[,index])) {
			return(FALSE)
		}
	}
	return(TRUE)
}
