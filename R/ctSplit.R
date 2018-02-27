#split by group and return list

ctSplit <- function(data, group)
UseMethod("ctSplit")

ctSplit.DataList <- function(data, group) {
	ctSplit(data$counts, group=group)
}

ctSplit.default <- function(data, group) {
	selection <- group$frame
	columns <- as.vector(selection[,1])
	matframe <- subset(data, select=eval(parse(text=list(columns))))
	subframe <- data.frame()
	sublist <- list()
	j <- 1
	k <- 0
	for (i in 1:(nrow(selection))) {
		if (i != nrow(selection)) {
			if (selection[i, 2] == selection[i+1, 2]) {
				j <- j+1
				k <- k+1
			}
			if (!(selection[i, 2] == selection[i+1, 2])) {
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
