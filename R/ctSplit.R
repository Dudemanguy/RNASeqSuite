#split by group and return list

ctSplit <- function(data, group)
UseMethod("ctSplit")

ctSplit.DataList <- function(data) {
	group <- data$samples$group
	ctSplit(data$counts, group=group)
}

ctSplit.default <- function(data, group) {
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
