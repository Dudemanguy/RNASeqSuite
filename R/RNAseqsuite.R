#countmatrixsubset
#create subset of count matrix based upon entered group

#TODO: currently this function will only consider pairs. Add support for abitrary n?
countmatrixsubset <- function (x, y, group1, group2) {
	if (class(x) != "data.frame") {
		stop("Error, count matrix must be a data frame")
	}
	if (class(y) != "data.frame") {
		stop("Error, list of groups must be a data frame")
	}
	if (!((group1 %in% y[,2]) & (group2 %in% y[,2]))) {
		stop("Error, invalid group names")
	}
	else {
		groupselection <- subset(y, V2 == group1 | V2 == group2)
		columnstokeep <- as.vector(groupselection[,1])
		matrixselection <- subset(x, select = eval(parse(text=list(columnstokeep))))
		return(matrixselection)
	}
}

#group
#obtains the group from a supplied tab-delimited list

#TODO: only supports two groups. Add support for abitrary n?
group <- function(x, group1, group2) {
	if (class(x) != "data.frame") {
		stop("Error, x must be a data frame")
	}
	if (!((group1 %in% x[,2]) & (group2 %in% x[,2]))) {
		stop("Error, invalid group names")
	}
	else {
		groupselection <- subset(x, V2 == group1 | V2 == group2)
		getgroup <- groupselection[,2]
		getgroup <- factor(getgroup, levels=unique(getgroup))
		return(list(getgroup, group1, group2))
	}
}


#CFilter
#Custom filter function based around standard deviation of normalized vectors

#TODO: Rewrite function to apply to matrix of arbitrary n
CFilter <- function(x, y, group) {
	if (class(x) != "data.frame") {
		stop("Error, please enter a data frame for x")
	}
	if ((class(y) != "numeric") & (y < 1)) {
		stop("Error, please enter an integer greater than one for y")
	}
	else {
		x_filter <- x[rowSums(x) != 0,]
		x_filter[,1] <- x_filter[,1]/norm(data.matrix(x[,1]), type="f")
		x_filter[,2] <- x_filter[,2]/norm(data.matrix(x[,2]), type="f")
		x_diff <- x_filter[,1, drop=FALSE] - x_filter[,2, drop=FALSE]
		mean_x <- mean(x_diff[,1])
		std_x <- sd(x_diff[,1])
		x_compressed <- x_diff[(mean_x - (y*std_x)) <= x_filter[,1, drop=FALSE] & 
				x_filter[,1, drop=FALSE] <= (mean_x + (y*std_x)), 1, drop=FALSE]
		print(paste(nrow(x) - nrow(x_filter), "zero elements discarded"))
		print(paste(nrow(x_diff - nrow(x_compressed), "outliers removed")))
		return(x_compressed)
	}
}

#countmatrix


#countmatrixFilter
#reads the supplied count matrix of reads and group data; filters data according to group

countmatrixFilter <- function(x, y, HTSFilter, Cfilter) {

	if(missing(HTSFilter)) {
		HTSFilter = TRUE
	}
	if (HTSFilter != (TRUE | FALSE)) {
		stop("Error, HTSFilter is a boolean value.")
	}

	if(missing(Cfilter)) {
		Cfilter = FALSE
	}
	if (Cfilter != (TRUE | FALSE)) {
		stop("Error, CFilter is a boolean value.")
	}

	else {
		datagroup <- group(x, y[[2]], y[[3]])
		if (HTSFilter == TRUE) {
			htsfilter <- HTSFilter(countmatrixsubset, group, s.min=1, s.max=200, s.len=25)
			countmatrixhtsfilter <- filter$filteredData

			if (Cfilter == TRUE) {
				countmatrixcfilter <- Cfilter(x, 1, group)
				return(countmatrixcfilter)
			}
			else {
				return(countmatrixhtsfilter)
			}
		}
	}
}
