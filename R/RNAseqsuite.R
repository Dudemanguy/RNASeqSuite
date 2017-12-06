#argumentinvalid
#general function to check argument input types and ensure they are valid

argumentinvalid <- function (dataframe, number, string1, string2) {
	if (is.null(dataframe)) {
		dataframe = FALSE
	}
	else if (class(dataframe) != "data.frame") {
		return(paste("Error,", deparse(substitute(dataframe)), "is not a dataframe"))
		return(TRUE)
	}

	if (is.null(number)) {
		number = FALSE
	}
	else if (class(number) != "numeric") {
		return(paste("Error,", deparse(substitute(number)), "is not numeric"))
		return(TRUE)

	}

	if (is.null(string1) & is.null(string2) | is.null(string2)) {
		string1 = FALSE
	}
	else if (class(string1) != "character") {
		return(paste("Error,", deparse(substitute(string1)), "is not a character"))
		return(TRUE)
	}
	else if (class(string2) != "character") {
		return(paste("Error,", deparse(substitute(string2)), "is not a character"))
		return(TRUE)
	}
}


#countmatrixsubset
#create subset of count matrix based upon entered group

#TODO: currently this function will only consider pairs. Add support for abitrary n?
countmatrixsubset <- function (x, y, group1, group2) {
	#if (inputinvalid(x, y)) {
	#	stop("Error, invalid input")
	#}
	#if (class(y) != "data.frame") {
	#	stop("Error, list of groups must be a data frame")
	#}
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
	#if (!((group1 %in% x[,2]) & (group2 %in% x[,2]))) {
	#	stop("Error, invalid group names")
	#}
	else {
		groupselection <- subset(x, V2 == group1 | V2 == group2)
		getgroup <- groupselection[,2]
		getgroup <- factor(getgroup, levels=unique(getgroup))
		return(list(getgroup, group1, group2))
	}
}


#CFilter
#Custom filter function based around standard deviation of normalized vectors

#TODO: Rewrite function to apply to abitrary n conditions?
#TODO: Finish porting to arbitrary dimension
CFilter <- function(x, y, z) {
	#if (class(x) != "data.frame") {
	#	stop("Error, please enter a data frame for x")
	#}
	if ((class(y) != "numeric") & (y < 1)) {
		stop("Error, please enter an integer greater than one for y")
	}
	else {
		x_filter <- x[rowSums(x) != 0,]
		x_avg_mat1 <- colMeans(x, dims=(1:nrow(subset(x, V2==grouping[[2]]))))
		x_avg_mat2 <- colMeans(x, dims=(nrow(subset(x, V2==grouping[[2]])):nrow(x)))
		x_filter[,1] <- x_filter[,1]/norm(data.matrix(x[,1]), type="f")
		x_filter[,2] <- x_filter[,2]/norm(data.matrix(x[,2]), type="f")
		x_diff <- x_filter[,1, drop=FALSE] - x_filter[,2, drop=FALSE]
		mean_x <- mean(x_diff[,1])
		std_x <- sd(x_diff[,1])
		x_compressed <- x_diff[(mean_x - (y*std_x)) <= x_filter[,1, drop=FALSE] & 
						x_filter[,1, drop=FALSE] <= (mean_x + (y*std_x)), 1, drop=FALSE]
		print(paste(nrow(x) - nrow(x_filter), "zero elements discarded"))
		print(paste(nrow(x_diff - nrow(x_compressed)), "outliers removed"))
		return(x_compressed)
	}
}

#countmatrix


#countmatrixFilter
#reads the supplied count matrix of reads and group data; filters data according to group

countmatrixFilter <- function(x, y, htsfilter, cfilter=z) {

	if(missing(htsfilter)) {
		htsfilter = TRUE
	}
	#if (HTSFilter != (TRUE | FALSE)) {
	#	stop("Error, HTSFilter is a boolean value.")
	#}

	#if(missing(Cfilter)) {
	#	Cfilter = FALSE
	#}
	#if (Cfilter != (TRUE | FALSE)) {
	#	stop("Error, CFilter is a boolean value.")
	#}

	else {
		#datagroup <- group(x, y[[2]], y[[3]])
		if (htsfilter == TRUE) {
			htsfilter <- HTSFilter(x, y[[1]], s.min=1, s.max=200, s.len=25)
			filter2 <- htsfilter$filteredData

			if (cfilter > 0) {
				countmatrixcfilter <- CFilter(filter2, 1, y[[1]])
				return(countmatrixcfilter)
			}
			#else {
			#	return(countmatrixhtsfilter)
			#}
		}
	}
}
