#argumentvalid
#general function to check argument input types and ensure they are valid

#TODO: Generalize arguments instead of explicitly stating them
#TODO: Make error messages cleaner
argumentvalid <- function (x) {
	type <- c("data.frame", "numeric", "character", "logical")
	possiblearguments <- data.frame(type)
	if (x %in% possiblearguments[,1]) {
			truthentry <- setNames(list(TRUE), x)
	}
	else if (!(x %in% possiblearguments[,1])) {
			truthentry <- setNames(list(FALSE), x)
	}
	return(truthentry)
}

argumentvalidreturn <- function (...) {
	argumentlist <-	lapply(c(...), argumentvalid)
	argumentframe <- data.frame(argumentlist)
	return(argumentframe)
}

errormessages <- function (x) {
	if (!(is.null(x$data.frame))) {
		if (x$data.frame == FALSE) {
			return(message(paste("Error, function does not accept objects with class 'data.frame'.")))
		}
	}
	if (!(is.null(x$numeric))) {
		if (x$numeric == FALSE) {
			return(message(paste("Error, function does not accept objects with class 'numeric'.")))
		}
	}
	if (!(is.null(x$character))) {
		if (x$character == FALSE) {
			return(message(paste("Error, function does not accept objects with class 'character'.")))
		}	
	}
	if (!(is.null(x$logical))) {
		if (x$logical == FALSE) {
			return(message(paste("Error, function does not accept objects with class 'logical'.")))
		}
	}
	if (!(is.null(x$function.))) {
		if (x$function. == FALSE) {
			return(message(paste("Error, function does not accept objects with class 'function'.")))
		}
	}
	else {
		return(TRUE)
	}
}


#stringmatch
#general function to match strings across dataframes

stringmatch <- function (dataframe, index, string1, string2) {
	if (!((string1 %in% dataframe[,index]) & string2 %in% dataframe[,index])) {
		return(FALSE)
	}
	else {
		return(TRUE)
	}
}


#countmatrixsubset
#create subset of count matrix based upon entered group

#TODO: currently this function will only consider pairs. Add support for abitrary n?
countmatrixsubset <- function (dataframe, groupframe, group1, group2) {
	argumentcheck <- c(class(dataframe),class(groupframe),class(group1),class(group2))
	if (!(errormessages(argumentvalidreturn(argumentcheck)))) {
		stop()
	}
	#if (!(stringmatch(groupframe, 2, group1, group2))) {
	#	stop(paste("Error, no entries in", deparse(substitute(groupframe)), "match the arguments."))
	#}
#	else {
#		groupselection <- subset(groupframe, V2 == group1 | V2 == group2)
#		columnstokeep <- as.vector(groupselection[,1])
#		matrixselection <- subset(dataframe, select = eval(parse(text=list(columnstokeep))))
#		return(matrixselection)
#	}
}


#group
#obtains the group from a supplied tab-delimited list

#TODO: only supports two groups. Add support for abitrary n?
grouping <- function(groupframe, group1, group2) {
	if (!(argumentvalid(groupframe, NULL, NULL, group1, group2, NULL, NULL))) {
		stop()
	}
	if (!(stringmatch(groupframe, 2, group1, group2))) {
		stop(paste("Error, no entries in", deparse(substitute(groupframe)), 
					"match the arguments."))
	}
	else {
		groupselection <- subset(groupframe, V2 == group1 | V2 == group2)
		getgroup <- groupselection[,2]
		getgroup <- factor(getgroup, levels=unique(getgroup))
		return(list(getgroup, group1, group2))
	}
}


#CFilter
#Custom filter function based around standard deviation of normalized vectors

#TODO: Rewrite function to apply to abitrary n conditions?
#TODO: Finish porting to arbitrary dimension
CFilter <- function(df, sd, groupfactor) {
	#if (class(x) != "data.frame") {
	#	stop("Error, please enter a data frame for x")
	#}
	if ((class(sd) != "numeric") & (sd < 1)) {
		stop("Error, please enter an integer greater than one for sd")
	}
	else {
		df_filter <- df[rowSums(df) != 0,]
		group1width <- table(groupfactor[[1]])[[1]]
		group2width <- table(groupfactor[[1]])[[2]]
		df_group1 <- df_filter[(1:group1width)]
		df_group2 <- df_filter[((group1width+1):(group1width+group2width))]
		df_avg1 <- data.matrix(rowMeans(df_group1))
		df_avg2 <- data.matrix(rowMeans(df_group2))
		df_norm1 <- data.matrix(df_avg1/norm(df_avg2, type="f"))
		df_norm2 <- data.matrix(df_avg2/norm(df_avg2, type="f"))
		df_diff <- data.matrix(df_norm1 - df_norm2)
		mean_df <- mean(df_diff)
		sd_df <- sd(df_diff)
		df_compressed <- df_diff[(mean_df - (sd*sd_df)) <= df_norm1 &
						df_norm2 <= (mean_df + (sd*sd_df)), 1, drop=FALSE]
		print(paste(nrow(df) - nrow(df_filter), "zero elements discarded"))
		print(paste(nrow(df_diff) - nrow(df_compressed), "outliers removed"))
		return(df_compressed)
	}
}


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
