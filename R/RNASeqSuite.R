#argumentValid
#general function to check argument input types and ensure they are valid

#TODO: Generalize arguments instead of explicitly stating them
#TODO: Make valid messages cleaner
argumentValid <- function (x) {
	type <- c("data.frame", "numeric", "character", "logical")
	possibleArguments <- data.frame(type)
	if (x %in% possibleArguments[,1]) {
			truthEntry <- setNames(list(TRUE), x)
	}
	else if (!(x %in% possibleArguments[,1])) {
			truthEntry <- setNames(list(FALSE), x)
	}
	return(truthEntry)
}

argumentValidreturn <- function (...) {
	argumentList <-	lapply(c(...), argumentValid)
	argumentFrame <- data.frame(argumentList)
	return(argumentFrame)
}

valid <- function (x) {
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


#stringMatch
#general function to match strings across dataframes

stringMatch <- function (data, index, string1, string2) {
	if (!((string1 %in% data[,index]) & string2 %in% data[,index])) {
		return(FALSE)
	}
	else {
		return(TRUE)
	}
}


#countSelection
#create subset of count matrix based upon entered group

#TODO: currently this function will only consider pairs. Add support for abitrary n?
countSelection <- function (data, frame, group) {
	check <- c(class(data),class(frame),class(group[[2]]),class(group[[3]]))
	if (!(valid(argumentValidreturn(check)))) {
		stop()
	}
	if (!(stringMatch(frame, 2, group[[2]], group[[3]]))) {
		stop(paste("Error, no entries in", deparse(substitute(groupframe)), "match the arguments."))
	}
	else {
		selection <- subset(frame, V2 == group[[2]] | V2 == group[[3]])
		columns <- as.vector(selection[,1])
		matframe <- subset(data, select = eval(parse(text=list(columns))))
		return(matframe)
	}
}


#group
#obtains the group from a supplied tab-delimited list

#TODO: only supports two groups. Add support for abitrary n?
grouplist <- function(frame, group1, group2) {
	check <- c(class(frame), class(group1), class(group2))
	if (!(valid(argumentValidreturn(check)))) {
		stop()
	}
	if (!(stringMatch(frame, 2, group1, group2))) {
		stop(paste("Error, no entries in", deparse(substitute(frame)), 
					"match the arguments."))
	}
	else {
		selection <- subset(frame, V2 == group1 | V2 == group2)
		getgroup <- selection[,2]
		getgroup <- factor(getgroup, levels=unique(getgroup))
		return(list(getgroup, group1, group2))
	}
}


#cFilter
#Custom filter function based around standard deviation of normalized vectors

#TODO: Rewrite function to apply to abitrary n conditions?
#TODO: Finish porting to arbitrary dimension
cFilter <- function(df, sd, group) {
	check <- c(class(df), class(sd), class(group))
	if (!(valid(argumentValidreturn(check)))) {
		stop()
	}
	if (sd < 0) {
		stop(paste("Error,", deparse(substitute(sd)), "must be positive."))
	}
	else {
		df_filter <- df[rowSums(df) != 0,]
		width1 <- table(group[[1]])[[1]]
		width2 <- table(group[[1]])[[2]]
		df_group1 <- data.matrix(df_filter[,(1:width1)])
		df_group2 <- data.matrix(df_filter[,((width1+1):(width1+width2))])
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


#countFilter
#reads the supplied count matrix of reads and group data; filters data according to group

countFilter <- function(data, frame, group, htsfilter, cfilter) {

	if(missing(htsfilter)) {
		htsfilter = TRUE
	}
	if(missing(cfilter)) {
		cfilter = 0
	}
	check <- c(class(data),class(frame),class(group),class(htsfilter),class(cfilter))
	if (!(valid(argumentValidreturn(check)))) {
		stop()
	}

	else {
		count <- countSelection(data, frame, group)
		if (htsfilter == TRUE) {
			htsfilter <- HTSFilter(count, group[[1]], s.min=1, s.max=200, s.len=25)
			htsfiltered <- htsfilter$filteredData

			if (cfilter > 0) {
				countcfilter <- CFilter(htsfiltered, cfilter, group[[1]])
				allfilter <- count[rownames(count) %in% rownames(countcfilter),]
				return(allfilter)
			}
			else {
				return(htsfiltered)
			}
		}
	}
}
