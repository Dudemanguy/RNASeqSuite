#general functions to check argument input types and ensure they are valid
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

validReturn <- function (...) {
	argumentList <-	lapply(c(...), argumentValid)
	x <- data.frame(argumentList)

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


#general function to match strings across dataframes

stringMatch <- function (data, index, string1, string2) {
	if (!((string1 %in% data[,index]) & string2 %in% data[,index])) {
		return(FALSE)
	}
	else {
		return(TRUE)
	}
}


#create subset of ct matrix based upon entered group
#TODO: currently this function will only consider pairs. Add support for abitrary n?

ctSelection <- function (data, frame, group) {
	check <- c(class(data),class(frame),class(group[[2]]),class(group[[3]]))
	if (!(validReturn(check))) {
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


#obtains the group from a supplied tab-delimited list
#TODO: only supports two groups. Add support for abitrary n?

grouplist <- function(frame, group1, group2) {
	check <- c(class(frame), class(group1), class(group2))
	if (!(validReturn(check))) {
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


#Custom filter function based around standard deviation of normalized vectors
#TODO: Rewrite function to apply to abitrary n conditions?

cFilter <- function(df, sd, group) {
	check <- c(class(df), class(sd), class(group))
	if (!(validReturn(check))) {
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


#reads the supplied ct matrix of reads and group data; filters data according to group

ctFilter <- function(data, frame, group, htsfilter, cfilter) {

	if(missing(htsfilter)) {
		htsfilter = TRUE
	}
	if(missing(cfilter)) {
		cfilter = 0
	}
	check <- c(class(data),class(frame),class(group),class(htsfilter),class(cfilter))
	if (!(validReturn(check))) {
		stop()
	}

	else {
		ct <- ctSelection(data, frame, group)
		if (htsfilter == TRUE) {
			htsfilter <- HTSFilter(ct, group[[1]], s.min=1, s.max=200, s.len=25)
			htsfiltered <- htsfilter$filteredData

			if (cfilter > 0) {
				ctcfilter <- cFilter(htsfiltered, cfilter, group[[1]])
				allfilter <- ct[rownames(ct) %in% rownames(ctcfilter),]
				return(allfilter)
			}
			else {
				return(htsfiltered)
			}
		}
	}
}

#uses edgeR to compute an exact test and find differentially expressed genes

edgeR <- function (data, frame, group, htsfilter, cfilter) {

	if(missing(htsfilter)) {
		htsfilter = TRUE
	}
	if(missing(cfilter)) {
		cfilter = 0
	}
	check <- c(class(data),class(frame),class(group),class(htsfilter),class(cfilter))
	if (!(validReturn(check))) {
		stop()
	}

	else {
		ct <- ctSelection(data, frame, group)	
		ct <- ctFilter(data, frame, group, htsfilter, cfilter)
		y <- DGEList(counts=ct, group=group[[1]])
		y <- calcNormFactors(y)
		y <- estimateDisp(y)
		et <- exactTest(y) 
		et_raw <- topTags(et, n=Inf, sort.by="none")
		et_frame <-  et_raw[[1]]
		width <- table(group[[1]])[[1]]
		a <- ct[,1:width]
		b <- ct[,(width+1):ncol(ct)]
		c <- data.frame(rowMeans(a))
		d <- data.frame(rowMeans(b))
		et_frame["Avg Ct A"] <- c
		et_frame["Avg Ct B"] <- d
		et_frame <- et_frame[,c(5,6,1,2,3,4)]
		et_frame <- et_frame[order(et_frame$FDR, decreasing=FALSE),]
		return(et_frame)
	}
}

#use DESeq2 to compute a Wald test and find differentially expressed genes
DESeq2 <- function (data, frame, group, htsfilter, cfilter) {

	if(missing(htsfilter)) {
		htsfilter = TRUE
	}
	if(missing(cfilter)) {
		cfilter = 0
	}
	check <- c(class(data),class(frame),class(group),class(htsfilter),class(cfilter))
	if (!(validReturn(check))) {
		stop()
	}
	
	else {
		ct <- ctSelection(data, frame, group)
		ct <- ctFilter(data, frame, group, htsfilter, cfilter)
		groupframe <- data.frame(group[[1]])
		colnames(groupframe) <- c("groupframe")
		dds <- DESeqDataSetFromMatrix(countData=ct, colData=groupframe, design=~groupframe)
		dds <- DESeq(dds)
		res <- data.frame(results(dds))
		width <- table(group[[1]])[[1]]
		a <- ct[,1:width]
		b <- ct[,(width+1):ncol(ct)]
		c <- data.frame(rowMeans(a))
		d <- data.frame(rowMeans(b))
		res["Avg Ct A"] <- c
		res["Avg Ct B"] <- d
		resOrder <- res[,c(7,8,1,2,3,4,5,6)]
		resOrder <- resOrder[order(resOrder$padj, decreasing=FALSE),]
		return(resOrder)
	}
}
