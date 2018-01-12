#general functions to check argument input types and ensure they are valid
#TODO: Make valid messages cleaner

argumentValid <- function (x, y) {
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

stringMatch <- function (data, index, strings) {
	for (i in 1:length(strings)) {
		if (!(strings[[i]] %in% data[,index])) {
			return(FALSE)
		}
	}
	return(TRUE)
}


#create subset of ct matrix based upon entered group

ctSelection <- function (data, frame, group) {
	check <- list(data, frame, group)
	ref <- list("data.frame","data.frame","list")
	argumentValid(check, ref)
	if (!(stringMatch(frame, 2, group[["factors"]]))) {
		stop(paste("Error, no entries in", deparse(substitute(groupframe)), "match the arguments."))
	}
	else {
		selection_grep <- lapply(group[1:(length(group)-1)], '==', frame$V2)
		selection <- data.frame(V1=character(), V2=character())
		for (i in 1:length(selection_grep)) {
			selection <- rbind(selection, frame[selection_grep[[i]],])
		}
		columns <- as.vector(selection[,1])
		matframe <- subset(data, select = eval(parse(text=list(columns))))
		return(matframe)
	}
}


#obtains the group from a supplied tab-delimited list

grouplist <- function(frame, groupselect) {
	check <- list(frame, groupselect)
	ref <- list("data.frame","list")
	argumentValid(check, ref)
	if (!(stringMatch(frame, 2, groupselect))) {
		stop(paste("Error, some entries in", deparse(substitute(frame)), 
					"do not match the arguments."))
	}
	else {
		selection_grep <- lapply(groupselect, '==', frame$V2)
		selection <- data.frame(V1=character(), V2=character())
		for (i in 1:length(selection_grep)) {
			selection <- rbind(selection, frame[selection_grep[[i]],])
		}
		getgroup <- selection[,2]
		getgroup <- factor(getgroup, levels=unique(getgroup))
		groupselect[["factors"]] <- getgroup
		return(groupselect)
	}
}


#Custom filter function based around standard deviation of normalized vectors
#TODO: Rewrite function to apply to abitrary n conditions?

cFilter <- function(df, sd, group) {
	check <- list(df, sd, group)
	ref <- list("data.frame", "numeric", "list")
	argumentValid(check, ref)
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
		df_norm1 <- data.matrix(df_avg1/norm(df_avg1, type="f"))
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
	else {
		check <- list(data, frame, group, htsfilter, cfilter)
		ref <- list("data.frame","data.frame","list","logical","numeric")
		argumentValid(check, ref)
		ct <- ctSelection(data, frame, group)
		if (htsfilter == TRUE) {
			htsfilter <- HTSFilter(ct, group[["factors"]], s.min=1, s.max=200, s.len=25)
			htsfiltered <- htsfilter$filteredData

			if (cfilter > 0) {
				ctcfilter <- cFilter(htsfiltered, cfilter, group[["factors"]])
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
	else {
		check <- list(data, frame, group, htsfilter, cfilter)
		ref <- list("data.frame","data.frame","list","logical","numeric")
		argumentValid(check, ref)
		ct <- ctSelection(data, frame, group)
		ct <- ctFilter(data, frame, group, htsfilter, cfilter)
		y <- DGEList(counts=ct, group=group[["factors"]])
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
	else {
		check <- list(data, frame, group, htsfilter, cfilter)
		ref <- list("data.frame","data.frame","list","logical","numeric")
		argumentValid(check, ref)	
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
