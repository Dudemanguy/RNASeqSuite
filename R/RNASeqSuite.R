#create subset of ct matrix based upon entered group

ctSelection <- function(data, frame, group) {
	check <- list(data=data, frame=frame, group=group)
	ref <- c("data.frame","data.frame","factor")
	.argumentValid(check, ref)
	if (!(.stringMatch(frame, 2, group))) {
		stop(paste("Error, no entries in", deparse(substitute(groupframe)), "match the arguments."))
	}
	else {
		selection_grep <- sapply(levels(group), '==', frame[,2])
		selection_list <- apply(selection_grep, 2, .select, frame)
		selection <- do.call("rbind", selection_list)
		columns <- as.vector(selection[,1])
		matframe <- subset(data, select = eval(parse(text=list(columns))))
		return(matframe)
	}
}

#obtains the group from a supplied tab-delimited list

grpSelection <- function(frame, groupselect) {
	check <- list(frame=frame, groupselect=groupselect)
	ref <- c("data.frame","character")
	.argumentValid(check, ref)
	if (!(.stringMatch(frame, 2, groupselect))) {
		stop(paste("Error, some entries in", deparse(substitute(frame)), 
					"do not match the arguments."))
	}
	else {
		selection_grep <- sapply(groupselect, '==', frame[,2])
		selection_list <- apply(selection_grep, 2, .select, frame)
		selection <- do.call("rbind", selection_list)
		getgroup <- selection[,2]
		getgroup <- factor(getgroup, levels=unique(getgroup))
		return(getgroup)
	}
}

#Custom filter function based around standard deviation of normalized vectors

cFilter <- function(dflist, sd, group) {
	check <- list(dflist=dflist, sd=sd, group=group)
	ref <- c("list", "numeric", "factor")
	.argumentValid(check, ref)
	if (sd < 0) {
		stop(paste("Error,", deparse(substitute(sd)), "must be positive."))
	}
	else {
		combine <- do.call("cbind", dflist)
		combine_nozero <- .removeZeros(combine)
		gene_vector <- rownames(combine_nozero)
		list_filter <- list()
		list_norm <- list()
		for (i in seq_along(dflist)) {
			df_filter <- .select(gene_vector, dflist[[i]])
			list_filter[[i]] <- df_filter
		}
		for (i in seq_along(list_filter)) {
			norms <- apply(list_filter[[i]], 2, .norm_vector, 'f')
			list_norm[[i]] <- norms
		}
		list_avgs <- lapply(list_norm, rowMeans)
		avg_mat <- do.call("cbind", list_avgs)
		avg_row <- data.frame(rowMeans(avg_mat))
		diff <- apply(avg_mat, 1, max) - apply(avg_mat, 1, min)
		df_diff <- data.frame(diff)
		mean_df <- mean(df_diff[,1])
		sd_df <- sd(df_diff[,1])
		df_compressed <- rownames(df_diff[(mean_df - (sd*sd_df)) <= avg_row 
									& avg_row <= (mean_df + (sd*sd_df)), 1, drop=FALSE])
		return(df_compressed)
	}
}

#reads the supplied ct matrix of reads and group data; filters data according to group

ctFilter <- function(data, frame, group, htsfilter=TRUE, cfilter=0, cutoff=0) {
	check <- list(data=data, frame=frame, group=group, htsfilter=htsfilter, cfilter=cfilter, cutoff=cutoff)
	ref <- c("data.frame","data.frame","factor","logical","numeric","numeric")
	.argumentValid(check, ref)
	ct <- ctSelection(data, frame, group)
	name_list <- list()
	if (cutoff > 0) {
		ct_split <- .ctSplit(data, frame, group)
		ct_splitavgs <- lapply(ct_split, rowMeans)
		ct_filter <- list()
		for (i in seq_along(ct_splitavgs)) {
			ct_filter[[i]] <- names(ct_splitavgs[[i]])[ct_splitavgs[[i]] > cutoff]
		}
		cutoff_names <- Reduce(intersect, ct_filter)
	}
	if (htsfilter == TRUE) {
		htsfilter <- HTSFilter(ct, group, s.min=1, s.max=200, s.len=25)
		hts_names <- rownames(htsfilter[[1]])
	}
	if (cfilter > 0) {
		dflist <- .ctSplit(data, frame, group)
		cfilter_names <- cFilter(dflist, cfilter, group)
	}
	if (exists('cutoff_names')) {
		name_list[["cutoff"]] <- cutoff_names
	}
	if (exists('hts_names')) {
		name_list[["hts_names"]] <- hts_names
	}
	if (exists('cfilter_names')) {
		name_list[["cfilter_names"]] <- cfilter_names
	}
	total_names <- Reduce(intersect, name_list)
	allfilter <- .select(total_names, ct)
	return(allfilter)
}

#uses edgeR to compute an exact test and find differentially expressed genes
#TODO: Remove hardcoded .idConvert options

edgeRclassic <- function(data, frame, group, htsfilter=TRUE, cfilter=0, cutoff=0) {
	check <- list(data=data, frame=frame, group=group, htsfilter=htsfilter, cfilter=cfilter, cutoff=cutoff)
	ref <- c("data.frame","data.frame","factor","logical","numeric","numeric")
	.argumentValid(check, ref)
	ct <- ctFilter(data, frame, group, htsfilter, cfilter, cutoff)
	y <- DGEList(counts=ct, group=group, genes=rownames(ct))
	y <- calcNormFactors(y)
	y <- estimateDisp(y)
	et <- exactTest(y) 
	et_FDR <- topTags(et, n=Inf, sort.by="none")
	et$table["FDR"] <- et_FDR$table$FDR
	width <- table(group)[[1]]
	a <- ct[,1:width]
	b <- ct[,(width+1):ncol(ct)]
	c <- data.frame(rowMeans(a))
	d <- data.frame(rowMeans(b))
	biomart <- .idConvert(et, rownames(et$genes), 'mouse', 'refseq_mrna', 'mgi_symbol')
	et$table["Symbol"] <- biomart$genes$Symbol
	et$table["Description"] <- biomart$genes$Description
	et$table["Avg Ct A"] <- c
	et$table["Avg Ct B"] <- d
	et$table <- et$table[,c(5,6,7,8,1,2,3,4)]
	et$table <- et$table[order(et$table$FDR, decreasing=FALSE),]
	y$results <- et
	return(y)
}

#preliminary wrapper for using the glmQLFTest

edgeRGLM <- function(data, frame, group, htsfilter=TRUE, cfilter=0, cutoff=0) {
	check <- list(data=data, frame=frame, group=group, htsfilter=htsfilter, cfilter=cfilter, cutoff=cutoff)
	ref <- c("data.frame","data.frame","factor","logical","numeric","numeric")
	.argumentValid(check, ref)
	ct <- ctFilter(data, frame, group, htsfilter, cfilter, cutoff)
	y <- DGEList(counts=ct, group=group, genes=rownames(ct))
	y <- calcNormFactors(y)
	design <- model.matrix(~group)
	y <- estimateDisp(y, design, robust=TRUE)
	fit <- glmQLFit(y, design, robust=TRUE)
	qlf <- glmQLFTest(fit, coef=c(2,ncol(design)))
	qlf_raw <- topTags(qlf, n=Inf)
	qlf_frame <- qlf_raw[[1]]
	return(qlf_frame)
}

#use DESeq2 to compute a Wald test and find differentially expressed genes

DESeq2 <- function(data, frame, group, htsfilter=TRUE, cfilter=0, cutoff=0) {
	check <- list(data=data, frame=frame, group=group, htsfilter=htsfilter, cfilter=cfilter, cutoff=cutoff)
	ref <- c("data.frame","data.frame","factor","logical","numeric","numeric")
	.argumentValid(check, ref)	
	ct <- ctFilter(data, frame, group, htsfilter, cfilter, cutoff)
	groupframe <- data.frame(group)
	colnames(groupframe) <- c("groupframe")
	dds <- DESeqDataSetFromMatrix(countData=ct, colData=groupframe, design=~groupframe)
	dds <- DESeq(dds)
	res <- data.frame(results(dds))
	width <- table(group)[[1]]
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

#make directory and output results

write.output <- function(dge, directory, p) {
	check <- list(dge=dge, directory=directory, p=p)
	ref <- c("DGEList","character","numeric")
	.argumentValid(check, ref)
	dir.create(directory)
	setwd(directory)
	sink("dgelist_output")
	print(dge)
	sink()
	write.table(dge$results$table, file="full_results", sep="\t", quote=FALSE)
	dge_cutoff <- dge$results$table[which(dge$results$table$FDR<p),]
	write.table(dge_cutoff, file="significant_results", sep="\t", quote=FALSE)
	setwd('..')
}
