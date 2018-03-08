#obtains the group from a supplied tab-delimited list

grpSelection <- function(frame, groupselect) {
	check <- list(frame=frame, groupselect=groupselect)
	ref <- c("data.frame", "character")
	.argumentValid(check, ref)
	if (!(.stringMatch(frame, 2, groupselect))) {
		stop(paste("Error, some entries in", deparse(substitute(frame)), 
					"do not match the arguments."))
	}
	else {
		selection_grep <- sapply(groupselect, '==', frame[,2])
		selection_list <- apply(selection_grep, 2, .select, frame)
		selection <- do.call("rbind", selection_list)
		rownames(selection) <- selection[,1]
		selection <- selection[,2, drop=FALSE]
		colnames(selection) <- c('group')
		getgroup <- list()
		getgroup[["factor"]] <- selection[,1]
		getgroup[["factor"]] <- factor(getgroup[["factor"]], levels=unique(getgroup[["factor"]]))
		getgroup[["frame"]] <- selection
		getgroup
	}
}

#create subset of ct matrix based upon entered group

ctSelection <- function(data, group) {
	check <- list(data=data, group=group)
	ref <- c("data.frame", "list")
	.argumentValid(check, ref)
	if (!(.stringMatch(group$frame, 1, group$factor))) {
		stop(paste("Error, no entries in", deparse(substitute(group)), "match the arguments."))
	}
	else {
		columns <- rownames(group$frame)
		matframe <- subset(data, select=eval(parse(text=list(columns))))
		matframe
	}
}

#Custom filter function based around standard deviation of normalized vectors

cFilter <- function(dflist, sd) {
	check <- list(dflist=dflist, sd=sd)
	ref <- c("list", "numeric")
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
		df_compressed
	}
}

#reads the supplied ct matrix of reads and group data; filters data according to group

ctFilter <- function(data, group, htsfilter=TRUE, cfilter=0, cutoff=0) {
	check <- list(data=data, group=group, htsfilter=htsfilter, cfilter=cfilter, cutoff=cutoff)
	ref <- c("data.frame", "list", "logical", "numeric", "numeric")
	.argumentValid(check, ref)
	ct <- ctSelection(data, group)
	name_list <- list()
	if (cutoff > 0) {
		ct_split <- .ctSplit(data, group)
		ct_splitavgs <- lapply(ct_split, rowMeans)
		ct_filter <- list()
		for (i in seq_along(ct_splitavgs)) {
			ct_filter[[i]] <- names(ct_splitavgs[[i]])[ct_splitavgs[[i]] > cutoff]
		}
		cutoff_names <- Reduce(intersect, ct_filter)
	}
	if (htsfilter == TRUE) {
		htsfilter <- HTSFilter(ct, group$factor, s.min=1, s.max=200, s.len=25)
		hts_names <- rownames(htsfilter[[1]])
	}
	if (cfilter > 0) {
		dflist <- .ctSplit(data, group)
		cfilter_names <- cFilter(dflist, cfilter)
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
	allfilter
}

#uses wrapper function around the exactTest to find differentially expressed genes and store them in a DataList object

exactWrapper <- function(data, group, htsfilter=TRUE, cfilter=0, cutoff=0) {
	check <- list(data=data, group=group, htsfilter=htsfilter, cfilter=cfilter, cutoff=cutoff)
	ref <- c("data.frame", "list", "logical", "numeric", "numeric")
	.argumentValid(check, ref)
	y <- DataList(counts=data, group=group)
	y <- calcNormFactors(y)
	y <- estimateDisp(y)
	y <- exactTest(y)
	y
}

#preliminary wrapper for using the glmQLFTest

edgeRGLM <- function(data, group, htsfilter=TRUE, cfilter=0, cutoff=0) {
	check <- list(data=data, group=group, htsfilter=htsfilter, cfilter=cfilter, cutoff=cutoff)
	ref <- c("data.frame", "list", "logical", "numeric", "numeric")
	.argumentValid(check, ref)
	y <- DataList(counts=data, group=group)
	y <- calcNormFactors(y)
	design <- model.matrix(~group$factor)
	y <- estimateDisp(y, design, robust=TRUE)
	fit <- glmQLFit(y, design, robust=TRUE)
	qlf <- glmQLFTest(fit, coef=c(2, ncol(design)))
	qlf_raw <- topTags(qlf, n=Inf)
	qlf_frame <- qlf_raw[[1]]
	qlf_frame
}

#use DESeq2 to compute a Wald test and find differentially expressed genes

DESeq2 <- function(data, group, htsfilter=TRUE, cfilter=0, cutoff=0) {
	check <- list(data=data, group=group, htsfilter=htsfilter, cfilter=cfilter, cutoff=cutoff)
	ref <- c("data.frame", "list", "logical", "numeric", "numeric")
	.argumentValid(check, ref)	
	ct <- ctFilter(data, group, htsfilter, cfilter, cutoff)
	groupframe <- data.frame(group$factor)
	colnames(groupframe) <- c("groupframe")
	dds <- DESeqDataSetFromMatrix(countData=ct, colData=groupframe, design=~groupframe)
	dds <- DESeq(dds)
	res <- data.frame(results(dds))
	width <- table(group$factor)[[1]]
	a <- ct[,1:width]
	b <- ct[,(width+1):ncol(ct)]
	c <- data.frame(rowMeans(a))
	d <- data.frame(rowMeans(b))
	res["Avg Ct A"] <- c
	res["Avg Ct B"] <- d
	resOrder <- res[,c(7, 8, 1, 2, 3, 4, 5, 6)]
	resOrder <- resOrder[order(resOrder$padj, decreasing=FALSE),]
	resOrder
}

#add annotations to DataList
#TODO: Write a better way to rearrange table order

idAdd <- function(dl, species, input_id, output_id, qlf=TRUE) {
	check <- list(dl=dl, species=species, input_id=input_id, output_id=output_id)
	ref <- c("DataList", "character", "character", "character")
	.argumentValid(check, ref)
	biomart <- .idConvert(rownames(dl), species, input_id, output_id)
	if (is.null(dl$genes)) {
		dl$genes <- data.frame(Symbol=rownames(dl))
		dl$genes$Symbol <- NULL
	}
	m <- match(rownames(dl), biomart[,1])
	dl$genes$Symbol <- biomart$mgi_symbol[m]
	dl$genes$Description <- biomart$description[m]

	#check for the existence of each statistical test and add to the data frame
	if (!(is.null(dl$et_results))) {
		dl$et_results["Symbol"] <- dl$genes$Symbol
		dl$et_results["Description"] <- dl$genes$Description
		dl$et_results <- dl$et_results[,c(7, 8, 1, 2, 3, 4, 5, 6)]
	}
	if (!(is.null(dl$lrt_results))) {
		dl$lrt_results["Symbol"] <- dl$genes$Symbol
		dl$lrt_results["Description"] <- dl$genes$Description
		dl$lrt_results <- dl$lrt_results[,c(6, 7, 1, 2, 3, 4, 5)]
	}
	if (!(is.null(dl$qlf_results))) {
		dl$qlf_results["Symbol"] <- dl$genes$Symbol
		dl$qlf_results["Description"] <- dl$genes$Description
		dl$qlf_results <- dl$qlf_results[,c(6, 7, 1, 2, 3, 4, 5)]
	}
	dl
}

#integrate goana analysis with DataList

gopathway <- function(dl, species, fdr=0.05) {
	check <- list(dl=dl, species=species, fdr=fdr)
	ref <- c("DataList", "character", "numeric")
	.argumentValid(check, ref)
	rownames(dl) <- dl$genes$EntrezGene
	go_list <- list()
	go <- goana.DataList(dl, species=species)
	go_up <- topGO(go, sort="Up", n=Inf)
	go_list[["Up"]] <- go_up[go_up$P.Up <= fdr,]
	go_down <- topGO(go, sort="Down", n=Inf)
	go_list[["Down"]] <- go_down[go_down$P.Down <= fdr,]
	go_list
}

#make directory and output results

write.output <- function(dl, directory, fdr=0.05) {
	check <- list(dl=dl, directory=directory, fdr=fdr)
	ref <- c("DataList", "character", "numeric")
	.argumentValid(check, ref)
	dir.create(directory)
	setwd(directory)
	sink("datalist_output")
	print(dl)
	sink()
	write.table(dl$et_results, file="full_results", sep="\t", quote=FALSE)
	dl_cutoff <- dl$et_results[which(dl$et_results$FDR<fdr),]
	write.table(dl_cutoff, file="significant_results", sep="\t", quote=FALSE)
	setwd('..')
}
