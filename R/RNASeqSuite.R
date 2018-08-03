#Custom filter function based around standard deviation of normalized vectors
cFilter <- function(dflist, sd) {
	check <- list(dflist=dflist, sd=sd)
	ref <- c("list", "numeric")
	.argumentValid(check, ref)
	if (sd < 0) {
		stop(paste("Error,", deparse(substitute(sd)), "must be positive."))
	} else {
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

#create subset of ct matrix based upon entered group
ctSelection <- function(data, group) {
	check <- list(data=data, group=group)
	ref <- c("data.frame", "list")
	.argumentValid(check, ref)
	columns <- rownames(group$frame)
	matframe <- subset(data, select=eval(parse(text=list(columns))))
	matframe
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

#wrapper function around the exactTest to find differentially expressed genes and return them
exactWrapper <- function(data, group, htsfilter=TRUE, cfilter=0, cutoff=0, adjust.method="BH", sort.by="FDR", decreasing=FALSE) {
	check <- list(data=data, group=group, htsfilter=htsfilter, cfilter=cfilter, cutoff=cutoff)
	ref <- c("data.frame", "list", "logical", "numeric", "numeric")
	.argumentValid(check, ref)
	y <- quickDGE(data, group=group)
	y <- estimateDisp(y)
	dge <- exactTest(y)
	adj.p.val <- p.adjust(dge$table$PValue, method=adjust.method)
	dge$table$FDR <- adj.p.val
	o <- switch(sort.by,
		"logFC" = order(dge$table$logFC, decreasing=decreasing),
		"logCPM" = order(dge$table$logCPM, decreasing=decreasing),
		"F" = order(dge$table$F, decreasing=decreasing),
		"PValue" = order(dge$table$PValue, decreasing=decreasing),
		"FDR" = order(dge$table$FDR, decreasing=decreasing),
		"none" = 1:nrow(dge)
	)
	dge <- dge[o,]
	dge
}

#integrate goana analysis with DGEList
gopathway <- function(dl, species, fdr=0.05, separate=TRUE) {
	check <- list(dl=dl, species=species, fdr=fdr, separate=separate)
	ref <- c("DGEList", "character", "numeric", "logical")
	.argumentValid(check, ref)
	go_list <- list()
	go <- goana.DGELRT(dl, species=species, separate=separate)
	if (isTRUE(separate)) {
		go_up <- topGO(go, sort="Up", n=Inf)
		go_up
		go_list[["Up"]] <- go_up[go_up$P.Up <= fdr,]
		go_down <- topGO(go, sort="Down", n=Inf)
		go_list[["Down"]] <- go_down[go_down$P.Down <= fdr,]
		return (go_list)
	}
	if (!isTRUE(separate)) {
		go <- topGO(go, n=Inf)
		go <- go[-(7)]
		colnames(go) <- c("Term","Ont","N","Up","Down","P")
		go <- go[go$P <= fdr,]
		go <- go[order(go$P),]
		return (go)
	}
}

#obtains the group from a supplied tab-delimited list
grpSelection <- function(frame, groupselect, column=1, multi=FALSE) {
	check <- list(frame=frame, groupselect=groupselect)
	ref <- c("data.frame", "character")
	.argumentValid(check, ref)
	if (groupselect == "all") {
		groupselect <- unique(as.character(frame[,column]))
	}
	if (!(.stringMatch(frame, column, groupselect))) {
		stop(paste("Error, some entries in", deparse(substitute(frame)), 
					"do not match the arguments."))
	}
	selection_grep <- sapply(groupselect, '==', frame[,column])
	selection_list <- apply(selection_grep, 2, .select, frame)
	names(selection_list) <- NULL
	selection <- do.call("rbind", selection_list)
	getgroup <- list()
	if (multi == TRUE) {
		output <- do.call(paste, selection)
		output <- factor(gsub(" ", ".", output))
		getgroup$factor <- output
	} else {
		getgroup$factor <- selection[,column]
	}
	getgroup$frame <- selection
	getgroup
}

#add annotations to DGEList
idAdd <- function(dl, species, input_id, output_id) {
	check <- list(dl=dl, species=species, input_id=input_id, output_id=output_id)
	ref <- c("DGEList", "character", "character", "character")
	.argumentValid(check, ref)
	genes <- rownames(dl)
	values <- convert(genes, species, input_id, output_id)
	m <- match(rownames(dl), values[[input_id]])
	dl$genes <- values[output_id][m,]

	if (identical(class(dl)[1], "DGEExact")) {
		dl$table <- data.frame(values[output_id][m,], dl$table)
	}
	if (identical(class(dl)[1], "DGELRT")) {
		dl$table <- data.frame(values[output_id][m,], dl$table)
	}
	rownames(dl) <- genes
	dl
}

#wrapper function around glmQLFTest to find differentially expressed genes and return them
#TODO: Make the select argument work properly
qlfWrapper <- function(data, group, select=NULL, htsfilter=TRUE, cfilter=0, cutoff=0, robust=TRUE, coef=ncol(design),
					   contrast=NULL, poisson.bound=TRUE, adjust.method="BH", sort.by="FDR", decreasing=FALSE) {
	check <- list(data=data, group=group, select=select, htsfilter=htsfilter, cfilter=cfilter, cutoff=cutoff)
	ref <- c("data.frame", "list", "character", "logical", "numeric", "numeric")
	.argumentValid(check, ref)
	y <- quickDGE(data, group=group)
	design <- model.matrix(~0+group$factor)
	colnames(design) <- gsub(".*factor", "", colnames(design))
	y <- estimateDisp(y, design)
	glm <- glmQLFit(y, design)
	if (!is.null(select)) {
		if (!identical(length(select), as.integer(2))) {
			stop("Selection vector must be length 2 for parsing")
		} else {
			contrast <- makeContrasts(paste0(select[1], "-", select[2]), levels=design)
		}
	}
	lrt <- glmQLFTest(glm, coef=coef, contrast=contrast, poisson.bound=poisson.bound)
	adj.p.val <- p.adjust(lrt$table$PValue, method=adjust.method)
	lrt$table$FDR <- adj.p.val
	o <- switch(sort.by,
		"logFC" = order(lrt$table$logFC, decreasing=decreasing),
		"logCPM" = order(lrt$table$logCPM, decreasing=decreasing),
		"F" = order(lrt$table$F, decreasing=decreasing),
		"PValue" = order(lrt$table$PValue, decreasing=decreasing),
		"FDR" = order(lrt$table$FDR, decreasing=decreasing),
		"none" = 1:nrow(lrt)
	)
	lrt <- lrt[o,]
	lrt
}

#wrapper function that integrates ctFilter function with DGEList
quickDGE <- function(data, group, htsfilter=TRUE, cfilter=0, cutoff=0) {
	check <- list(data=data, group=group, htsfiter=htsfilter, cfilter=cfilter, cutoff=cutoff)
	ref <- c("data.frame", "list", "logical", "numeric", "numeric")
	.argumentValid(check, ref)
	ct <- ctFilter(data, group, htsfilter, cfilter, cutoff)
	y <- DGEList(ct, group=group$factor)
	y <- calcNormFactors(y)
	y
}

#make directory and output results
write.output <- function(dl, directory, fdr=0.05) {
	check <- list(dl=dl, directory=directory, fdr=fdr)
	ref <- c("DGEList", "character", "numeric")
	.argumentValid(check, ref)
	dir.create(directory, recursive=TRUE)
	origin <- getwd()
	setwd(directory)
	sink("datalist_output")
	print(dl)
	sink()
	write.table(dl$table, file="full_results", sep="\t", quote=FALSE)
	dl_cutoff <- dl$table[which(dl$table$FDR<fdr),]
	write.table(dl_cutoff, file="significant_results", sep="\t", quote=FALSE)
	setwd(origin)
}
