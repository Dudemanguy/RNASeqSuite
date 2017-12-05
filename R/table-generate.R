#library(id-conversion.R)
print("Create table of mgi symbols and gene description from edgeR and DESeq2 results.")
query <- readline("Press 1 to read from edgeR results. Press 2 to read from DESeq2 results. \n")

if (query==1) {
	x <- read.delim("edgeR/significant_results")
	y <- ID_convert(rownames(x), "mouse", "refseq_mrna", "mgi_symbol")
	z <- ID_convert(rownames(x), "mouse", "refseq_ncrna", "mgi_symbol")
	
	#remove empty entries and combine the two tables
	y_noempty <- y[y$mgi_symbol != "",]
	z_noempty <- z[z$mgi_symbol != "",]

	rownames(y_noempty) <- 1:nrow(y_noempty)
	rownames(z_noempty) <- (nrow(y_noempty)+1):(nrow(y_noempty)+nrow(z_noempty))
	colnames(y_noempty)[1] <- "refseq"
	colnames(z_noempty)[1] <- "refseq"
	full_conversion <- rbind(y_noempty, z_noempty)

	#clean up some names and add the diff gene data to the matrix
	rownames(full_conversion) <- full_conversion$refseq
	full_conversion <- full_conversion[-(1)]
	x_filter <- x[rownames(x) %in% rownames(full_conversion),]
	complete_frame <- cbind(x_filter, full_conversion)

	#reorder columns and write to table
	refcols <- c("mgi_symbol", "description")
	complete_frame <- complete_frame[,c(refcols, setdiff(names(complete_frame), refcols))]
	write.table(complete_frame, file="filtered_edgeR_results", sep="\t", quote=FALSE)
}

if (query==2) {
	x <- read.delim("DESeq2/significant_results")
	y <- ID_convert(rownames(x), "mouse", "refseq_mrna", "mgi_symbol")
	z <- ID_convert(rownames(x), "mouse", "refseq_ncrna", "mgi_symbol")
	
	#remove empty entries and combine the two tables
	y_noempty <- y[y$mgi_symbol != "",]
	z_noempty <- z[z$mgi_symbol != "",]

	rownames(y_noempty) <- 1:nrow(y_noempty)
	rownames(z_noempty) <- (nrow(y_noempty)+1):(nrow(y_noempty)+nrow(z_noempty))
	colnames(y_noempty)[1] <- "refseq"
	colnames(z_noempty)[1] <- "refseq"
	full_conversion <- rbind(y_noempty, z_noempty)

	#clean up some names and add the diff gene data to the matrix
	rownames(full_conversion) <- full_conversion$refseq
	full_conversion <- full_conversion[-(1)]
	x_filter <- x[rownames(x) %in% rownames(full_conversion),]
	complete_frame <- cbind(x_filter, full_conversion) 

	#reorder columns and write to table
	refcols <- c("mgi_symbol", "description")
	complete_frame <- complete_frame[,c(refcols, setdiff(names(complete_frame), refcols))]
	write.table(complete_frame, file="filtered_DESeq2_results", sep="\t", quote=FALSE)
	}

#else {
	#print("Error. Please enter 1 for edgeR or 2 for DESeq2.")
