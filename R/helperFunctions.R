#general functions to check argument input types and ensure they are valid

.argumentValid <- function(classnames, refnames) {
	classframe <- data.frame(names(classnames[]))
	classlist <- list()
	for (i in seq_along(classnames)) {
		classlist[[i]] <- class(classnames[[i]])
	}
	classvector <- unlist(classlist)
	classframe["class"] <- classvector
	classframe["ref"] <- refnames
	for (i in nrow(classframe)) {
		if (classframe[i, 2] != classframe[i, 3]) {
			stop(paste(classframe[i, 2], " is an invalid class for ", classframe[i, 1], ". It must be a ", classframe[i, 3], sep=""), call.=FALSE)
		}
	}
}

#general function to match strings across dataframes

.stringMatch <- function(data, index, strings) {
	for (i in 1:length(strings)) {
		if (!(strings[[i]] %in% data[,index])) {
			return(FALSE)
		}
	}
	return(TRUE)
}

#helper function for normalizing vectors

.norm_vector <- function(vec, norm_type) {
	normalized <- vec/norm(data.matrix(vec), type=norm_type)
	return(normalized)
}

#abstract selection operations for apply

.select <- function(sub, whole) {
	output <- whole[sub,]
	return(output)
}

#abstract operation to remove rows with rowsums equal to zeros

.removeZeros <- function(df) {
	df_nozero <- df[rowSums(df) != 0,]
}

#split by group and seperate dataframe into lists

.ctSplit <- function(data, frame, group) {
	selection_grep <- sapply(levels(group), '==', frame[,2])
	selection_list <- apply(selection_grep, 2, .select, frame)
	selection <- do.call("rbind", selection_list)
	columns <- as.vector(selection[,1])
	matframe <- subset(data, select=eval(parse(text=list(columns))))
	subframe <- data.frame()
	sublist <- list()
	j <- 1
	k <- 0
	for (i in 1:(nrow(selection))) {
		if (i != nrow(selection)) {
			if (selection[i, 2] == selection[i+1, 2]) {
				j <- j+1
				k <- k+1
			}
			if (!(selection[i, 2] == selection[i+1, 2])) {
				subframe <- matframe[,(j-k):j]
				sublist[[i]] <- subframe
				j <- j+1
				k <- 0
			}
		}
		if (i == nrow(selection)) {
			subframe <- matframe[,(j-k):j]
			sublist[[i]] <- subframe
		}
	}
	sublist[sapply(sublist, is.null)] <- NULL
	return(sublist)
}

#use biomaRt to convert gene ids

.idConvert <- function(gene, org, attr_in, attr_out) {
	ensembl <- useMart("ensembl")
	if (org == "mouse") {
		ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
	}
	else if (org == "human") {
		ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
	}
	else if (org == "rat") {
		ensembl <- useMart("ensembl", dataset="rnorvegicus_gene_ensembl")
	}
	else if (org %in% listDatasets(ensembl)$dataset) {
		ensembl <- useMart("ensembl", dataset=org)
	}
	else {
		print("Invalid organism input. Please choose 'mouse', 'human', 'rat' or the exact dataset name from biomaRt.")
		stop()
	}

	attr_matrix <- listAttributes(ensembl)

	if (attr_in %in% attr_matrix$name) {
		id_input <- attr_in
	}
	else {
		attr_in_grep <- grep(attr_in, attr_matrix$name, ignore.case=TRUE)
		attr_in_filtered <- attr_matrix[attr_in_grep,]
		if (nrow(attr_in_filtered) == 0) {
			print("No results found for attribute input. Please try again.")
			stop()
		}
		else {
			if (nrow(attr_in_filtered) > nrow(head(attr_in_filtered))) {
				print(head(attr_in_filtered))
				id_input <- readline("Type in the name of the input id from the list. Type 'show' to see full matrix. \n")
				while (id_input == "show") {
					print(attr_in_filtered)
					id_input <- readline("Type in the name of the input id from the list. \n")
				}
			}
			else {
				print(attr_in_filtered)
				id_input <- readline("Type in the name of the input id from the list. \n")
			}
		}
	}

	#handle the output attribute
	if (attr_out %in% attr_matrix$name) {
		id_output <- attr_out
	}
	else {
		attr_out_grep <- grep(attr_out, attr_matrix$name, ignore.case=TRUE)
		attr_out_filtered <- attr_matrix[attr_out_grep,]
		if (nrow(attr_out_filtered) == 0) {
			print("No results found for attribute output. Please try again.")
			stop()
		}
		else {
			if (nrow(attr_out_filtered) > nrow(head(attr_out_filtered))) {
				print(head(attr_in_filtered))
				id_output <- readline("Type in the name of the output id from the list. Type 'show' to see full matrix. \n")
				while (id_output == "show") {
					print(attr_out_filtered)
					id_output <- readline("Type in the name of the output id from the list. \n")
				}
			}
			else {
				print(attr_out_filtered)
				id_output <- readline("Type in the name of the output id from the list. \n")
			}
		}
	}

	converted <- getBM(attributes=c(id_input, id_output, "description"), filters=id_input, values=c(gene), mart=ensembl)
	return(converted)
}

#fetch additional annotation from org R packages and insert into DGEList
#TODO: remove hardcoded options for species and annotations

.annotationFetch <- function(count, species) {
	selection <- paste("org.", species, ".eg.db", sep="")
	library(selection, character.only=TRUE)
	if (selection == 'org.Mm.eg.db') {
		idfound <- count$genes$genes %in% mappedRkeys(org.Mm.egREFSEQ)
		count <- count[idfound,]
		}
		egREFSEQ <- toTable(org.Mm.egREFSEQ)
	if (selection == 'org.Hs.eg.db') {
		idfound <- count$genes$genes %in% mappedRkeys(org.Hs.egREFSEQ)
		count <- count[idfound,]
		egREFSEQ <- toTable(org.Hs.egREFSEQ)
		}
	if (species == 'org.Rn.eg.db') {
		idfound <- count$genes$genes %in% mappedRkeys(org.Rn.egREFSEQ)
		count <- count[idfound,]
		egREFSEQ <- toTable(org.Rn.egREFSEQ)
		}
	m <- match(count$genes$genes, egREFSEQ$accession)
	count$genes$EntrezGene <- egREFSEQ$gene_id[m]
	egSYMBOL <- toTable(org.Mm.egSYMBOL)
	m <- match(count$genes$EntrezGene, egSYMBOL$gene_id)
	count$genes$Symbol <- egSYMBOL$symbol[m]
	o <- order(count$table$PValue)
	count <- count[o,]
	d <- duplicated(count$genes$EntrezGene)
	count <- count[!d,]
	return(count)
}
