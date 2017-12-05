#simple R function to convert gene IDs from one format to another; based on biomaRt
library(biomaRt)

print("R function to easily convert gene IDs across different formats/databases.")
print("Usage: ID_convert(gene_list, 'organism', 'input_attribute', 'output_attribute', desc=TRUE)")
ID_convert <- function(x, org, attribute_in, attribute_out, desc){
	
	if(missing(desc)) {
		desc = TRUE
	}

	else if (desc != (TRUE | FALSE)) {
		print("Error, desc is a boolean value.")
		stop()
	}

	ensembl <- useMart("ensembl")

	#handle the organism argument
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

	attribute_matrix <- listAttributes(ensembl)

	#handle the input attribute
	if (attribute_in %in% attribute_matrix$name) {
		id_input <- attribute_in
	}

	else {
		attribute_in_grep <- grep(attribute_in, attribute_matrix$name, ignore.case=TRUE)
		attribute_in_filtered <- attribute_matrix[attribute_in_grep,]
		
		if (nrow(attribute_in_filtered) == 0) {
			print("No results found for attribute input. Please try again.")
			stop()
		}

		else {
			print(head(attribute_in_filtered))
			id_input <- readline("Type in the name of the input id from the list. Type 'show' to see full matrix. \n")

				while (id_input == "show") {
					print(attribute_in_filtered)
					id_input <- readline("Type in the name of the input id from the list. \n")
				}
		}
	}

	#hande the output attribute
	if (attribute_out %in% attribute_matrix$name) {
		id_output <- attribute_out
	}

	else {
		attribute_out_grep <- grep(attribute_out, attribute_matrix$name, ignore.case=TRUE)
		attribute_out_filtered <- attribute_matrix[attribute_out_grep,]

		if (nrow(attribute_out_filtered) == 0) {
			print("No results found for attribute output. Please try again.")
			stop()
		}

		else {
			print(head(attribute_out_filtered))
			id_output <- readline("Input the desired output of the gene list. Type 'show' to see full matrix. \n")

				while (id_output == "show") {
					print(attribute_out_filtered)
					id_output <- readline("Input the desired output of the gene list. \n")
				}
		}				
	}

	#return output
	if (desc == TRUE) {
		converted <- getBM(attributes=c(id_input, id_output, "description"), filters=id_input, values=c(x), mart=ensembl)
		return(converted)
	}

	else if (desc == FALSE) {
		converted <- getBM(attributes=c(id_input, id_output), filters=id_input, values=c(x), mart=ensembl)
		return(converted)
	}
}
