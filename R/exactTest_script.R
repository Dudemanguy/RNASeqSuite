library(HTSFilter)
library(edgeR)
suppressPackageStartupMessages(library(DESeq2))
options(max.width=9999999, width=10000)

#define functions

#user input
print("Script for performing an exact test between two groups")
count_raw <- readline("Enter the file containing the count matrix \n")
count_matrix <- read.table(count_raw)
sample_query <- readline("Press 1 to enter samples manually. Press 2 to read from a file (each entry on a separate line) \n")

#inputting sample groups; provide manual option and reading from a file
if (sample_query==1){
	columns_to_keep <- readline("Enter the desired samples as numbers separated by a comma prefaced by a lowercase h (e.g. h001,h002,h003,etc.) \n")
	sample_split <- strsplit(columns_to_keep, ",")
	count_subset <- subset(count_matrix, select=eval(parse(text=sample_split)))

	#ask for condition entry
	condition_query2 <- readline("Press 1 to enter the condition manually. Press 2 to read from a file (each entry on a separate line) \n")
	if (condition_query2==1){
		define_group <- readline("Enter the condition for each sample separated by commas (e.g. 1,1,2,2). Note that the order must be the same as the samples. \n")
		group_split <- strsplit(define_group, ",")
		group <- unlist(group_split)
		group <- factor(group, levels=unique(unlist(group_split)))
	}
	
	else if (condition_query2==2){
		define_group <- readline("Enter the file name containing the condition \n")
		group <- readLines(define_group)
		group <- factor(group, levels=unique(readLines(define_group)))
	}
}

#sample entry
if (sample_query==2){
	sample_query2 <- readline("Press 1 if the file only contains the desired samples. Press 2 if the file also defines the condition (tab delimited) \n")
	if (sample_query2==1){
		print("File only contains sample information")
		sample_file <- readline("Enter the file name \n") 
		columns_to_keep <- readLines(sample_file)
		count_subset <- subset(count_matrix, select = eval(parse(text=list(columns_to_keep))))
		condition_query <- readline("Press 1 to enter the condition manually. Press 2 to read from a file (each entry on a separate line) \n")

		#ask for condition entry
			if (condition_query==1){
				define_group <- readline("Enter the condition for each sample separated by commas (e.g. 1,1,2,2). Note that the order must be the same as the samples. \n")
				group_split <- strsplit(define_group, ",")
				group <- unlist(group_split)
				group <- factor(group, levels=unique(unlist(group_split)))
			}

			else if (condition_query==2){
				define_group <- readline("Enter the file name containing the condition \n")
				group <- read.table(define_group)
				group <- factor(group, levels=unique(read.table(define_group)))
			}
	}

	else if (sample_query==2){
		print("File contains sample and condition information")
		complete_file <- readline("Enter the file name \n") 
		complete_table <- read.table(complete_file) 
		columns_to_keep <- as.vector(complete_table[,1])
		count_subset <- subset(count_matrix, select = eval(parse(text=list(columns_to_keep))))
		group <- complete_table[,2]
		group <- factor(group, levels=unique(complete_table[,2]))
	}
}


#error handling

if (exists("count_subset")==FALSE){
	print("Error. No count matrix found. Retry entry.")
	stop()
}

#get the number of samples in the first group
group_table <- table(group)
n <- group_table[[1]][1]

#output directory for results

dir_output <- readline("Enter the name of the output directory \n")
dir.create(dir_output) 

#apply HTSFilter
filter <- HTSFilter(count_subset, group, s.min=1, s.max=200, s.len=25)
count_filter <- filter$filteredData

#apply edgeR
y <- DGEList(counts=count_filter, group=group)
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)

#apply DESeq2
group_dataframe <- data.frame(group)
dds <- DESeqDataSetFromMatrix(countData=count_filter, colData=group_dataframe, design=~group)
dds <- DESeq(dds)
res <- results(dds)

#create data frame for edgeR results
et_results <- topTags(et, n=Inf, sort.by="none")
write.table(et_results, file="et_output", sep="\t")
et_table <- read.table("et_output")
et_table <- et_table[-(2)]
file.remove("et_output")

#create data frame for DESeq2 results
resOrdered <- res[order(res$padj),]
write.table(resOrdered, file="DESeq_output", sep="\t")
res_order <- read.table("DESeq_output")
file.remove("DESeq_output")

#grab counts, average them, add to data frame, and sort

a <- count_filter[,1:n]
b <- count_filter[,(n+1):ncol(count_filter)]
c <- data.frame(rowMeans(a))
d <- data.frame(rowMeans(b))
et_table["Avg Ct A"] <- c
et_table["Avg Ct B"] <- d
et_table <- et_table[,c(4,5,1,2,3)]
et_table <- et_table[order(et_table$FDR, decreasing=FALSE),]

res_order["Avg Ct A"] <- c
res_order["Avg Ct B"] <- d
res_order <- res_order[,c(7,8,1,2,3,4,5,6)]

#print results and save to file

setwd(dir_output)

dir.create("edgeR")
setwd("edgeR")
sink("dgelist_output")
print(y)
sink()
write.table(et_table, file="full_results", sep="\t", quote=FALSE)
et_sig <- et_table[which(et_table$FDR<0.05),]
write.table(et_sig, file="significant_results", sep="\t", quote=FALSE)
setwd("..")

dir.create("DESeq2")
setwd("DESeq2")
sink("deseq2_summary")
print(summary(res))
sink()
write.table(res_order, file="full_results", sep="\t", quote=FALSE)
res_sig <- res_order[which(res_order$padj<0.05),]
write.table(res_sig, file="significant_results", sep="\t", quote=FALSE)
setwd("..")
setwd("..")
