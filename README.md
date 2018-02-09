## RNASeqSuite

**Note**: This is a WIP and is only here for convenience.

RNASeqSuite is intended to be a handy cli wrapper in R around common RNASeq processing libraries to ease the workflow, simplify user input, and allow for quick, on-the-fly switching of algorithms and statistical tests. The main focus is on differential gene expression. Currently depends on edgeR, DESeq2, and HTSFilter for backend calcuations. This will likely be phased out in the future.

## Installation
Using the R package, ``devtools`` run ``devtools::install_github("Dudemanguy911/RNASeqSuite")``.

## Preliminary Inputs
RNASeqSuite needs only a few inputs in order to make use of its internal functions. At this time, RNASeqSuite does not create a data frame of read counts, so that must be supplied in R (`subread` is recommended for obtaining the data frame). In addition, RNASeqSuite needs a data frame with one column containing the column names of the matrix of read counts, and the other column containing the corresponding group/condition. Currently, RNASeqSuite only supports pairwise differential expression.

## Available Functions

``grpSelection(frame, groupselect)``

Creates a list selecting the desired groups to use for comparisons. "frame" is a data frame containing all group information while "groupselect" is a character vector the names of the groups selected for comparison.

``ctSelection(data, frame, group)``

Simple function to easily subset read counts based on selected groups. "data" denotes the matrix of read counts, "frame" is the data frame containing group and column information, and "group" is the is the factor object of groups obtained from the grpSelection function.

``ctFilter(data, frame, group, htsfilter, cfilter, cutoff)``

Creates a data frame of filtered counts from the data frame of read counts. "data" is the read counts, "frame" is the data frame containing groups, and "group" is the factor object of groups obtained by the grpSelection function. "htsfilter" is a boolean value that enables the use of HTSFilter, an R package that filters counts based on the Jaccard index. This argument is optional and if omitted, defaults to TRUE. "cfilter" is a positive, numeric value that is used to remove outliers. Outlier detection is based on an internal function that calculates normalized vectors and removes those that are outside the standard deviation. The cfilter value is simply a multiplier of the standard deviation, so a cfilter=1 input will preserve all vectors whose values lie within a standard deviation of the mean. A lower value of cfilter, will be more stringent. A value of 0 will disable this filter. This argument is optional and if omitted, defaults to 0 (disabling this filter). "cutoff" is an optional hard numeric cutoff on the amount of counts. The average count value for each group is calculated and if the average value falls below the specified threshold for every group, the gene is removed. A value of 0 disables this filter.

``edgeRclassic(data, frame, group, htsfilter, cfilter, cutoff)``

Inputs are exactly the same as `ctFilter`. The `edgeR` function calls a DGEList from the edgeR library and returns a data frame containing differential genes ordered by their FDR along with count numbers. Currently, this function is limited to classic mode and only pairwise comparisons. It will only handle the default values of edgeR's internal functions and will likely be phased out in the future.


``DESeq2(data, frame, group, htsfilter, cfilter, cutoff)``

Inputs are exactly the same as `ctFilter`. The `DESeq2` function calls the DESeq function from the DESeq2 library and returns a data frame containing differential genes ordered by their padj along with count numbers. Currently, this function is limited to only pairwise comparisons. It will only handle the default values of DESeq2's internal functions and will likely be phased out in the future.

``write.table(dge, directory, fdr)``

A convenience wrapper function to output results from the edgeRclassic analysis. `dge` is the DGEList generated frome edgeRclassic. `directory` is the name of the desired directory for output. `fdr` denotes the cutoff fdr for the significant results output. `write.table` outputs 3 text files: the output from the DGEList, the full table of results from the exactTest, and the table of results with the fdr cutoff.

## License
GPLv2 or later.
