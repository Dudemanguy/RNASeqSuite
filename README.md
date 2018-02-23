## RNASeqSuite

**Note**: This is a WIP and essentially an alpha.

RNASeqSuite is intended to be a handy cli wrapper in R around common RNASeq processing libraries to ease the workflow, simplify user input, and allow for quick, on-the-fly switching of algorithms and statistical tests. The main focus is on differential gene expression. Currently depends HTSFilter for backend calcuations. Statistical backends are imported from edgeR with some slight modifications to accommodate for the ``DataList`` class which is based heavily on edgeR's ``DGEList``. Currently, RNASeqSuite only supports pairwise comparisons.

## Installation
Using the R package, ``devtools`` run ``devtools::install_github("Dudemanguy911/RNASeqSuite")``.

## Quick Usage
To use RNASeqSuite, one only needs a data frame of read counts and data frame of sample/group information. `featureCounts` in the `subread` package is recommended for generating a read count matrix. Creating a simple, tab-delimited list with one column containing the columns in the read count matrix and the other column containing the corresponding group is also recommended. To quickly perform a Fisher's Exact Test, one only needs to do the following.

```
data <- read.table('/path/to/count/matrix')
frame <- read.table('/path/to/colnames/and/group')
group <- grpSelection(frame, c('name of group A', 'name of group B'))
y <- edgeRclassic(data, frame, group)
```

``y`` is the DataList object containing all the information on the count matrix, samples, pvalues, comparison, etc. By default, edgeRclassic will filter data for you based on ``HTSFilter`` which uses a Jaccard index to remove low, constant expression genes. There are more options available in the manpages. To quickly write results to a file, use the write.output wrapper.

``write.output(DataList, Directory, fdr)``

This outputs the console output of the Datalist, the full results table from the exact test, and the significant results table from the exact test to separate files within the specified directory. The fdr is merely a cutoff for the significant results file. It defaults to 0.05, but can be specified to something else.

## License
GPLv2 or later.
