## RNASeqSuite

RNASeqSuite is a set of handy cli wrapper functions in R around common RNASeq processing libraries to ease the workflow, simplify user input, and allow for quick, on-the-fly switching of algorithms and statistical tests. The main focus is on differential gene expression. RNASeqSuite depends on edgeR, DESeq2, and HTSFilter for statistical backends. Additionally, RNASeqSuite uses [geneConvert](https://github.com/Dudemanguy911/geneConvert), an integrated XML parsing and SQL database builder, for internal conversion of gene annotations for further analysis.

## Installation
Using the R package, ``devtools`` run ``devtools::install_github("Dudemanguy911/RNASeqSuite")``.

## Quick Usage
To use RNASeqSuite, one only needs a data frame of read counts and data frame of sample/group information. `featureCounts` in the `subread` package is recommended for generating a read count matrix. Creating a simple, tab-delimited list with one column containing the columns in the read count matrix and the other column containing the corresponding group is also recommended. The rownames for the group data frame need to be the names of the samples for ``grpSelection`` to work properly. To quickly perform a Fisher's Exact Test, one only needs to do the following.

```
data <- read.delim('/path/to/count/matrix')
frame <- read.delim('/path/to/colnames/and/group')
group <- grpSelection(frame, c('name of group A', 'name of group B'))
y <- exactWrapper(data, group)
```

``y`` is the object containing the adjusted p-values for the genes. By default, exactWrapper will filter data for you based on ``HTSFilter`` which uses a Jaccard index to remove low, constant expression genes. There are more options available in the manpages. A ``qlfWrapper`` function is also available for edgeR's ``glmQLFTest``. Additional gene annotation can be added to the results object by the use of the ``idAdd`` function.  

`` y <- idAdd(y, species, input, output)``

``idAdd`` uses geneConvert for gene conversions. Out of the box, it has support for human, mouse, and rat genes. However, other organisms can be added. See the geneConvert [page](https://github.com/Dudemanguy911/geneConvert) for more details. To write results to a file, use the write.output wrapper.

``write.output(y, Directory, fdr)``

This outputs the console output of the object, the full results table from the exact test, and the significant results table from the exact test to separate files within the specified directory. The fdr is merely a cutoff for the significant results file. It defaults to 0.05, but can be specified to something else.

## License
GPLv3
