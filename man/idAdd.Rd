\name{idAdd}
\alias{idAdd}

\title{Adds additional columns to qlf or et results containing gene descriptions.}
\description{
Adds a column of gene output with its description to the DataList using biomaRt.
}
\usage{
idAdd(dl, species, input_id, output_id)
}

\arguments{
  \item{dl}{A DataList containing results from either the exactTest, glmLRTTest, or glmQLFTest.}
  \item{species}{The species name. Must be compliant with biomaRt ensembl names. Shortcut names "mouse", "rat", and "human" are valid.}
  \item{input_id}{The valid biomaRt input id for the query. A simplistic search is built into the function.}
  \item{output_id}{The valid biomaRt output id for the query. A simplistic search is built into the function.}
}

\examples{
y <- idAdd(y, "mouse", "refseq_mrna", "mgi_symbol")
}
