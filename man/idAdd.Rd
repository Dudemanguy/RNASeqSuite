\name{idAdd}
\alias{idAdd}

\title{Adds additional gene annotiations to a DGEList.}
\description{
Adds a column of gene output with its description to the DGEList using geneConvert.
}
\usage{
idAdd(dl, species, input_id, output_id)
}

\arguments{
  \item{dl}{A DGEList with rownames that are valid gene annotations.}
  \item{species}{The species name. Must be compliant with NCBI names. Shortcut names "mouse", "rat", and "human" are hardcoded as valid alternatives for convenience.}
  \item{input_id}{The valid geneConvert input id for the query. A simplistic search is built into the function.}
  \item{output_id}{The valid geneConvert output id for the query. A simplistic search is built into the function.}
}

\examples{
y <- idAdd(y, "mouse", "geneid", c("symbol", "description"))
}
