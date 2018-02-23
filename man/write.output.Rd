\name{write.output}
\alias{write.output}

\title{A wrapper function to easily export results to directory}
\description{
Exports the results in the DataList object to a directory
}
\usage{
\method{write.output}(dl, directory, fdr=0.05)
}

\arguments{
  \item{dl}{A data list object containing results}
  \item{directory}{The name of the directory to output results}
  \item{fdr}{The fdr value to use for a cutoff on significant results}
}

\examples{
write.output(dl, 'output', fdr=0.05)
}
