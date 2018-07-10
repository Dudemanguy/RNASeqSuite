\name{write.output}
\alias{write.output}

\title{A wrapper function to easily export results to directory}
\description{
Exports the results in the DGEExact or DGELRT object to a directory
}
\usage{
write.output(dl, directory, fdr=0.05)
}

\arguments{
  \item{dl}{A valid DGEExact or DGELRT object}
  \item{directory}{The name of the directory to output results}
  \item{fdr}{The fdr value to use for a cutoff on significant results}
}

\examples{
write.output(dl, 'output', fdr=0.05)
}
