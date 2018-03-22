\name{gopathway}
\alias{gopathway}

\title{Integrates goana analysis with DataList}
\description{
Performs a pathway analysis using results from exactTest, glmLRTTest, or glmQLFTest. \code{idAdd} must be run first.
}
\usage{
gopathway(dl, species, fdr=0.05)
}

\arguments{
  \item{dl}{A DataList containing results from either the exactTest, glmLRTTest, or glmQLFTest as well as output from \code{idAdd}.}
  \item{species}{The species name. Must be a two character identifier. i.e. "Mm" for mouse.}
  \item{fdr}{The cutoff for pathways to be considered significant. The default is 0.05}
}

\examples{
go <- gopathway(y, species="Mm")
}
