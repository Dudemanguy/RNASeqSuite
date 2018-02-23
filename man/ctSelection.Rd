\name{ctSelection}
\alias{ctSelection}

\title{Create subset of count matrix based on supplied group}
\description{
Subsets the count matrix into selected columns based on group.
}

\usage{
ctSelection(data, frame, group)
}

\arguments{
  \item{data}{A data frame of the raw counts}
  \item{frame}{A data frame containing the selected samples and their respective groups}
  \item{group}{The factor created from grpSelection containing the different group names}
}

\examples{
ct <- ctSelection(data, frame, group)
}
