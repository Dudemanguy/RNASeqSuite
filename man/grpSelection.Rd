\name{grpSelection}
\alias{grpSelection}

\title{Create a factor containg the group names}
\description{
Creates a factor of group names from the selected character vector.
}
\usage{
grpSelection(frame, groupselect)
}

\arguments{
  \item{frame}{A data frame containing the selected samples and their respective groups}
  \item{groupselect}{A character vector containing the desired group names}
}

\examples{
group <- grpSelection(frame, c("group_A","group_B"))
}
