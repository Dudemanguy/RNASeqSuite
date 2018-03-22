\name{grpSelection}
\alias{grpSelection}

\title{Create a factor containg the group names}
\description{
Creates a list containing a factor of group names as well as a data frame of group names from the selected character vector
}
\usage{
grpSelection(frame, groupselect, column=1, multi=FALSE)
}

\arguments{
  \item{frame}{A data frame containing the selected samples and their respective groups}
  \item{groupselect}{A character vector containing the desired group names}
  \item{column}{Integer that selects which column of the data frame to use as the group names}
  \item{multi}{Boolean that if set true, will combine all entries in the data frame rowwise with "." as a separator. Intended for multifactoral analysis}
}

\examples{
group <- grpSelection(frame, c("group_A","group_B"))
}
