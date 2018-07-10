\name{exactWrapper}
\alias{exactWrapper}

\title{Wrapper function to easily perform an exactTest analysis on genes}
\description{
Performs an exactTest and stores the results in a DGEExact object.
}
\usage{
exactWrapper(data, frame, group, htsfilter=TRUE, cfilter=0, cutoff=0)
}

\arguments{
  \item{data}{A data frame of the raw counts}
  \item{frame}{A data frame containing the selected samples and their respective groups}
  \item{group}{The factor created from grpSelection containing the different group names}
  \item{htsfilter}{A boolean determining whether or not to use HTSFilter to filter count matrix}
  \item{cfilter}{A numeric value controlling the use of cfilter, a custom outlier detector}
  \item{cutoff}{Filters based on a hard cutoff of read counts}
}

\details{
\code{htsfilter} calls the HTSFilter library if set to TRUE. This will remove genes based on a similarity ranking via a Jaccard index. Genes with constant, low counts will be removed.

\code{cfilter} is a custom outlier detector. cfilter normalizes the count vectors and then takes the difference between the two groups. The standard deviation and mean of the difference vector is stored and the genes are subsequently filtered based on the value of cfilter. Gene count values that lie outside of the standard deviation times the cfilter will be filtered. The lower the value of cfilter, the more stringent the filter will be. The default value of 0 disables this filter.

\code{cutoff} is a hard count cutoff. The average of every group is calculated and if the average for all groups falls below the specified cutoff value, the gene is removed. The default value of 0 disables this filter.
}

\examples{
y <- exactWrapper(data, frame, group, htsfilter=TRUE, cfilter=2, cutoff=100)
}
