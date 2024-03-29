% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2_rankings.R
\name{ffs_latest_rankings}
\alias{ffs_latest_rankings}
\title{Download latest rankings from DynastyProcess GitHub}
\usage{
ffs_latest_rankings(type = c("draft", "week"))
}
\arguments{
\item{type}{one of "draft" or "week" - controls whether to pull preseason or inseason rankings.}
}
\value{
a dataframe with a copy of the latest FP rankings from DynastyProcess's data repository
}
\description{
Fetches a copy of the latest FantasyPros redraft positional rankings data from DynastyProcess.com's data repository.
}
\details{
If you have any issues with the output of this data, please open an issue in
the DynastyProcess data repository.
}
\examples{
\donttest{
try({ # try block to prevent CRAN-related issues
ffs_latest_rankings()
})
}

}
\seealso{
\url{https://github.com/dynastyprocess/data}

\code{vignette("custom")} for example usage
}
