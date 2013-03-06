\name{get.counts}
\alias{get.counts}
\title{
  get counts
}
\description{
  subroutine that counts transitions
}
\usage{
get.counts(x, state.num)
}
\arguments{
  \item{x}{
    vector of data
}
  \item{state.num}{
    scalar indicating number of states 
}
}
\value{
  returns a state.num by state.num matrix of transition counts
}

\examples{

x = sample( 1:5, size = 20, replace = TRUE)
get.counts(x, 5)

}
