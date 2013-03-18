\name{removezeros}
\alias{removezeros}
\title{
  remove zeros
}
\description{
  remove zeros from a vector
}
\usage{
removezeros(X)
}
\arguments{
  \item{X}{
    a vector
}
}

\value{
  the same vector but with zeros removed
}

\examples{

X = c(1,0,2,0,3)

removezeros(X)

## equivalently
 X[X!=0]

}
