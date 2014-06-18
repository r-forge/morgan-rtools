\name{ibdhap.make.true}
\alias{ibdhap.make.true}
\title{
  ibdhap make true states
}
\description{
  This function reads ibd_haplo output of true (simulated) ibd states
  into an R matrix data structure. The true states are compared to the
  state calls.
}
\usage{
ibdhap.make.true( true.filename )
}
\arguments{
  \item{true.filename}{ The filename of the true pairwise ibd states
    from, for example, outfifteenibd.txt or outnineibd.txt}
}
\value{
  An R matrix data structure with a column for each pairwise comparison
  in the input file. The pair names are thrown away.
}
\examples{

## For an existing file called "outfifteenibd.txt" use
## ibdhap.make.true( "outfifteenibd.txt" )

## An example of a data set already read into R is
## data( trueibd_phased )

}
