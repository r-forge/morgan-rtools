\name{ibdhap.names}
\alias{ibdhap.names}
\title{ Get chromosome names }
\description{
  Get identifying infomation of chromosomes making up each
set, this consists of a number identifying the individual and 0/1 indicator of which of their
two chromosomes was used.
}
\usage{
ibdhap.names(ids.file=NULL, ids.filename = NULL)
} 
\arguments{
  \item{ids.file}{ A matrix of the .ids file used to run ibd_haplo.If
    this is left as the default NULL, ids.filename should be specified
  } 
\item{ids.filename}{ The filename (location) of the ids file used to run
  ibd_haplo. This is to be input as a character string. The ids file
  will be loaded from this location if ids.file is not specified.
}
}
\value{
 Matrix with a row for each set of chromosomes and two columns for each
 constituent chromosome. The first column is the individual number and
 the second column indicates which of their two chromosomes is used. 
}
\author{
F L Grimson
}
\examples{
## this example is taken from the package vignette.
##See vignette(IBDhaploRtools_tutorial)

data(ids_phased)
ibdhap.names( ids.file = ids_phased )
}