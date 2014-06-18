\name{ibdhap.transitions}
\alias{ibdhap.transitions}

\title{
  create transition matrix
}
\description{ Creates a matrix of transition counts from when ibd
state switches along the chromosome.
}
\usage{
ibdhap.transitions(calls, data.type = c("h", "g", "r"))
}

\arguments{
  \item{calls}{The data.frame created from running ibdhap.make.calls on ibd_haplo
    output. 
}
  \item{data.type}{
    "h" : haplotypic data
    "g" : genotypic data (or hap data ran as genotypic)
    "r" : reduced data (output from running ibdhap.reduce.states and
    then ibdhap.make.states) 
  }
  }
\details{
  To create this matrix, all no calls are ignored. This is because, when
  transitioning out of being relatively certain of an ibd state, the
  marginal probabilities of ibd state by marker
  usually move into a segment of uncertainty (hence no calls)
  before it becomes relatively certain of an ibd state and therefore switches states. 
}
\value{

A matrix of size 15 x 15 (haplotypic) or 9 x 9 (genotypic) that shows
counts of ibd state transitions. Element [i,j] of the output is the
number of times state i changed to state j in the data.
}

\author{
  MD Brown
}

\examples{
## this example is taken from the package vignette.
##See vignette(IBDhaploRtools_tutorial)

data(qibd_phased)
data(ids_phased)
data(trueibd_phased)

phased.gold <- ibdhap.make.calls( qibd.file  = qibd_phased,
                   ids.file= ids_phased, cutoff = 0.8)
transitions.phased <- ibdhap.transitions(phased.gold, data.type="h")
}
