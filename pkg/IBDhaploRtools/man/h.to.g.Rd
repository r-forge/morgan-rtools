\name{h.to.g}
\alias{h.to.g}
\title{
 Convert haplotype states to genotype states
}
\description{
Convert haplotype states (1 - 15) to genotype states (1 - 9)
}
\usage{
h.to.g(hap.states)
}
\arguments{
  \item{hap.states}{
vector or matrix of haplotype states. Values should be integers from
1-9
}
}

\value{
 vector or matrix of same dimension as hap.states but with the
 corresponding genotypic states. 
}
\author{
F L Grimson
}
\examples{
## this example is taken from the package vignette.
##See vignette(IBDhaploRtools_tutorial)

data(trueibd_phased)
trueibd_unphased <- h.to.g( trueibd_phased )
}