\name{ibdhap.compare.loci}
\alias{ibdhap.compare.loci}

\title{ ibdhap compare.loci }
\description{Compares inferred ibd state with
  simulated "true" states. Calculates the proportion of markers called
  in the correct state, false positives (i.e. inferring ibd shared when
  none is shared), false negatives (i.e. inferring no ibd shared when
  ibd sharing is present) and the proportion of no calls.  }

\usage{
ibdhap.compare.loci(calls, true, data.type)
}

\arguments{ \item{calls}{ The data.frame created from running
  ibdhap.make.calls on ibd_haplo output.  }

\item{true}{ The data frame of same dimension as calls, but with true (probably simulated) ibd states. 
}

\item{data.type}{
    "h" : haplotypic data
    "g" : genotypic data (or hap data ran as genotypic)
    "r" : reduced data (output from running ibdhap.reduce.states and
    then ibdhap.make.states) 
}



}
  


\value{

Returns a list consisting of three matrices called "all", "ibd",
"nonibd" and "categories".

The matrix for "all" contain a row for each of the
following quantities, for all the loci:
\item{Number of sites}{}
\item{Called Correctly}{}
\item{Called IBD Incorrectly}{}
\item{Called noIBD Incorrectly}{}

The matrix for "ibd" contains a row for each of the following
quantities, based on just the loci that are truly in an IBD state:
\item{Number of sites}{}
\item{Called Correctly}{}
\item{Called as wrong IBD}{}
\item{Called as no IBD}{}
\item{No-call}{}

The matrix for "nonibd" contains a row for each of the follwing
quantities, based on just the loci that are truly in an non-ibd state:
\item{Number of sites}{}
\item{Called Correctly}{}
\item{Called as IBD}{}
\item{No-call}{}

Finally, the matrix for "categories" tabulates the percentage of the
loci that are truly in, and called to be in, each IBD state. 
}
\author{
  Fiona Grimson
}

\examples{

## this example is taken from the package vignette.
##See vignette(IBDhaploRtools_tutorial)

data(qibd_phased)
data(ids_phased)
data(trueibd_phased)

phased.gold <- ibdhap.make.calls( qibd.file  = qibd_phased,
                   ids.file= ids_phased, cutoff = 0.8)

ibdhap.compare.loci(phased.gold, trueibd_phased, "h")
}




