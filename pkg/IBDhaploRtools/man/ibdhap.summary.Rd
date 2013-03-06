\name{ibdhap.summary}
\alias{ibdhap.summary}

\title{
  ibdhap summary
}
\description{
  Summarizes the data created by ibdhap.make.states by
giving mean lengths of ibd segments, mean proportions of ibd shared, and
counts on ibd segments. This gives information on the group of sets of
haplotypes/genotypes or on an individual pairing.
}

\usage{
ibdhap.summary(states.dat, data.type = c("h", "g", "r"), position = NA)
}

\arguments{
  \item{states.dat}{ The data.frame created from running ibdhap.make.states on ibd_haplo
    output. 
}
  \item{data.type}{
    "h" : haplotypic data
    "g" : genotypic data (or hap data ran as genotypic)
    "r" : reduced data (output from running ibdhap.reduce.states and
    then ibdhap.make.states) 
}

  \item{position}{ A position vector, with the same length as
    nrow(states.dat) describing the positions (in cM, M, or any other
    metric) of each marker. If positions is not included, "segment
    length" refers to number of SNPs making up a segment.  } } \details{
    ibdhap.summary analyzes all the data it is given.If states.dat
    consists of more than one column, means are calculated across all
    sets of haplotypes. For summaries on only one set of haplotypes/pair
    of genotypes, just pass this function the column of data
    corresponding to the specific set on which you want summary
    statistics.  }

\value{

Returns a list consisting of:
\item{mean.prop }{ mean proportion of chromosomes called in any ibd
  shared, no ibd shared, and no calls }
\item{mean.length}{mean lengths of segments of chromosome called in any ibd
  shared, no ibd shared, and no calls. Length is measured by the
  position vector (see above).}
\item{seg.counts}{ total counts of segments of chromosome called in any ibd
  shared, no ibd shared, and no calls}

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

phased.gold <- ibdhap.make.states( qibd.file  = qibd_phased,
                   ids.file= ids_phased, cutoff = 0.8)

summary.phased <- ibdhap.summary( phased.gold, data.type="h")
}
