\name{ibdhap.make.states}
\alias{ibdhap.make.states}

\title{ ibdhap.make.states }
\description{
 Stores and simplifies the qibd files created by IBD Haplo by "calling" a marker to be in an ibd state if it's marginal
probability meets some "cutoff" value, a zero or "no call" is assigned
to a marker in which no single state meets the value assigned to "cutoff". The R data.frame that this
function creates is expected by other functions in this package.
}

\usage{
ibdhap.make.states(qibd.filename = NULL, ids.filename = NULL, qibd.file
= NULL, ids.file = NULL, cutoff = 0.8)
}

\arguments{
  \item{qibd.file}{ A matrix of the contents of the qibd.out file the is
  produced from running ibd_haplo. If this is left as the default NULL,
  qibd.filename should be specified}
  \item{qibd.filename}{ The filename (location) of the qibd.out file the
  is produced from running ibd_haplo. This is to be input as a character
  string. The qibd file will be loaded from this location if qibd.file
  is not specified.
  }
  \item{ids.file}{ A matrix of the .ids file used to run ibd_haplo.If
    this is left as the default NULL, ids.filename should be specified
  } 
\item{ids.filename}{ The filename (location) of the ids file used to run
  ibd_haplo. This is to be input as a character string. The ids file
  will be loaded from this location if ids.file is not specified.
}
  \item{cutoff}{A scalar value between 0 and 1 to act as the "cutoff"
    value. This is the value which, if the maximum marginal probability
  of an ibd state is greater than, that marker will be called that
  state. Otherwise, the marker is called as zero, which in this context,
  means that there was  not evidence enough to determine the specific
  ibd state of that marker. Default value is 0.8.
}
}
\value{
  Returns a data.frame with
  ncol = # of sets of haplotypes/ pairs of
  genotypes in the qibd file.
  nrow = # of markers/SNPs
  Each column of this data.frame consists of integers (0 - 15 for
  haplotypes, 0-9 for genotypic data) corresponding to the ibd state at
  that marker (if the probability of that state for the marker is
  maximal and exceeds "cutoff" value, or a 0 value (for no call).  
}

\author{
M.D. Brown
}
\note{  
    This data.frame is required for all other functions in this package,
  so calling this function first is required. 
}

\examples{

## this example is taken from the package vignette.
##See vignette(IBDhaploRtools_tutorial)

data(qibd_phased)
data(ids_phased)

phased.gold <- ibdhap.make.states( qibd.file  = qibd_phased,
                   ids.file= ids_phased, cutoff = 0.8)

## alternatively, specify the file location, e.g.
## qibd.filename <- '~/Documents/qibd_unphased_2011.gold'
## ids.filename <- '~/Documents/ids_unphased_2011.gold'
## phased.gold <- ibdhap.make.states( qibd.filename = qibd.filename, 
##                   ids.filename = ids.filename, cutoff = 0.8)

}