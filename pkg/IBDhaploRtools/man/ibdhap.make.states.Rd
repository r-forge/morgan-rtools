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
ibdhap.make.states(qibd.filename, dat.filename, cutoff = 0.8)
}

\arguments{
  \item{qibd.filename}{ The filename of the qibd.out file the is produced from running
    ibd_haplo. This is to be input as a character string. In the
    examples, this file is called "qibd_h.gold". 
}
\item{dat.filename}{The filename of the .dat file used to run ibd_haplo. This is the
  complicated parameter file that consists of three lines of values /
  indicators that tells ibd_haplo what to do. In the examples, this
  files is called "compu_4hap.dat".  
}
  \item{cutoff}{A scalar value between 0 and 1 to act as the "cutoff" value. This is
  the value which, if the maximum marginal probability of an ibd state
  is greater than, that marker will be called that state. Otherwise, the
  marker is called as zero, which in this context, means that there was
  not evidence enough to determine the specific ibd state of that
  marker. Default value is 0.8. }


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


}