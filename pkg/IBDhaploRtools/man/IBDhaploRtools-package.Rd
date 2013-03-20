\name{IBDhaploRtools-package}
\alias{IBDhaploRtools-package}
\alias{IBDhaploRtools}
\docType{package}
\title{
Functions for the Analysis of IBD Haplo output
}
\description{
IBDhaploRtools consists of several functions to analyze, plot, and store the output of the IBD_Haplo software package. More information regarding IBD_Haplo can be found at www.stat.washington.edu/thompson/Genepi/pangaea.shtml.
}
\details{
\tabular{ll}{
Package: \tab IBDhaploRtools\cr
Type: \tab Package\cr
Version: \tab 1.5\cr
Date: \tab 2013-03-18\cr
License: \tab GPL-3\cr
LazyLoad: \tab yes\cr
}

Main Functions:

ibdhap.make.states : stores and simplifies the qibd files created by
IBD Haplo by "calling" a marker to be in an ibd state if it's marginal
probability meets some "cutoff" value. The R data.frame that this
function creates is expected by other functions in this package.

ibdhap.reduce.states : reduces the number of ibd states by grouping the  15 or 9(depending on running ibd_haplo on haplotypes or genotypes resp.) to an
easier to manage set of 4 ibd groupings. Output is a .txt file that
resembles that of output of ibd_haplo.

ibdhap.summary : summarizes the data created by ibdhap.make.states by
giving mean lengths of ibd segments, mean proportions of ibd shared, and
counts on ibd segments. This gives information on the group of sets of
haplotypes/genotypes or on an individual pairing.

ibdhap.seg.lengths : given the ibd states from a set of haplotypes/pair of genotypes (taken
from a column of the output of ibdhap.make.states), this function creates a
data.frame consisting of all segments of differing ibd state, paired
with their respective length.

ibdhap.transitions: creates a matrix of transition counts from when ibd
state switches along the chromosome.

ibdhap.barplot : Graphically displays regions of any ibd sharing, no ibd
sharing and no calls along a chromosome for a set of haplotypes / pair of genotypes.

ibdhap.correct.bylen : if simulated data(true ibd state) is
available, this function returns a data frame of how well ibd_haplo infers
segments based on segment length. 
}

\author{
  Marshall D. Brown

  Fiona Grimson
  <fgrimson@uw.edu>

}
\references{
  none
}


\examples{
# See the tutorial that  outlines the use of all the
# functions found in this package: 
# vignette( IBDhaploRtools_tutorial )
}
