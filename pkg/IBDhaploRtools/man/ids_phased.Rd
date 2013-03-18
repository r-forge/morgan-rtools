\name{ids_phased}
\alias{ids_phased}
\alias{qibd_phased}
\alias{ids_unphased}
\alias{qibd_unphased}
\alias{trueibd_phased}
\alias{posvec}
\docType{data}

\title{Example MORGAN IBDhaplo output}
\description{
 These files are example output from the MORGAN IBDhaplo program. There
 is an ids and qibd file produced in each IBDhaplo run. This output is
 from two runs, each on the same original data one of which was
 genotypic (unphased) and one of which was haplotypic (phased). The true
 ibd states of the data are given in trueibd_phased, in terms of
 haplotypic IBD states.

 There is also the vector posvec which gives the positions of each of the
 2000 SNPs used in the data files. 

 Please refer to the MORGAN IBDhaplo documentation for information on
 how to generate these files and what they contain.  }


 \usage{ids_phased} 

\format{ids files contain a 4x11 matrix, qibd files contain an 8000x18 matrix.}
\source{MORGAN IBDhaplo}
\references{
  http://www.stat.washington.edu/thompson/Genepi/pangaea.shtml
}
\keyword{datasets}