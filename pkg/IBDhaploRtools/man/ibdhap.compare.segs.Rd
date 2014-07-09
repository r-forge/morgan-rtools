\name{ibdhap.compare.segs}
\alias{ibdhap.compare.segs}
\title{ ibdhap compare segments }
\description{ Compares inferred ibd states with simulated "true"
  states. Calculates segments of ibd in the true data, giving
  descriptions of the segments and the proportion of correct and
  incorrect calls within the segments.}
\usage{
ibdhap.compare.segs(calls, true, data.type, seg.cutoff, pos=NA)
}
\arguments{
\item{calls}{The data.frame created from running
  ibdhap.make.calls on ibd_haplo output.  }
\item{true}{The data frame of same dimension as calls, but with true (probably simulated) ibd states. 
}
\item{data.type}{
    "h" : haplotypic data
    "g" : genotypic data (or hap data ran as genotypic)
    "r" : reduced data (output from running ibdhap.reduce.states and
    then ibdhap.make.states) 
}
\item{seg.cutoff}{A scalar value between 0 and 1 to act as the cutoff
  value, that is, the percentage of loci in the segment whose calls must agree
  to call the segment for a particular ibd state.}
\item{pos}{A position vector with the same length as the number of loci
  (rows of calls or true) describing the positions in cM, M or any other
metric of each marker. If positions are not included, the segment
lengths reported will be the number of markers making up a segment}
}
\value{

  Returns a list containing two elememts "seg.stats", "seg.info".

  The "seg.stats" matrix has a row for each of the following statistics
  \item{Number of segments}{}
  \item{Number of IBD segments}{}
  \item{IBD segs called correctly}{}
  \item{IBD segs called no-IBD}{}
  \item{IBD segs called wrong IBD}{}
  \item{IBD segs with no call}{}
  
  The "seg.info" matrix has a row for each segment of IBD in the true
  data and a column for each of the following statistics for each
  segment.
  \item{seg.length}{}
  \item{true.state}{}
  \item{seg.call}{}
  \item{mode.call}{}
  \item{prop.corr}{}
}
\author{
  Fiona Grimson
}
\examples{

## this example is taken from the package vignette,
## See vignette(IBDhaploRtools_tutorial)

data(qibd_phased)
data(ids_phased)
data(trueibd_phased)

phased.gold <- ibdhap.make.calls( qibd.file  = qibd_phased,
                   ids.file= ids_phased, cutoff = 0.8)

ibdhap.compare.segs(phased.gold, trueibd_phased, "h", 0.8, pos=NA)
}
