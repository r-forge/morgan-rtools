\name{ibdhap.barplot}
\alias{ibdhap.barplot}

\title{ create ibd barplot
}
\description{
   Graphically displays regions of any ibd sharing, no ibd
sharing and no calls along a chromosome for a set of haplotypes / pair of genotypes.
}

\usage{
ibdhap.barplot(x, data.type = c("h", "g", "r"), col = c("grey", "red",
"white"), position = NA, top = 1, bottom = 0, density = 50, ...)
}

\arguments{
  \item{x}{ A vector of ibd states (with values 0 - 15 for haplotypic, 0-9
    genotypic, one for each marker). This
    is expected to be a single column taken from the output of
    ibdhap.make.states.
}
  \item{data.type}{
    "h" : haplotypic data
    "g" : genotypic data (or hap data ran as genotypic)
    "r" : reduced data (output from running ibdhap.reduce.states and
    then ibdhap.make.states) 
}
  \item{col}{ Vector consisting of three colors. The colors represent no calls,
   any ibd shared, and no ibd shared, respectively. 
}
  \item{position}{ A position vector, with the same length as x describing the
    positions (in cM, M, or any other metric) of each marker. If
    positions is not included, "segment length" refers to number of SNPs
    making up a segment.
}
  \item{top}{top value of rectangles in barplot
}
  \item{bottom}{ bottom value of rectangles in barplot
}
  \item{density}{
    density of diagonal lines filling each segment in lines per sq
    inch. 
}
 
}

\author{

  MD Brown
}

\examples{

}
