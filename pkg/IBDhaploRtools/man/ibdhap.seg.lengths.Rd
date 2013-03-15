\name{ibdhap.seg.lengths}
\alias{ibdhap.seg.lengths}

\title{ ibdhap.seg.lengths
}

\description{
  Given the ibd states from a set of haplotypes/pair of genotypes (taken
from a column of the output of ibdhap.make.states), this function creates a
data.frame consisting of all segments of differing ibd state, paired
with their respective length.
}
\usage{
ibdhap.seg.lengths(x, position = NA)
}

\arguments{
  \item{x}{A vector of ibd states (with values 0 - 15 for haplotypic, 0-9
    genotypic, one for each marker). This
    is expected to be a single column taken from the output of
    ibdhap.make.states.  
}
  \item{position}{ A position vector, with the same length as x describing the
    positions (in cM, M, or any other metric) of each marker. If
    positions is not included, "segment length" refers to number of SNPs
    making up a segment. 
}
}


\value{
A data.frame with 2 columns and nrow = nrow(x)(number of markers). column1 is the integer value corresponding
to ibd state, column2 is the length of the segment for that state as
measured by the positions vector. 
}

\author{
  MD Brown
}

\examples{

## The function is currently defined as
function( x, position=NA ){
#x is a single column from ibdhap.make.states output
## NOT INCLUDING THE haplotype/genotype names!,
# thus it is the ibd.states from one pair of individuals (genotypes)
# or set of 4 haplotypes


n.marker<- length(x)  #number of markers

#if positions are given, we use them, otherwise "length" refers to
# number of SNPs
if(is.element(position,NA)){position <- 1:(n.marker) }



# obtain a vector of ibd state change points --index where ibd state changes
	change.points<-c(1)

	for(imarker in 2:n.marker)
	{
		prev.val<-x[imarker-1]
		val <- x[imarker]

		if( prev.val!=val)
		{
		   change.points=c(change.points, imarker)
		}

	}

	# tidy up the end of change.points
  if(change.points[length(change.points)]!= n.marker){change.points=c(change.points, n.marker)}

	change.points.pos<-position[change.points]
	seg.lengths<-diff(change.points.pos)
	ibd.state<-x[change.points[1:length(seg.lengths)] ]

  return( as.data.frame(cbind(ibd.state = ibd.state, seg.lengths = seg.lengths)))

	}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
