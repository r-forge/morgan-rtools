\name{ibdhap.correct.bylen}
\alias{ibdhap.correct.bylen}

\title{
  create data frame comparing inference vs simulated data by segment length
}
\description{
 
   If simulated data(true ibd state) is
available, this function creates a dataframe of how well ibd_haplo infers segments based on segment length.
}
\usage{
ibdhap.infer.vs.simu.plot(states.dat, simu.dat, data.type=c("h", "g","r"),position = NA, ibd.only=TRUE)
}

\arguments{
  \item{states.dat}{ The data.frame created from running ibdhap.make.states on ibd_haplo
    output. 
}
  \item{simu.dat}{ The data.frame created from running ibdhap.make.states on
    simu.c/ibd_states.c output
}
  \item{data.type}{
    "h" : haplotypic data
    "g" : genotypic data (or hap data ran as genotypic)
    "r" : reduced data (output from running ibdhap.reduce.states and
    then ibdhap.make.states) 
}
  \item{position}{A position vector, with the same length as nrow(states.dat) describing the
    positions (in cM, M, or any other metric) of each marker. If
    positions is not included, "segment length" refers to number of SNPs
    making up a segment.
}
  \item{ibd.only}{ When true, the result consists of only segments that are some form of ibd shared. 
}
}
\value{
A dataframe where each row represents a segment of true (simulated) ibd
(shared or not shared).  The first column is the length of the segment, the second column is the proportion snps called correctly within the true segment, and the third column is the simulated (true) ibd state. 
}
\author{
 MD Brown 
}

\examples{


}

