\name{ibdhap.compare}
\alias{ibdhap.compare}

\title{
  ibdhap compare
}
\description{Compares inferred ibd state with simulated "true" states. Calculates the proportion of markers called in the correct state, false positives (i.e. inferring ibd shared when none is shared), false negatives (i.e. inferring no ibd shared when ibd sharing is present) and the proportion of no calls. 
}

\usage{
ibdhap.summary(states.dat, simu.data, data.type=c("h", "g","r")
}

\arguments{
  \item{states.dat}{ The data.frame created from running ibdhap.make.states on ibd_haplo
    output. 
}

\item{simu.dat}{ The data frame of same dimension as states.dat, but with true (probably simulated) ibd states. 
}

\item{data.type}{
    "h" : haplotypic data
    "g" : genotypic data (or hap data ran as genotypic)
    "r" : reduced data (output from running ibdhap.reduce.states and
    then ibdhap.make.states) 
}



}
  


\value{

Returns a list consisting of:
\item{prop.correct}{ proportion of markers called in exactly the correct state}
\item{false.positive}{proportion of markers in which ibdhaplo inferred ibd sharing when none is present}
\item{false.negative}{ proportion of markers in which ibdhaplo inferred no ibd sharing when sharing is present}
\item{no.calls}{proportion of markers that did not meet the cutoff selected when using ibdhap.make.states}
}

\author{
  MD Brown
}

\examples{

}




