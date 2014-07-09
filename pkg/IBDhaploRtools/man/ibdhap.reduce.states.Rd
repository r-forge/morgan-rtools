\name{ibdhap.reduce.states}
\alias{ibdhap.reduce.states}

\title{
  ibdhap reduce states
}
\description{
  This function reduces the columns of a qibd.out file created by
  ibd_haplo by summing probabilities over certain columns. When reducing
  an output file that was run on haplotypes, this script will sum over
  columns 1-8 (other), 9-10 (2 pairs chrs. ibd), 11-14 (one pair of
  chrs. ibd),and  15 (no ibd). When reducing a file from genotypic date,
  the corresponding columns are summed over so that they reflect the
  same values ( other, one pair ibd, two pairs ibd, not ibd). A file will then be output with these
  probabilites displayed just like the original ibd_haplo output qibd
  file. 
}

\usage{
ibdhap.reduce.states(qibd.filename, dat.filename, output.filename)
}

\arguments{
  \item{qibd.filename}{ The filename of the qibd.out file the is produced from running
    ibd_haplo. This is to be input as a character string. In the
    examples, this file is called "qibd_g.gold". 
}
  \item{dat.filename}{ The filename of the .dat file used to run ibd_haplo. This is the
  complicated parameter file that consists of three lines of values /
  indicators that tells ibd_haplo what to do. In the examples, this
  files is called "compu_4haps.dat".
  
}
  \item{output.filename}{The name of the file that this function will print to. In the
    examples, this filename is "test.file.txt". 
}
}

\value{
A txt document the looks exactly like the qibd file input, but the ibd
state probabilities per marker are summed over in the manner described
above. 
}


\author{
MD Brown
}

\examples{
## See vignette("IBDhaploRtools_tutorial" )
}
