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

## The function is currently defined as
function(qibd.filename, dat.filename, output.filename){


   #glean information from the par file
   par.file<-read.table(dat.filename, fill =TRUE, colClasses = "numeric")

   n.snps<- par.file[3,1]  #number of markers
   n.sets<-par.file[2,1]   #number of set of haplotypes
   
   temp<-FALSE
   

   for( iset in 1:n.sets){

      #make the text above the (sets of) marker information in both files the same
      line1<-scan(qibd.filename, skip=((n.snps+2)*(iset-1)), nlines=1, what='raw')
      write(line1, file=output.filename, ncolumns=length(line1), append=temp)
      temp<-TRUE
  
      line2<-scan(qibd.filename, skip=((n.snps+2)*(iset-1)+1), nlines=1, what='raw')
      write(line2, file=output.filename, ncolumns=length(line2), append=TRUE)

      rm(line1,line2)
      
      #read the qibd.filename, just the information for the markers in "iset"
      ibd.dat<-scan(qibd.filename,what = "numeric",skip = (2*(iset)+(n.snps)*(iset-1)), nlines = n.snps)
      ibd.dat<-t(matrix(as.numeric(unlist(ibd.dat)), nrow=17))




     reduced.dat <- t( apply(ibd.dat,1,sumcol) )

     reduced.dat <- cbind( ibd.dat[,1:2], reduced.dat)

     write.table( reduced.dat,file=output.filename,  append = TRUE, col.names=FALSE, row.names=FALSE)
     }
     
     return(paste("file", output.filename, "written in", getwd()) )
  }
}
