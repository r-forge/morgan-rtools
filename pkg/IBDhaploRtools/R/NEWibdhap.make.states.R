## modified 3 Feb 2013

## example inputs
qibd.filename <- "qibd_ten_ss.gold"
dat.filename <- "ids_ten_ss.gold"
## cutoff <- 0.8

## source("meet.cutoff.R")
## source("return.ibd.val.R")

## qibd.gold <- ibdhap.make.states(qibd.filename = "qibd_phased_2011.gold", dat.filename = "ids_phased_2011.gold", cutoff = 0.9)

ibdhap.make.states <-function( qibd.filename, dat.filename, cutoff = .8){

  ## information from dat file
  n.col.par <- max(count.fields(dat.filename)) ## max number of columns in the input file
  par.file <- read.table(dat.filename, fill =TRUE, colClasses = "numeric", col.names=1:n.col.par)
  name.sets <- par.file[,1] ## names of the sets
  
  ## information form qibd file
  n.col.rawdat <- count.fields(qibd.filename) ## number of coulmns of each row in rawdat file
  rawdat <- read.table(qibd.filename, fill =TRUE, colClasses = "numeric", col.names=1:max(n.col.rawdat))
  print("Data read in successfully")
  n.snps <- dim( rawdat[name.sets==name.sets[1],] )[1] ## number of snps

  
  ##geno.indicator, n.names etc were defined in original but in the
  ##end not used.  intended to get the names of the
  ##genotypes/haplotypes.

  ## empty matrix that will eventually be the output (assuming all sets have the same number of markers
  ibd.states <- matrix(NA, ncol=length(name.sets), nrow=n.snps)

  count <- 1
  for( iset in name.sets){
    
    ## subset just the data for the markers in the set.
    xind <- which( rawdat[,1]==iset)
    n.col.i <- n.col.rawdat[xind][1] ## number of columns in the set i
    ibd.dat <- t(rawdat[xind, 2:n.col.i])

    ## truncate based on "cutoff" : put 1 if cutoff is met, 0 otherwise
    ibd.trunc <- apply(ibd.dat[ 3:(n.col.i-1), ], c(1,2), meet.cutoff, cutoff)

    ## for each column (state probs for a marker) list ibd.state, or 0 if no state > cutoff
    ibd.states.temp<-apply(ibd.trunc, 2, return.ibd.val)

    ## add the marker names -- was not implemented in original

    ## save in the final matrix
    ibd.states[,count]=ibd.states.temp
    ##print(paste("pair", iset, "of", n.sets, "complete...", sep=" "))

    count <- count + 1
  }
  return(as.data.frame(ibd.states))   
}

  ############################################################################
    ##                                     #glean information from the par file
    ## par.file<-read.table(dat.filename, fill =TRUE, colClasses = "numeric")
    
    ## n.snps<- par.file[3,1]  #number of markers
    ## n.sets<-par.file[2,1]   #number of set of haplotypes
    ## geno.indicator<-par.file[1,3] #data input as genotypic or not
    ##                                     #information to expect individual names (genotypic data) 
    ##                                     #or haplotype names (hap data)
    ## if(geno.indicator){
    ##   c.classes<-c(rep("character", 2), rep("numeric", 2))
    ##   names.col<-c(3,4)
    ##   n.names<-2
    ## }
    ## else{ c.classes<-c("character", rep("numeric",4))
    ##       names.col<-(2:5)
    ##       n.names<-4
    ##     }
  ###############################################################################

    ####################################################################################
    ##                                     #find how many columns to plan for for different data (genotype, haplotype, reduced)
    ## n.col<-length(read.table(qibd.filename, skip = 2, nrow=1))
    
    ##                                     #define the output dataframe "ibd.states"
    ## ibd.states<-matrix( 0, ncol=n.sets, nrow =(n.snps))
    ## rawdat<-read.table(qibd.filename, fill=TRUE, header=FALSE)
    ## print("Data read in successfully")
    ####################################################################################
  
  ###############################################################################################
   ##  for(iset in 1:n.sets){
      
   ##                                      #read the qibd.filename, just the information for the markers in "iset"
   ##    jumpind<-2*(iset)+(n.snps)*(iset-1)+1
   ##    ibd.dat<-t(rawdat[jumpind:(jumpind+n.snps-1),])
   ##    ibd.dat<-matrix(as.numeric(unlist(ibd.dat)), nrow=(n.col))
      
   ##                                      #get haplotype names in this set of haplotypes
   ##                                      # names<-read.table(qibd.filename, colClasses = c.classes, skip = (2*(iset)+(n.snps)*(iset-1)-1), nrow=1) 
      
   ##                                      #truncate the data based on "cutoff", leave a 1 if the cutoff is met, 0 othwse
   ##    ibd.trunc<-apply(ibd.dat[ 3:(n.col), ], c(1,2), meet.cutoff, cutoff) 
      
   ##                                      #now for each column (set of state probabilities for a marker), list the ibd.state
   ##                                      #(or 0 if no state probablity is > cutoff
   ##    ibd.states.temp<-apply(ibd.trunc, 2, return.ibd.val) 
      
   ##                                      #add the marker names
   ##                                      #ibd.states.temp<-as.numeric(unlist(c(names[names.col], ibd.states.temp)))
      
   ## #add to the final product
   ##    ibd.states[,iset]=ibd.states.temp
   ##    print(paste("pair", iset, "of", n.sets, "complete...", sep=" "))
   ##  }
    
  ##   return(as.data.frame(ibd.states))
    
  ## }

######################################################################################################
