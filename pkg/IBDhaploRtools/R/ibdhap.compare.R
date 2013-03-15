

###
#Main Function

#function to compare inferred results to simulated truths

ibdhap.compare<-function( states.dat, simu.dat, data.type = c("h", "g","r")){

#define parameters based on haplotype/genotype/or reduced data
  if(length(data.type)>1){ stop("data.type improperly indicated")}
  else if( is.element("h", data.type)){ no.ibd.ind = 15}
  else if( is.element("g", data.type)){no.ibd.ind = 9}
  else if( is.element("r", data.type)){no.ibd.ind = 4}
  else{ stop("data.type improperly indicated")}  

  n.pairs<-ncol(states.dat)
  n.snps<-nrow(states.dat)
  
  prop.correct<- sum(states.dat==simu.dat)/(n.pairs*n.snps)

  #inferrring ibd when there is none
  false.positives<-sum( simu.dat==no.ibd.ind & states.dat!=simu.dat & states.dat!=0)/ (n.pairs*n.snps)

  #inferring no ibd when there is some

  false.negs<-sum( simu.dat!=no.ibd.ind & states.dat==no.ibd.ind)/(n.pairs*n.snps)

  #no calls

  no.calls<-sum(states.dat==0)/(n.pairs*n.snps)

  return(list(prop.correct=prop.correct, false.positives=false.positives,false.negatives= false.negs, no.calls=no.calls))
}
