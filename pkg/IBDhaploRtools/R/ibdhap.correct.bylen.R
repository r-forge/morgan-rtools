ibdhap.correct.bylen<-
function( states.dat, simu.dat,data.type = c("h", "g","r"),ibd.only=TRUE, position = NA){


n.sets<-ncol(states.dat)
n.snp<-nrow(states.dat)
if(is.element(position[1],NA)){position <- 1:(n.snp) }
result=NULL

for( iset in 1:n.sets) #for each simulation
{

  #define parameters based on haplotype/genotype/or reduced data
if(length(data.type)>1){ stop("data.type improperly indicated")}
else if( is.element("h", data.type)){ no.ibd.ind = 15}
else if( is.element("g", data.type)){no.ibd.ind = 9}
else if( is.element("r", data.type)){no.ibd.ind = 4}
else{ stop("data.type improperly indicated")}


  
  seg.lengths<-c(1)	
  
  temp.sim<-as.numeric(unlist(simu.dat[,iset])) #take the simulation in question
  temp.estimates<-as.numeric(unlist(states.dat[,iset]))
  
  
  for(snp.index in 2:n.snp)
    {
      prev.val<-temp.sim[snp.index-1]
      val=temp.sim[snp.index]
      
      if( prev.val!=val) 
        {
          seg.lengths=c(seg.lengths, snp.index)
        }
      
    }	
  if(length(seg.lengths)==0){seg.lengths=c(1, n.snp)}
  if(seg.lengths[length(seg.lengths)]!= n.snp){seg.lengths=c(seg.lengths, n.snp)}
  
  
                                        #now for each segment of ibd we record its length in cM and the proportion of ibd that we called correctly
                                        #store this data in "result"
  
  
  for( seg.index in 1:(length(seg.lengths)-1))
    {
      
      seg.beg<-seg.lengths[seg.index]
      seg.end<-seg.lengths[seg.index+1]-1		
      
      true.ibd.val<-temp.sim[seg.beg]
      
      if(true.ibd.val==no.ibd.ind){
        proportion.right<- sum(temp.estimates[seg.beg:seg.end]==true.ibd.val)/(seg.end-seg.beg+1)
      }
      else{

        temp<-temp.estimates[seg.beg:seg.end]
        proportion.right<-1-(sum(temp==no.ibd.ind)+sum(temp==0))/(seg.end-seg.beg+1)
        
      }
      
      length.cM<-position[seg.end+1]-position[seg.beg]
      result=rbind(result, c(length.cM, proportion.right, true.ibd.val))
      
      
    }
  

}

#take out the segments that are no ibd shared
if(ibd.only){
  result=result[result[,3]<no.ibd.ind,]
}
  
result<-as.data.frame(result)
names(result)=c("seg.len", "prop.corr", "ibd.state")
return(result)
}

