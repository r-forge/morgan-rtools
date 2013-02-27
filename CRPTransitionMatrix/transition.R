### Create a transition matrix for IBD states among n chromosomes

# ******************************************************************#
### Functions for switching from label to state
###    these functions are based on
###    /mounts/gaw16/MORGAN-3_EXPORT_Feb2011/Markers/index_ibd.c
###    functions index_ibd() and invert_label()
# ******************************************************************#

ibd2label <- function(s){
  ## go from ibd vector to label
  ## s =vector of ibd states e.g. s = c(1,2,3,4)

  n = length(s)
  label = 0
  for( i in 0:(n-1)){
    label = label*i + (s[i+1] - 1)
  }
  return(label)
}

label2ibd <- function(label,n){
  ## go from label to ibd vector

  ##initialise
  cc = label
  flag = 0
  s = numeric(n)
  for( i in (n-1):0 ){
    ##f cc = 0 break (flag=0)
    if (cc==0) break 
    ##find gg
    gg = cc - (cc %/% i)*i 
    ##continue if non-zero or cc < CDgg
    if (gg !=0 | cc< maxlabel(i+1) ){ 
      s[i+1] = gg+1
      cc = cc %/% i
    }else{
      flag =1; break
    }
    i=i-1 
  }
  ##fill rest of s
  if(flag == 0) s[1:(i+1)] = 1
  if(flag == 1) s[1:(i+1)] = 1:(i+1)
  return(s)
}


maxlabel <- function(n){
  ## find maximum number of labels for n chromosome partitions
  ibd2label( s = 1:n)
}

# ******************************************************************#
### Functions to calculate possible transitions from s
# ******************************************************************#

rplc <- function(x){
  
  ## function to turn a state vector to the minimum numbering, e.g. s
  ## = c(2,3,1,1) to s = c(1,2,3,3).

  ## This ensures that when going from state -> label, the label is
  ## always the first label representing that state.
  
  y <- numeric( length(x))
  xs <- unique(x)
  for( i in 1:length(xs)){
    y <- replace(y, x==xs[i], i)
  }
  return(y)
}


transitions <- function(s, beta){

  ## given a state s (e.g. s=c(1,2,3,3)) and a beta (prob of joining a
  ## group size), this function returns a matrix with three columns:
  ## the labels of all the possible new states; the frequency of
  ## transitions to each state; and the probability of getting that
  ## state.
  
  ## Note that some states will be redundant( more than one label per
  ## state) but this function will calculate transitions to only the
  ## label that is the first representation of that state.

  ## s = state coming from
  ## beta = prob. that joins a group
  n = length(s) #No of chromosomes
  d = maxlabel(n) #maximum label

  gps = unique(s)
  gpsz = tabulate(s)

  newguy = c(unique(s), n+1)
  probs = c( beta*gpsz, 1- beta )
 
  ##calculate all the new transitions and their rate
  t2 <- NULL
  probs2 <- NULL
  for( i in 1:length(newguy)){##for each arrival
    
    for( j in 1:n ){##for each exit of existing chrom
      trans <- s
      trans[j] <- newguy[i]
      t2 <- rbind(t2, trans)
    }
    t2 <- rbind(t2, s)##new guy leaves straight away
    ##probability of this event
    probs2 <- c( probs2, rep(probs[i], n+1))
  }

  ##change to std notation
  t3 <- t(apply(t2, 1, rplc))

  ##change to labels
  t4 <- apply(t3, 1, ibd2label)
  
  ##transition frequencies
  t5 <- matrix( ncol = 3, nrow = maxlabel(n)+1)
  colnames(t5) = c("label", "freq", "prob")
  for( i in 0:maxlabel(n)){
    labi = which(t4== i)
    nlabi = length(labi)
    probi = sum( probs2[labi] )

    t5[i+1, ] = c(i, nlabi, probi)
  }
  return(t5)
}
  
# ******************************************************************#
### Functions to calculate transition prob matrix (Q)
# ******************************************************************#
library(mgcv)#for uniquecombs function


createmat <- function(n, beta){
  ## function to create the transition matrix from the number of
  ## chromosomes and the beta probability
  
  
  ## Inputs:
  ##  n = no of chromosomes
  ##  beta = prob of joining new group
  ## Outputs: returns a list of
  ##  full = matrix with rows and cols of zeros for "dead" labels
  ##  reduced = matrix with zero rows and cols removed
  ##  ind = vector to go from reduced to full,
  ##     e.g. to recover the full matrix from the reduced run
  ##     full = red[ind,ind]; full[ind,]=0; full[,ind]=0
  
  d = maxlabel(n)

  smat = matrix(ncol = n, nrow = d+1) #may have multiple labels for same state
  for( i in 1:(d+1)) smat[i,] = label2ibd( i-1, n)
  smat1 = t(apply(smat, 1, rplc))#change to std numbering
  smat2 = uniquecombs(smat1)#reduce to unique rows

  ns = nrow(smat2)#number unique states
  ind <- attr(smat2, "index")
  
  tmat = matrix(ncol = d+1, nrow = ns)##cols for all possible labels, rows for unique
  for( i in 1:ns) tmat[i,]= transitions(smat2[i,], beta)[,3]

  ##remove excess rows cols
  zcols = which( colSums(tmat) == 0 )#zero columns
  redtmat = tmat[,-zcols] #reduced down
  dt = dim(redtmat)[1]

  ##fix the diagonals
  redtmat0 = redtmat - diag( diag(redtmat), ncol = dt, nrow = dt)
  redtmat0[cbind(1:dt, 1:dt)] = - apply( redtmat0, 1, sum)

  ##filled out with zero rows and cols
  fulltmat = redtmat0[ind,ind];
  fulltmat[zcols,]=0; fulltmat[,zcols]=0 #filled out

  return(list(full = fulltmat, reduced = redtmat0, ind = ind))
}

# ******************************************************************#
### Function to produce 4 chromosome matrix with standard order
# ******************************************************************#

fourchrom <- function(beta){
  ## function that, given beta, returns the 15*15 transition matrix
  ## with standard ordering

  fullmat <- createmat(4, beta)$full
  order <- c(0, 4, 1, 3, 5, 6, 10, 14, 7, 9, 8, 12, 11, 13, 15)+1
  stdmat <- fullmat[order, order]
  return(stdmat)
}


#### TO DO:
# check the 4 chrom matrix against what Chris emails.
