
source("transition.R")


integer2bin <-function(x){
  ## a function from the internet that changes integers to binary
  ## numbers. Returns a vector of zeros and ones (but vector may need
  ## to have zeros added to the front to bring it up to the correct no
  ## of SNPs)
  if(x==0){ return(0) }
  b=2
  xi <- as.integer(x)
  if(any(is.na(xi) | ((x-xi)!=0)))
    print(list(ERROR="x not integer", x=x))
  N <- length(x)
  xMax <- max(x)	
  ndigits <- (floor(logb(xMax, base=2))+1)
  Base.b <- array(NA, dim=c(N, ndigits))
  for(i in 1:ndigits){#i <- 1
    Base.b[, ndigits-i+1] <- (x %% b)
    x <- (x %/% b)
  }
  if(N ==1) Base.b[1, ] else Base.b
}

bin2integer <- function(x){
  ## only for testing
  sum(x * 2^(rev(seq_along(x)) - 1))
}
maxSNP <- function(n){
  ## only for testing
  bin2integer( rep(1, n))
}

pen <- function(SNP.integer, freq, eps, nSNP){
  ## A function that will calculate the penetrance proabilities for
  ## all the different IBD states
  ## i.e. P(observe input SNPs | IBD state)

  ## INPUTS :
  ## SNP = integer representing SNP values
  ## freq = vector length 2 = (prob snp=1, prob snp=2)
  ## eps = probability of an error
  ## nSNP = how many SNPs there should be

  ## OUTPUTS :
  ## list of two elements
  ## 1. matrix of all the ibd states
  ## 2. the corresponding probabilities
  
  SNP.bin <- integer2bin(SNP.integer)##change integer to biary
  rem <- nSNP - length(SNP.bin) ##how many zeros to put in front
  SNP <- c( rep(0, times=rem), SNP.bin) ## a SNP vector of the correct length
  
  ## calculate all the states
  mlab <- maxlabel( n = nSNP )##get all possible states
  states <- sapply( 0:mlab, label2ibd, n = nSNP) 
  states <- t(apply(states, 2, rplc)) ##remove invalid and duplicate states
  states <- uniquecombs(states)

  ##for each state, get the probability
  nStates = nrow(states)
  probs <- rep(0, nStates)
  for( i in 1:nStates){

    classes = unique( states[i,] ) ##get the classes
    cprobs = rep(0, length = length(classes))

    for( j in classes ){ ##for each class, get the class probability
      
      snp = SNP[ states[i,]==j ] ##the snps in the class
      p1 = eps^sum(snp!=0)*(1-eps)^sum(snp==0)   ##prob data given true state is 1
      p2 = eps^sum(snp!=1)*(1-eps)^sum(snp==1)   ##prob data given true state is 2  
      cprobs[j] = freq[1]*p1 + freq[2]*p2 ##the class probability
      
    }
  
    probs[i] = prod( cprobs ) ##sum over the class probs to get state prob
  }
  return( list( states = states, probs = probs) )
}

#### TESTING ####

## ## Four chromosome model from chris - should be the same

## ## the genotypes (integers)
## geno <- 0:15

## mat26 <- matrix(ncol = 15, nrow = 16)#15 states and 16 genotypes

## for( i in 1:16){
##   mat26[i,] = pen( SNP.integer = geno[i],
##          freq = c(0.26, 0.74), eps = 0, nSNP=4)$probs
## }

## mat74 <- matrix(ncol = 15, nrow = 16)#15 states and 16 genotypes

## for( i in 1:16){
##   mat74[i,] = pen( SNP.integer = geno[i],
##        freq = c(0.26, 0.74), eps = 1, nSNP=4)$probs
## }

## Six Chroms
## n=6
## geno <- 0:maxSNP(n)

## mat <- matrix( ncol= 203, nrow = maxSNP(n)+1)

## for( i in 1:(maxSNP(n)+1)){
##   mat[i,] <- pen( geno[i], freq = c(0.4, 0.6), eps =0.1, nSNP=6)$probs
## }
