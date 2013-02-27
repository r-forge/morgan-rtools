#### Example of use of transition matrix code
#### See also the commented code in transition.R


source("transition.R")
#examp <- read.table("matrix_15x15.txt", header = FALSE, skip = 1)
examp <- read.csv("tmat_15x15.txt", header = FALSE)

### Switch ibd state to label
ibd2label( s = c(1,1,2,2) )
ibd2label( s = c(1,2,3,3,4,4))

### Switch label back to ibd state
label2ibd( label = 4, n = 4)
label2ibd( label = 298, n = 6)

### Find maximum label for n chromosomes
maxlabel( n = 4)
maxlabel( n = 10)

### Change an invalid label to a valid one
rplc( c(2,3,4,5) )
rplc( c(3,1,2,7) )

### Transition Matrix
m1 <- createmat( n = 4, beta = 0.1)
m1$full    #with zero rows/cols for invalid labels
m1$reduced #with the zero rows/cols removed
m1$ind     #indices used to remove rows/cols

m2 <- createmat( n = 4, beta = 0.5)
m2$full
m2$reduced
m2$ind

### 4 Chromosome Transition Matrix with Standard Ordering
m3 <- fourchrom( beta = 0.15)
m3

m4 <- fourchrom( beta = 0.5)
m4


### testing if this code gives the same as examp?
examp
m3

