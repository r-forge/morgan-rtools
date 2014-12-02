\name{State Vectors}
\alias{allVec}
\alias{vec2label}
\alias{label2vec}
\alias{rplc}
\title{
IBD state vectors
}
\description{
  
These vectors are complete descriptions of the IBD state. There is one
vector element per gamete. The first gamete is always numbered one, with
subsequent gametes being given the same number if they are IBD with a
previous gamete, or a new number if they are not.

For example in (1,2,2,3), the second and third gametes are IBD.

Note that label starts at ZERO.
}
\usage{
allVec(ngam)

vec2label(vec, renum=TRUE)

label2vec(label, ngam) 

rplc(vec)
}

\arguments{
  \item{ngam}{Number of gametes if a label is supplied, otherwise number
    of gametes will be determined by the length of the vector.}
  \item{vec}{IBD state description vector. Must be a single vector.}
  \item{label}{IBD label. Must be an integer.}
  \item{renum}{Do you want to run rplc on the vector before converting to label.}
}
\value{

For allVec, a matrix of the vectors for each label is produced, with a
vector on each row and a row for each label. The labels are the row
names. 

For vec2label and label2vec an integer label or vector is prouduced
respectively.

  If renum=TRUE, vec2label function calls rplc, so any vector supplied
  will be renumbered. If (1,1,1,3) is given, it will be renumbered
  (1,1,1,2) and given the label 1, not 2. Some labels that refer to the
  same state as another label cannot be reached in this case.
 

rplc returns a vecor of the same length as the input, renumbered.

}
\author{
Fiona Grimson
}
\examples{

## Vectors for all labels
allVec( ngam = 4 )

## Convert vector to label, with renumbering or without
vec2label( c(1,1,1,3), renum = TRUE )
vec2label( c(1,1,1,2), renum = TRUE ) 

## Convert label to vector
label2vec( 1, ngam = 4 )
label2vec( 2, ngam = 4 )
label2vec( 2, ngam = 6 ) 

## renumbering
rplc( c( 1,1,1,3) )
rplc( c(5,1,5,6) ) 

}
\keyword{ Vec }
\keyword{ Vector }
