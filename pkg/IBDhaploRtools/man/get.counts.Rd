\name{get.counts}
\alias{get.counts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
  get counts
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  subroutine
}
\usage{
get.counts(x, state.num)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    %%     ~~Describe \code{x} here~~
    vector
}
  \item{state.num}{
    %%     ~~Describe \code{state.num} here~~
    scalar
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function( x, state.num){

count.mat<-matrix(0, nrow=state.num, ncol=state.num)

   row.len<-length(x)

   for(i in 2:row.len){

      count.mat[x[i-1], x[i]] <- count.mat[x[i-1], x[i]] + 1

   }
return(count.mat)
  }
}
