\name{cumul.sum}
\alias{cumul.sum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
  cumulative.sums
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cumul.sum(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    %%     ~~Describe \code{x} here~~
    x is a vector
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
Returns a vector y where y[i]=sum(x[<i])
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
function(x){
y<-rep(0,length(x))
#cumulatively sums a vector x and returns the result
for(i in 1:length(x)){
y[i]<-sum(x[1:i])
  }
return(y)
  }
}
