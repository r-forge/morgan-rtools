\name{return.ibd.val}
\alias{return.ibd.val}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
  return ibd val
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  returns the index of a 1, if there is one in this vector
}
\usage{
return.ibd.val(col.dat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{col.dat}{
    %%     ~~Describe \code{col.dat} here~~
    a vector of all zeros with at most a single 1
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


## The function is currently defined as
function(col.dat){
   #expects a vector of zeros with at most one 1
   #returns the index of where the 1 is, or zero if 1
   # is not present in col.dat    
return(which.max(col.dat)*max(col.dat))
  }
}
