\name{sumcol}
\alias{sumcol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
  sum over specific columns
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sumcol(rowdat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rowdat}{
    a vector of length 15
%%     ~~Describe \code{rowdat} here~~
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
a vector of length 4, since some elements were summed over
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
function(rowdat){

  output<-c(sum(rowdat[3:10]), sum(rowdat[11:12]), sum(rowdat[13:16]), sum(rowdat[17]))
  return(t(output))
  }
}
