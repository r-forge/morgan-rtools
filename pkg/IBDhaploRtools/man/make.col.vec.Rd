\name{make.col.vec}
\alias{make.col.vec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  make color vector
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make.col.vec(x, colors)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    vector of 0, 1, 2's
%%     ~~Describe \code{x} here~~
}
\item{colors}{
  vector of three colors.
  
%%     ~~Describe \code{colors} here~~
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
function(x, colors){
#input a vector (x) of 0, 1, 2 outputs a vector of "col1", "col2" , "col3"
# as defined by the vector of three colors (col)
result<-NULL
for( i in 1:length(x)){result = c(result, colors[x+1])}
return(result)
  }
}
