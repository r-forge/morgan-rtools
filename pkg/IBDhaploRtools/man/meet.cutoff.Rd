\name{meet.cutoff}
\alias{meet.cutoff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  meet cutoff
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
meet.cutoff(num, cutoff)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{num}{
    number to be tested
%%     ~~Describe \code{num} here~~
}
\item{cutoff}{
  cutoff to be tested against
%%     ~~Describe \code{cutoff} here~~
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
function( num, cutoff) {
  #returns a 0 (false) or 1 (true) to tell if num is larger than cutoff
	as.numeric(num > cutoff);
  }
}
