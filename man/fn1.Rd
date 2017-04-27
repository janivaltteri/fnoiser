\name{fn1}
\alias{fnoisegenerator}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Coloured noise generating routine
}
\description{
Generates 1/f noise.
}
\usage{
fn1(t, g)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{
    The length of time series of noise.
  }
  \item{g}{
    The colour (exponent) of noise.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  good
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
Jani Anttila
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

joo1 <- fn1(100,1.0)
joo2 <- fn1(100,0.0)
joo3 <- fn1(100,-1.0)

plot(joo1,type='l',col='red',ylim=c(-4,4),lwd=3)
lines(joo2,lwd=2)
lines(joo3,col='blue')

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
