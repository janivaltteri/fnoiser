\name{fnn}
\alias{fnoisegenerator-n}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Coloured noise generating routine
}
\description{
Generates n 1/f noises with a given correlation.
}
\usage{
fn1(time, gamma, correlation)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
    The number of noises.
  }
  \item{t}{
    The length of time series of noise.
  }
  \item{g}{
    The colour exponent of noise.
  }
  \item{r}{
    Correlation between the two generated noises.
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
  Cohen et al. (1998) Proc. R. Soc. B. DOI: 10.1098/rspb.1998.0257
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

nsn <- fnn(3,100,1.0,0.5)
plot(nsn[,1],type='l')
lines(nsn[,2])
lines(nsn[,3])

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
