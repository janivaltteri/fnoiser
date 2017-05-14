\name{fnoiser-package}
\alias{fnoiser-package}
\alias{fnoiser}
\docType{package}
\title{
  1/f Noise generation routine for R
}
\description{
  Generates coloured noise by summing sine waves.
}
\details{
  \tabular{ll}{
    Package: \tab fnoise\cr
    Type: \tab Package\cr
    Version: \tab 0.2.1\cr
    Date: \tab 2017-05-14\cr
    License: \tab GPL-2\cr
  }
  This package contains currently two functions: 'fn1(t,g)' and
  'fn2(t,g,r)'.
}
\author{
  Jani Anttila

  Maintainer: Jani Anttila <jani.anttila@helsinki.fi>
}
\references{
  Cohen et al. (1998) Proc. R. Soc. B. DOI: 10.1098/rspb.1998.0257
}
\keyword{ noise }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
  noise1 <- fn1(50,0.5)
  plot(noise1,type='l')
}
