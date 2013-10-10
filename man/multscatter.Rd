\name{multscatter}
\alias{multscatter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to Compute Several Scatter Matrices for the Same Data 
}
\description{
The function can be used to compute several scatter matrices for the same data.
}
\usage{
multscatter(scatterlist, X, toshape = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scatterlist}{
a vector with the names of the scatter matrices to be computed. Note that each of these functions should only return a matrix of size p times p.
}
  \item{X}{
the n times p data matrix for which the scatter should be computed.
}
  \item{toshape}{
logical, whether scatter matrices should be converted to shape matrices. If TRUE, all matrices will have determinant 1.
}
}
\details{
It is important that the functions do not need any additional imput and that they return only the p times p scatter matrix. Hence it might be sometimes necessary
to write wrappers for some of the functions. See examples.
}
\value{
An array of dimension c(p,p,k) where k is the number of scatter matrices.
}

\author{Klaus Nordhausen}

\examples{
# example requires the packages ICS and ICSNP
library(ICSNP)
X <- cbind(rexp(1000), rt(1000,6), runif(1000)) 

my.tM1 <- function(X,df=1) tM(X,)$V
my.tM2 <- function(X,df=2) tM(X,)$V

multscatter(c("cov","cov4","HP1.shape","my.tM1", "my.tM2"), X)
multscatter(c("cov","cov4","HP1.shape","my.tM1", "my.tM2"), X, toshape=FALSE)
}
\keyword{ multivariate }
