\name{coef.bss}
\alias{coef.bss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Coefficients of a bss Object
}
\description{
Extracts the estimated unmixing matrix from an object of class bss.
}
\usage{
\method{coef}{bss}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{object of class bss.}
  \item{\dots}{further arguments to be passed to or from methods.}
}

\author{Klaus Nordhausen, \email{klaus.nordhausen@uta.fi}}

\examples{
A<- matrix(rnorm(9),3,3)
s1 <- arima.sim(list(ar=c(0.3,0.6)),1000)
s2 <- arima.sim(list(ma=c(-0.3,0.3)),1000)
s3 <- arima.sim(list(ar=c(-0.8,0.1)),1000)

S <- cbind(s1,s2,s3)
X <- S \%*\% t(A)

res1<-AMUSE(X)
coef(res1)
coef(res1) \%*\% A # should be a matrix with one dominant element in each row and column
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
