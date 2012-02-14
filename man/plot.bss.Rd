\name{plot.bss}
\alias{plot.bss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plotting an Object of Class bss}
\description{
Plots the estimated sources resulting from an bss method. If the bss method is based on second order assumptions and returned the sources as a time series object it will plot the sources 
using \code{plot.ts}, otherwise it will plot a scatter plot matrix using \code{pairs}. 
}
\usage{
\method{plot}{bss}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{object of class bss.}
  \item{\dots}{further arguments to be passed to or from methods.}
}

\author{Klaus Nordhausen, \email{klaus.nordhausen@uta.fi}}

\seealso{
\code{\link{plot.ts}}, \code{\link{pairs}}
}
\examples{
A<- matrix(rnorm(9),3,3)
s1 <- arima.sim(list(ar=c(0.3,0.6)),1000)
s2 <- arima.sim(list(ma=c(-0.3,0.3)),1000)
s3 <- arima.sim(list(ar=c(-0.8,0.1)),1000)

S <- cbind(s1,s2,s3)
X <- S \%*\% t(A)

res1 <- AMUSE(X)
plot(res1)
# not so useful:
plot(res1, plot.type = "single", col=1:3)

# not meaningful for this data
res2 <- JADE(X)
plot(res2)
}
\keyword{methods}
