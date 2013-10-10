\name{bss.components}
\alias{bss.components}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to Extract Estimated Sources from an Object of Class bss
}
\description{
Extracts the sources estimated by an bss method.
}
\usage{
bss.components(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{object of class bss}
}

\author{Klaus Nordhausen}

\examples{
A<- matrix(rnorm(9),3,3)
s1 <- arima.sim(list(ar=c(0.3,0.6)),1000)
s2 <- arima.sim(list(ma=c(-0.3,0.3)),1000)
s3 <- arima.sim(list(ar=c(-0.8,0.1)),1000)

S <- cbind(s1,s2,s3)
X <- S \%*\% t(A)

res1<-AMUSE(X)
head(bss.components(res1))
colMeans(bss.components(res1))
cov(bss.components(res1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
