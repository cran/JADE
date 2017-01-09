\name{NSS.SD}
\alias{NSS.SD}
\alias{NSS.SD.default}
\alias{NSS.SD.ts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
NSS.SD Method for Nonstationary Blind Source Separation
}
\description{
The NSS.SD method for nonstationary blind source separation. The function estimates
the unmixing matrix in a nonstationary source separation model by simultaneously diagonalizing
two covariance matrices computed for different time intervals.
}
\usage{
NSS.SD(X, ...)

\method{NSS.SD}{default}(X, n.cut=NULL, ...)
\method{NSS.SD}{ts}(X, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a numeric matrix or a multivariate time series object of class \code{\link{ts}}. Missing values are not allowed.}
  \item{n.cut}{either an integer between 1 and nrow(X) or an vector of length 3 of the form c(1,n.cut,nrow(X)) to specify where to split the time series. If NULL, then
  c(1,floor(nrow(X)/2),nrow(X)) is used.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{
The model assumes that the mean of the p-variate time series is constant but the variances change over time.
}
\value{
 A list with class 'bss' containing the following components:
  \item{W}{estimated unmixing matrix.}
  \item{EV}{eigenvalues from the eigenvalue-eigenvector decomposition.}
  \item{n.cut}{specifying the intervals where data is split}
  \item{S}{estimated sources as time series objected standardized to have mean 0 and that the sources in the first interval are 1.}
}


\references{
\cite{Choi S. and Cichocki A. (2000), Blind separation of nonstationary sources in noisy mixtures, 
            \emph{Electronics Letters}, 36, 848--849.}

\cite{Choi S. and Cichocki A. (2000), Blind separation of nonstationary and temporally correlated
             sources from noisy mixtures, Proceedings of the 2000
            \emph{IEEE Signal Processing Society Workshop Neural Networks for Signal Processing X}, 1, 405--414.}
            
\cite{Nordhausen K. (2014), On robustifying some second order blind source separation methods for nonstationary time series, 
            \emph{Statistical Papers}, 55, 141--156.}

\cite{Miettinen, J., Nordhausen, K. and Taskinen, S. (2017), Blind Source Separation Based on Joint Diagonalization in R: The Packages JADE and BSSasymp, \emph{Journal of Statistical Software}, \bold{76}, 1--31, <doi:10.18637/jss.v076.i02>.}

}

\author{Klaus Nordhausen}


\seealso{
\code{\link{ts}, \link{NSS.JD}, \link{NSS.TD.JD}}
}


\examples{
n <- 1000
s1 <- rnorm(n)
s2 <- 2*sin(pi/200*1:n)* rnorm(n)
s3 <- c(rnorm(n/2), rnorm(100,0,2), rnorm(n/2-100,0,1.5))
S <- cbind(s1,s2,s3)
plot.ts(S)
A<-matrix(rnorm(9),3,3)
X<- S\%*\%t(A)

NSS1 <- NSS.SD(X)
NSS1
MD(coef(NSS1),A)
plot(NSS1)
cor(NSS1$S,S)

NSS1b <- NSS.SD(X, n.cut=400)
MD(coef(NSS1b),A)

NSS1c <- NSS.SD(X, n.cut=c(1,600,1000))
MD(coef(NSS1c),A)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ ts }
