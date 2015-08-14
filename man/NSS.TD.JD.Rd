\name{NSS.TD.JD}
\alias{NSS.TD.JD}
\alias{NSS.TD.JD.default}
\alias{NSS.TD.JD.ts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
NSS.TD.JD Method for Nonstationary Blind Source Separation
}
\description{
The NSS.TD.JD method for nonstationary blind source separation. The method first whitens the complete data and then divides it into K time intervals.
It is then assumed that within each interval the time series is approximately second order stationary and within each interval L autocovariance are computed.
The underlying sources are then found by jointly diagonalizing the K*L autocovariance matrices using  \code{\link{frjd}}.
}
\usage{
NSS.TD.JD(X, ...)

\method{NSS.TD.JD}{default}(X, K=12, Tau=0:11, n.cuts=NULL, eps = 1e-06, maxiter = 100, ...)
\method{NSS.TD.JD}{ts}(X, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a numeric matrix or a multivariate time series object of class \code{\link{ts}}. Missing values are not allowed.}
  \item{K}{number of intervals to be used.}
  \item{Tau}{Lags for the autovariance matrices to be computed within each interval.}
  \item{n.cuts}{if NULL, then the time series is divided into K equally long intervals. To specify intervals n.cuts should be given in the form c(1,n.cut.1,...,n.cut.k, nrow(X)) to specify where to split the time series.}
  \item{eps}{maximum number of iterations for \code{\link{frjd}}.}
  \item{maxiter}{convergence tolerance for \code{\link{frjd}}.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{
The model assumes that the mean of the p-variate time series is constant but the variances change over time.
}
\value{
 A list with class 'bss' containing the following components:
  \item{W}{estimated unmixing matrix.}
  \item{k}{the lags used for the autocovariance matrix used in each interval.}
  \item{n.cut}{specifying the intervals where data is split}
  \item{K}{the number of intervals used}
  \item{S}{estimated sources as time series objected standardized to have mean 0 and that the sources 1.}
}

\references{
\cite{Choi S. and Cichocki A. (2000), Blind separation of nonstationary sources in noisy mixtures, 
            \emph{Electronics Letters}, 36, 848--849.}

\cite{Choi S. and Cichocki A. (2000), Blind separation of nonstationary and temporally correlated
             sources from noisy mixtures, Proceedings of the 2000
            \emph{IEEE Signal Processing Society Workshop Neural Networks for Signal Processing X}, 1, 405--414.}
            
\cite{Nordhausen K. (2013), On robustifying some second order blind source separation methods for nonstationary time series, 
            \emph{to appear in Statistical Papers}, ??, ???--???.}

}

\author{Klaus Nordhausen}


\seealso{
\code{\link{ts}, \link{NSS.JD}, \link{NSS.JD}, \link{SOBI}}
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

NSS3 <- NSS.TD.JD(X)
NSS3
MD(coef(NSS3),A)
plot(NSS3)
cor(NSS3$S,S)

NSS3b <- NSS.TD.JD(X, Tau=c(0,3,7,12), K=6)
MD(coef(NSS3b),A)

NSS3c <- NSS.TD.JD(X, n.cuts=c(1,300,500,600,1000))
MD(coef(NSS3c),A)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ ts }
