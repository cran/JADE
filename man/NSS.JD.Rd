\name{NSS.JD}
\alias{NSS.JD}
\alias{NSS.JD.default}
\alias{NSS.JD.ts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
NSS.JD Method for Nonstationary Blind Source Separation
}
\description{
The NSS.JD method for nonstationary blind source separation. The method first whitens the complete data and then divides it into K time intervals.
Then \code{\link{rjd.fortran}} is used to jointly diagonalize the covariance matrices computed for the individual time intervals to find the sources.
}
\usage{
NSS.JD(X, ...)

\method{NSS.JD}{default}(X,  K=12, Tau=0, n.cuts=NULL, eps = 1e-06, maxiter = 100, ...)
\method{NSS.JD}{ts}(X, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a numeric matrix or a multivariate time series object of class \code{\link{ts}}. Missing values are not allowed.}
  \item{K}{number of intervals to be used.}
  \item{Tau}{By default 0 which means covariance are computed of each time interval, if Tau is an integer > 0 then rather autocovariance matrices at lag Tau are used for the joint diagonaliation.}
  \item{n.cuts}{if NULL, then the time series is divided into K equally long intervals. To specify intervals n.cuts should be given in the form c(1,n.cut.1,...,n.cut.k, nrow(X)) to specify where to split the time series.}
  \item{eps}{maximum number of iterations for \code{\link{rjd.fortran}}.}
  \item{maxiter}{convergence tolerance for \code{\link{rjd.fortran}}.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{
The model assumes that the mean of the p-variate time series is constant but the variances change over time.
}
\value{
 A list with class 'bss' containing the following components:
  \item{W}{estimated unmixing matrix.}
  \item{k}{the lag used for the autocovariance matrix.}
  \item{n.cut}{specifying the intervals where data is split}
  \item{K}{the number of intervals used}
  \item{S}{estimated sources as time series objected standardized to have mean 0 and that the variance of the sources are 1.}
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
\code{\link{ts}, \link{NSS.SD}, \link{NSS.TD.JD}}
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

NSS2 <- NSS.JD(X)
NSS2
MD(coef(NSS2),A)
plot(NSS2)
cor(NSS2$S,S)

NSS2b <- NSS.JD(X, Tau=1)
MD(coef(NSS2b),A)

NSS2c <- NSS.JD(X, n.cuts=c(1,300,500,600,1000))
MD(coef(NSS2c),A)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ ts }
