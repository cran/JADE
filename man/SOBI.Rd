\name{SOBI}
\alias{SOBI}
\alias{SOBI.default}
\alias{SOBI.ts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
SOBI Method for Blind Source Separation
}
\description{
The SOBI method for the second order blind source separation problem. The function estimates
the unmixing matrix in a second order stationary source separation model by jointly diagonalizing
the covariance matrix and several autocovariance matrices at different lags.
}
\usage{
SOBI(X, ...)

\method{SOBI}{default}(X, k=12, method="frjd", eps = 1e-06, maxiter = 100, ...)
\method{SOBI}{ts}(X, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a numeric matrix or a multivariate time series object of class \code{\link{ts}}. Missing values are not allowed.}
  \item{k}{if a single integer, then the lags 1:k are used, if an integer vector, then these are used as the lags.}
  \item{method}{method to use for the joint diagonalization, options are \code{\link{djd}}, \code{\link{rjd}} and \code{\link{frjd}}}. 
  \item{eps}{convergence tolerance.}
  \item{maxiter}{maximum number of iterations.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{The order of the estimated components is fixed so that the sums of squared autocovariances are in the decreasing order.
}
\value{
 A list with class 'bss' containing the following components:
  \item{W}{estimated unmixing matrix.}
  \item{k}{lags used.}
  \item{method}{method used for the joint diagonalization.}
  \item{S}{estimated sources as time series objected standardized to have mean 0 and unit variances.}
}
\references{
\cite{Belouchrani, A., Abed-Meriam, K., Cardoso, J.F. and Moulines, R. (1997), A blind source separation technique using second-order statistics, 
            \emph{IEEE Transactions on Signal Processing}, 434--444.}
            
\cite{Miettinen, J., Nordhausen, K., Oja, H. and Taskinen, S. (2014),  Deflation-based Separation of Uncorrelated Stationary Time Series, 
            \emph{Journal of Multivariate Analysis}, 123, 214--227.} 

\cite{Miettinen, J., Illner, K., Nordhausen, K., Oja, H., Taskinen, S. and Theis, F.J. (2016),  Separation of Uncorrelated Stationary Time Series Using Autocovariance Matrices, 
            \emph{Journal of Time Series Analysis}, 37, 337--354.} 

\cite{Miettinen, J., Nordhausen, K. and Taskinen, S. (2017), Blind Source Separation Based on Joint Diagonalization in R: The Packages JADE and BSSasymp, \emph{Journal of Statistical Software}, \bold{76}, 1--31, <doi:10.18637/jss.v076.i02>.}

}
\author{Klaus Nordhausen}

\seealso{
\code{\link{ts}}
}
\examples{
# creating some toy data
A<- matrix(rnorm(9),3,3)
s1 <- arima.sim(list(ar=c(0.3,0.6)),1000)
s2 <- arima.sim(list(ma=c(-0.3,0.3)),1000)
s3 <- arima.sim(list(ar=c(-0.8,0.1)),1000)

S <- cbind(s1,s2,s3)
X <- S \%*\% t(A)

res1<-SOBI(X)
res1
coef(res1)
plot(res1) # compare to plot.ts(S)
MD(coef(res1),A)

# input of a time series
X2<- ts(X, start=c(1961, 1), frequency=12)
plot(X2)
res2<-SOBI(X2, k=c(5,10,1,4,2,9,10))
plot(res2)
}
\keyword{ multivariate }
\keyword{ ts }
