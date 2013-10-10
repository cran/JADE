\name{AMUSE}
\alias{AMUSE}
\alias{AMUSE.default}
\alias{AMUSE.ts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  AMUSE Method for Blind Source Separation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
AMUSE method for the second order blind source separation problem. The function estimates
the unmixing matrix in a second order stationary source separation model by jointly diagonalizing
the covariance matrix and an autocovariance matrix at lag k.
}
\usage{
AMUSE(x, ...)

\method{AMUSE}{default}(x, k = 1, \dots)
\method{AMUSE}{ts}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric matrix or a multivariate time series object of class \code{\link{ts}}. Missing values are not allowed.}
  \item{k}{integer lag for the autocovariance matrix, must be larger than 0. Default is 1.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{
The lag \code{k} has a huge effect on the performance and it should be chosen so that the eigenvalues of autocovariance matrix are distinct. The function assumes always as many sources as there are time series.
}
\value{
 A list with class 'bss' containing the following components:
  \item{W}{estimated unmixing matrix.}
  \item{EV}{eigenvectors of autocovariance matrix.}
  \item{k}{lag of the autocovariance matrix used.}
  \item{S}{estimated sources as time series objected standardized to have mean 0 and unit variances.}
}

\references{
\cite{Tong, L., Soon, V.C., Huang, Y.F. and Liu, R. (1990), AMUSE: a new blind identification algorithm, in Proceedings of {IEEE} International Symposium on
            \emph{Circuits and Systems 1990}, 1784--1787.}
            
\cite{Miettinen, J., Nordhausen, K., Oja, H. and Taskinen, S. (2012), Statistical properties of a blind source separation estimator for stationary time series, 
            \emph{Statistics & Probability Letters}, 82, 1865--1873.}
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

res1<-AMUSE(X)
res1
coef(res1)
plot(res1) # compare to plot.ts(S)
MD(coef(res1),A)

# input of a time series
X2<- ts(X, start=c(1961, 1), frequency=12)
plot(X2)
res2<-AMUSE(X2, k=2)
plot(res2)
}
\keyword{ multivariate }
\keyword{ ts }
