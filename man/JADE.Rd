\name{JADE}
\alias{JADE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ JADE Algorithm for ICA }
\description{
  This is an \pkg{R} version of Cardoso's JADE ICA algorithm (for real data) ported from matlab. The ported version is 1.5. Some minor changes compared to the matlab function are explained in the details section. The matlab code can be found for example on the ICA central homepage.

The function uses \code{\link{frjd}} for the joint diagonalization.
}
\usage{
JADE(X, n.comp = NULL, eps = 1e-06, maxiter = 100, na.action = na.fail)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Numeric data matrix or dataframe. }
  \item{n.comp}{Number of components to extract.}
  \item{eps}{Convergence tolerance.}
  \item{maxiter}{Maximum number of iterations.}
  \item{na.action}{A function which indicates what should happen when the data
          contain 'NA's.  Default is to fail.}
}


\details{
  Some minor modifications were done when porting the function to \pkg{R}, and they are:
  \describe{
  \item{1}{the model assumed here is \eqn{X=S A' +\mu}{X = S A' + Xmu}. Therefore \eqn{S} and \eqn{X} have one row per observation. Note that this still differs from 
  the model definition in \pkg{R} of \code{FastICA} and \code{PearsonICA} but agrees with \code{ics}.}
  \item{2}{The whitening covariance matrix is divided by n-1 and not n (n = number of observations).}
  \item{3}{The initial value for the joint diagonalisation is always I.}
  \item{4}{The original eps would be \eqn{\frac{1}{100\sqrt{n}}}{1/(100*sqrt(n))}.}
  }
  It is also worth mentioning that the estimated independent components \eqn{S} are scaled to unit variance and are ordered in such a way, that their fourth moments are in the decreasing order.
  The signs of the unmixing matrix \eqn{W} are fixed so that the sum of the elements on each row is positive.
  
  For further details see also the documention of the original matlab code ("MatlabjadeR.m") on the ICA central homepage (\url{http://www.tsi.enst.fr/icacentral/}). 
}
\value{
  A list with class 'bss' containing the following components:
  \item{A}{The estimated mixing matrix.}
  \item{W }{The estimated unmixing matrix.}
  \item{S }{Dataframe with the estimated independent components.}
  \item{Xmu }{The location of the original data.}
}
\references{
   \cite{Cardoso, J.-F. and Souloumiac, A., (1993), Blind beamforming for non Gaussian signals, \emph{IEE Proceedings-F}, \bold{140},  362--370.}  
}
\author{Jean-Francois Cardoso. Ported to \pkg{R} by Klaus Nordhausen}


\examples{
# 3 source and 3 signals

S <- cbind(rt(1000, 4), rnorm(1000), runif(1000))
A <- matrix(rnorm(9), ncol = 3)
X <- S \%*\% t(A)
res<-JADE(X,3)
res$A
res$W
res$S[1:10,]
(sweep(X,2,res$Xmu) \%*\% t(res$W))[1:10,]
round(res$W \%*\% A,4) 

# 2 sources and 3 signals

S2 <- cbind(rt(1000, 4), rnorm(1000)) 
A2 <- matrix(rnorm(6), ncol = 2)
X2 <- S2 \%*\% t(A2)
res2 <-JADE(X2,2)
res2$A
res2$W
res2$S[1:10,]
(sweep(X2,2,res2$Xmu) \%*\% t(res2$W))[1:10,]
SIR(S2,res2$S)
}
\keyword{ multivariate }
