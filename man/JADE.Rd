\name{JADE}
\alias{JADE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ JADE algorithm for ICA }
\description{
  This is an \pkg{R} version of Cardoso's JADE ICA algorithm for real data ported from matlab. The ported version is 1.5, some minor
  changes compared to the matlab function are explained in the details section. The matlab code can be found for example on the 
  ICA central homepage.
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
  Worth mentioning is also that the estimated independent components \eqn{S} are scaled to unit variance and are ordered in such a way, that the most energetic component comes first.
  The signs of the unmixing matrix \eqn{W} are fixed such, that the first column of \eqn{W} has positive elements.
  
  For further details see also the documention of the original matlab code ("MatlabjadeR.m") on the ICA central homepage (\url{http://www.tsi.enst.fr/icacentral/}). 
}
\value{
  A list with the components
  \item{A}{The estimated mixing matrix.}
  \item{W }{The estimated unmixing matrix.}
  \item{S }{Dataframe with the estimated independent components.}
  \item{Xmu }{The location of the original data.}
}
\references{
   \cite{Cardoso, J.-F. and Souloumiac, A., (1993), Blind beamforming for non Gaussian signals, \emph{IEE Proceedings-F}, \bold{140},  362--370. 
   \url{ftp://sig.enst.fr/pub/jfc/Papers/iee.ps.gz}.}  
}
\author{Jean-Francois Cardoso. Ported to \pkg{R} by Klaus Nordhausen, \email{klaus.nordhausen@uta.fi}}


\examples{
S <- cbind(rt(1000, 4), rnorm(1000), runif(1000))
A <- matrix(rnorm(9), ncol = 3)
X <- S \%*\% t(A)
res<-JADE(X,3)
res$A
res$W
res$S[1:10,]
(sweep(X,2,res$Xmu) \%*\% t(res$W))[1:10,]
round(res$W \%*\% A,4) 
}
\keyword{ multivariate }
