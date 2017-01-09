\name{k_JADE}
\alias{k_JADE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast Equivariant k-JADE Algorithm for ICA 
}
\description{
This algorithm generalizes the \code{\link{JADE}} algorithm, as it provides \code{\link{JADE}} when k is set to the number of dimensions. Otherwise k can be considered as a way to reduce the number of cumulant matrices to be jointly diagonalized. Hence small values of k speed up the method considerably in high-dimensional cases. In general, k can be considered as maximum number of underlying identical sources.  
 
The function uses \code{\link{FOBI}} as an initial estimate and \code{\link{frjd}} for the joint diagonalization.

}
\usage{
k_JADE(X, k = 1, eps = 1e-06, maxiter = 100, na.action = na.fail)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Numeric data matrix or dataframe.}
  \item{k}{integer value between 1 and the number of columns of X. Default is 1.}
  \item{eps}{Convergence tolerance.}
  \item{maxiter}{Maximum number of iterations.}
  \item{na.action}{A function which indicates what should happen when the data
          contain 'NA's.  Default is to fail.}
}
\details{The order of the estimated components is fixed so that their fourth moments are in the decreasing order.
}
\value{
 A list with class 'bss' containing the following components:
  \item{A}{The estimated mixing matrix.}
  \item{W }{The estimated unmixing matrix.}
  \item{S }{Matrix with the estimated independent components.}
  \item{Xmu }{The location of the original data.}
}

\references{
\cite{Miettinen, J., Nordhausen, K., Oja, H. and Taskinen, S. (2013), Fast Equivariant JADE, 
In the Proceedings of \emph{38th IEEE International Conference on Acoustics, Speech, and Signal Processing (ICASSP 2013)}, 6153--6157.}

\cite{Miettinen, J., Nordhausen, K. and Taskinen, S. (2017), Blind Source Separation Based on Joint Diagonalization in R: The Packages JADE and BSSasymp, \emph{Journal of Statistical Software}, \bold{76}, 1--31, <doi:10.18637/jss.v076.i02>.}

}
\author{
Jari Miettinen
}
\note{
The function uses \code{\link{FOBI}} as initial estimate and \code{\link{frjd}} for the joint diagonalization.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{JADE}}, \code{\link{FOBI}}, \code{\link{frjd}}
}
\examples{
# 3 source and 3 signals

S <- cbind(rt(1000, 4), rnorm(1000), runif(1000))
A <- matrix(rnorm(9), ncol = 3)
X <- S \%*\% t(A)
res_k1<-k_JADE(X,1)
res_k1$A
res_k1$W
res_k1$S[1:10,]

MD(coef(res_k1),A) 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }

