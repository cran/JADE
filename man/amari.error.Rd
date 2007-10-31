\name{amari.error}
\alias{amari.error}

\title{Amari Error}
\description{
 Computes the Amari Error to evaluate the performance of an ICA algorithm.
}
\usage{
amari.error(W.hat, A, standardize = F)
}

\arguments{
  \item{W.hat}{The estimated square unmixing matrix W.}
  \item{A}{The true square mixing matrix A.}
  \item{standardize}{Logical value if A and W.hat need to be standardized. Default is FALSE.}
}
\details{
  The Amari Error can be used in simulation studies to evaluate the performance of an
  ICA algorithm. The Amari error is permutation invariant but not scale invariant. Therefore if different
  algorithms should be compared the matrices should be scaled in the same way.
  If \code{standardize} is TRUE, this will be done by the function by standardizing 'W.hat' and the inverse of 'A'
  in such a way, that every row has length 1, the largest absolute value of the row has a positive sign
  and the rows are ordered decreasingly according to their largest values.
  
  Note that this function assumes the ICA model is \eqn{X = S A'}, as is assumed by \code{\link{JADE}} and \code{ics}. However \code{fastICA} and 
  \code{PearsonICA} assume \eqn{X = S A}. Therefore matrices from those functions have to be transposed first.
  
  The Amari Error is scaled in such a way, that it takes a value between 0 and 1. And 0 corresponds to an optimal separation.  
}
\value{
  The value of the Amari Error.
}
\references{
   \cite{Amari, S., Cichocki, A. and Yang, H.H. (1996), A new learning algorithm for blind signal separation, \emph{Advances in Neural Information Processing Systems}, \bold{8},  757--763. }  
}
\author{Klaus Nordhausen, \email{klaus.nordhausen@uta.fi}}

\seealso{\code{\link{ComonGAP}}, \code{\link{SIR}}}
\examples{
S <- cbind(rt(1000, 4), rnorm(1000), runif(1000))
A <- matrix(rnorm(9), ncol = 3)
X <- S \%*\% t(A)

W.hat <- JADE(X, 3)$W
amari.error(W.hat, A)
amari.error(W.hat, A, TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
