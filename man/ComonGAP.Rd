\name{ComonGAP}
\alias{ComonGAP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Comon's Gap}
\description{
  Comon's GAP criterion to evaluate the performance of an ICA algorithm.
}
\usage{
ComonGAP(A, A.hat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{The true square mixing matrix.}
  \item{A.hat}{The estimated square mixing matrix.}
}

\details{
  Comon's GAP criterion is permutation and scale invariant. It can take every positive value and 0 corresponds to an optimal separation.
  If \code{A} is however nearly singular the values of the criterion can be huge.
  
  Note that this function assumes the ICA model is \eqn{X = S A'}, as is assumed by \code{\link{JADE}} and \code{ics}. However \code{fastICA} and 
  \code{PearsonICA} assume \eqn{X = S A}. Therefore matrices from those functions have to be transposed first.
}
\value{
  The value of the Comon's GAP.
}
\references{
   \cite{Comon, P., (1994), Independent Component Analysis, A new concept?, \emph{Signal Processing}, \bold{36},  287--314. }  
}

\author{Klaus Nordhausen}

\seealso{\code{\link{amari.error}}, \code{\link{SIR}}}
\examples{
S <- cbind(rt(1000, 4), rnorm(1000), runif(1000))
A <- matrix(rnorm(9), ncol = 3)
X <- S \%*\% t(A)

A.hat <- JADE(X, 3)$A
ComonGAP(A, A.hat)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }


