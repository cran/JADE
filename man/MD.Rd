\name{MD}
\alias{MD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Minimum Distance index MD}

\description{Computes the Minimum Distance index MD to evaluate the performance of an ICA algorithm.}

\usage{
MD(W.hat, A)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{W.hat}{The estimated square unmixing matrix W.}
  \item{A}{The true square mixing matrix A.}
}
\details{
 \deqn{MD(\hat{W},A)=\frac{1}{\sqrt{p-1}} \inf_{P D}{||PD \hat{W} A-I||,}}{MD(W.hat,A) = 1/sqrt(p-1) inf_(P D) ||P D W.hat A - I||,}
 where \eqn{P}{P} is a permutation matrix and \eqn{D}{D} a diagonal matrix with nonzero diagonal entries.
 
 The step that minimizes the index of the set over all permutation matrix can be expressed as a linear sum assignment problem (LSAP)
 for which we use as solver the Hungarian method implemented as \code{solve_LASP} in the \pkg{coin} package.
 
 Note that this function assumes the ICA model is \eqn{X = S A'}, as is assumed by \code{\link{JADE}} and \code{ics}. However \code{fastICA} and 
 \code{PearsonICA} assume \eqn{X = S A}. Therefore matrices from those functions have to be transposed first.
  
  The MD index is scaled in such a way, that it takes a value between 0 and 1. And 0 corresponds to an optimal separation.  
}
\value{The value of the MD index.}

\references{
\cite{Ilmonen, P., Nordhausen, K., Oja, H. and Ollila, E. (2010): A New Performance Index for ICA: Properties, Computation and Asymptotic Analysis. 
In Vigneron, V., Zarzoso, V., Moreau, E., Gribonval, R. and Vincent, E. (editors) \emph{Latent Variable Analysis and Signal Separation}, 229--236, Springer.}
}
\author{Klaus Nordhausen, \email{klaus.nordhausen@uta.fi}}

\seealso{\code{\link{ComonGAP}}, \code{\link{SIR}}, \code{\link{amari.error}}, \code{solve_LASP}}
\examples{
S <- cbind(rt(1000, 4), rnorm(1000), runif(1000))
A <- matrix(rnorm(9), ncol = 3)
X <- S \%*\% t(A)

W.hat <- JADE(X, 3)$W
MD(W.hat, A)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }

