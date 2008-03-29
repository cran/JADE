\name{rjd}
\alias{rjd}

\title{Joint Diagonalization of Real Matrices  }
\description{
  This is an \pkg{R} version of Cardoso's rjd matlab function for joint diagonalization of k real-valued square matrices.
  
}
\usage{
rjd(X, eps = 1e-06, maxiter = 100, na.action = na.fail)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{A matrix of k stacked pxp matrices with dimension c(kp,p) or an array with dimension c(p,p,k).}
  \item{eps}{ Convergence tolerance.}
  \item{maxiter}{ Maximum number of iterations.}
  \item{na.action}{A function which indicates what should happen when the data
          contain 'NA's.  Default is to fail.}
}
\details{
  Denote the square matrices as \eqn{A_i}{A_i}, \eqn{i=1,\ldots,k}{i=1,...,k}. This algorithm searches then an orthogonal matrix V
  so that  \eqn{D_i=V'A_iV}{D_i=V'A_iV} is diagonal for all \eqn{i}. If the \eqn{A_i}{A_i} commute then there is an exact solution. If not, the function
  will perform an approximate joint diagonalization by trying to make the \eqn{D_i}{D_i} as diagonal as possible.
  
  Cardoso points out that notion of approximate joint diagonalization
  is ad hoc and very small values of \code{eps} make in that case not much sense since the diagonality 
  criterion is ad hoc itself.
}
\value{
  A list with the components
  \item{V }{An orthogonal matrix.}
  \item{D }{A stacked matrix with the diagonal matrices or an array with the diagonal matrices. The form of the output
            depends on the form of the input.}
}


\references{
   \cite{Cardoso, J.-F. and Souloumiac, A., (1996), Jacobi angles for simultaneous diagonalization, \emph{SIAM J. Mat. Anal. Appl.}, \bold{17},  161--164.}  
}
\author{Jean-Francois Cardoso. Ported to \pkg{R} by Klaus Nordhausen, \email{klaus.nordhausen@uta.fi}}

\examples{
Z <- matrix(runif(9), ncol = 3)
U <- eigen(Z \%*\% t(Z))$vectors
D1 <- diag(runif(3))
D2 <- diag(runif(3))
D3 <- diag(runif(3))
D4 <- diag(runif(3))

X.matrix <- rbind(t(U) \%*\% D1 \%*\% U, t(U) \%*\% D2 \%*\% U,
                  t(U) \%*\% D3 \%*\% U, t(U) \%*\% D4 \%*\% U)
res.matrix <- rjd(X.matrix)
res.matrix$V
round(U \%*\% res.matrix$V, 4) # should be a signed permutation matrix if V is correct.

round(res.matrix$D, 4)

X.array <- aperm(array(t(X.matrix), dim = c(3,3,4)), c(2,1,3))

res.array <- rjd(X.array)
round(res.array$D, 4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ array }

