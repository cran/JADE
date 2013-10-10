\name{cjd}
\alias{cjd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Joint Diagonalization of Complex Matrices}
\description{
This is an \pkg{R} version of Cardoso's joint_diag matlab function for joint diagonalization of k complex-valued square matrices.
}
\usage{
cjd(X, eps = 1e-06, maxiter = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{A matrix of k stacked pxp complex matrices with dimension c(kp,p) or an array with dimension c(p,p,k).}
  \item{eps}{ Convergence tolerance.}
  \item{maxiter}{ Maximum number of iterations.}
}

\value{
    \item{V }{An orthogonal matrix.}
    \item{D }{A stacked matrix with the diagonal matrices or an array with the diagonal matrices. The form of the output
            depends on the form of the input.}
}
\references{
   \cite{Cardoso, J.-F. and Souloumiac, A., (1996), Jacobi angles for simultaneous diagonalization, \emph{SIAM J. Mat. Anal. Appl.}, \bold{17},  161--164.}  
}
\author{Jean-Francois Cardoso. Ported to \pkg{R} by Klaus Nordhausen.}

\seealso{
\code{\link{rjd}}, \code{\link{rjd.fortran}}
}
\examples{
D1 <- diag(complex(real=runif(3,0,2), imaginary=runif(3)))
D2 <- diag(complex(real=runif(3,0,2), imaginary=runif(3)))
D3 <- diag(complex(real=runif(3,0,2), imaginary=runif(3)))
D4 <- diag(complex(real=runif(3,0,2), imaginary=runif(3)))

Z <- matrix(runif(9), ncol = 3)
V <- eigen(Z \%*\% t(Z))$vectors

M1 <- t(V)\%*\%D1\%*\%V
M2 <- t(V)\%*\%D2\%*\%V
M3 <- t(V)\%*\%D3\%*\%V
M4 <- t(V)\%*\%D4\%*\%V
MS <- rbind(M1,M2,M3,M4)
Ms <- array(0,dim=c(3,3,4))
Ms[,,1]<-M1
Ms[,,3]<-M3
Ms[,,2]<-M2
Ms[,,4]<-M4
res.array <- cjd(Ms)
res.mat <- cjd(MS)
Re(res.array$V)
V
round(V\%*\%Re(res.array$V),2)
round(V\%*\%Re(res.mat$V),2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ array }

