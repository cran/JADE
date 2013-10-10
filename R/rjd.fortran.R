rjd.fortran<-function(X, weight=NULL, maxiter=100, eps=1e-06, na.action = na.fail)
{
  X<-na.action(X)
  dim.X<-dim(X)
    
  if (length(dim.X)==2) type<-"Matrix"
  if (length(dim.X)==3) type<-"Array"
  if ((length(dim.X) %in% c(2,3))==F) stop("'X' must have two or three dimensions")
    
  if (type == "Matrix")
     {
      p<-dim.X[2]
      K<-dim.X[1]/p
      if (floor(K) != ceiling(K)) stop("'X' must be a matrix of k stacked pxp matrices")
      X<-array(t(X),c(p,p,K))
      dim.X<-dim(X)
     }
   
  if (dim.X[1] != dim.X[2]) stop("'X' must be an array with dim of the form c(p,p,K)")
  p<-dim(X)[1]
  K<-dim(X)[3]
  
  if(is.null(weight)) weight<-rep(1,K) 

  B<-diag(p)
  res<-.Fortran("FRJD", A=as.double(as.vector(X)),RN=as.double(weight),B=as.double(B),IP=as.integer(p),K=as.integer(K),IFAULT=as.integer(0),NF=as.integer(maxiter),EPSF=as.double(eps),result=double(1),PACKAGE="JADE") 
  if (res$IFAULT != 0) stop("maxiter reached without convergence")
  b<-res$B
  V<-t(matrix(b,p,p))
  D<-array(res$A,c(p,p,K))
  if(type == "Matrix"){
   D <- aperm(D, c(1,3,2))
   D <- matrix(D, ncol=p)
  }

  list(V=V, D=D, iter=res$NF)
}

