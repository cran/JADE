#include <cmath>
//#include <R.h>
//#include <Rdefines.h>
//#include <Rinternals.h>


extern "C" {

using namespace std;

  double **prepmat(double *X, int n, int k)
  {
    int i;
    int j;
    double **Y = new double* [n];
    for (i=0; i<n; i++) Y[i]=new double [k];
    for (i=0;i<n; i++) 
      for (j=0;j<k;j++)
    Y[i][j]=X[j*n+i];
    return Y;
  }

  double **mult(double **X, double **Y, int p)
  {
    int i;
    int j;
    int k;
    double **Z1 = new double* [p];
    for (i=0; i<p; i++) Z1[i] = new double [p];
    double **Z = new double* [p];
    for (i=0; i<p; i++) Z[i] = new double [p];
    for (i=0; i<p; i++){ 
     for (j=0; j<p; j++){
      Z1[i][j] = 0;
      for (k=0; k<p; k++)
       Z1[i][j] += Y[k][i]*X[k][j]; 
     }  
    }
    for (i=0; i<p; i++){ 
     for (j=0; j<(i+1); j++){
      Z[i][j]=0;
      for (k=0; k<p; k++){     
       Z[i][j] += Z1[i][k]*Y[k][j];
      } 
      Z[j][i] = Z[i][j];
     }  
    }
    
   return Z;
  }

  void rjdc(double *X, int *kpmaxit, double *w, double *eps, double *result)
  {
    int K = kpmaxit[0]; 
    int p = kpmaxit[1];
    int maxiter = kpmaxit[2]; 
    int i;
    int j;
    int k;
    int l;
    double iter = 0;
    double alpha;
    double beta; 
    double a1;
    double b1; 
    double a2;
    double b2; 
    double u;
    double c;
    double si;
    double co;
    double theta;     
    bool cc = 1;
    double **A = prepmat(X,K*p,p);  
    double **B = new double* [p];
    for (i=0; i<p; i++) B[i]=new double [p];
    for (i=0; i<(p-1); i++){
     B[i][i] = 1;
     for (j=(i+1); j<p; j++){ 
      B[i][j] = 0;
      B[j][i] = 0;
     }
    }
    B[p-1][p-1] = 1; 
    
    while(cc){ 
     iter++;
     if(iter>maxiter){
      B[0][0] = 2;   
      break;
     } 
     cc = 0;
     for(i=0; i<(p-1); i++){
      for(j=(i+1); j<p; j++){
       alpha = 0;
       beta = 0;
       for(k=0; k<K; k++){
	 u = A[k*p+i][i]-A[k*p+j][j];
         c = A[k*p+i][j];
         alpha += w[k]*(u*u-4*c*c);
         beta += w[k]*4*u*c; 
       }
       theta = atan(beta/(alpha+sqrt(alpha*alpha+beta*beta)))/2;
       si = sin(theta);
       co = cos(theta);
       if(abs(si)>eps[0]) cc = 1;
       for(l=0; l<p; l++){
	 b1 = B[i][l];
 	 b2 = B[j][l];
         B[i][l] = co*b1+si*b2;
         B[j][l] = co*b2-si*b1;
       }
       for(k=0; k<K; k++){
        for(l=0; l<p; l++){
         a1 = A[k*p+i][l];
         a2 = A[k*p+j][l];
	 A[k*p+i][l] = co*a1+si*a2;
	 A[k*p+j][l] = co*a2-si*a1;
        }
        for(l=0; l<p; l++){
         a1 = A[k*p+l][i];
         a2 = A[k*p+l][j];
         A[k*p+l][i] = co*a1+si*a2;
	 A[k*p+l][j] = co*a2-si*a1;
        }
       }
     }
    }
   }

   l=0;
   for(i=0;i<p;i++){   
    for(j=0;j<p;j++){
     result[l] = B[i][j];
     l++;
    }
   }
   result[p*p] = iter;
 
   for(i=0;i<(K*p);i++)   
     delete [] A[i]; 
   delete [] A;
  
   for(i=0;i<p;i++)   
     delete [] B[i]; 
   delete [] B;
 }

 
 void FG(double *X, double *b, int *kpmaxit, double *w, double *eps, double *result)
  {
    int K = kpmaxit[0]; 
    int p = kpmaxit[1];
    int maxiter = kpmaxit[2]; 
    int i;
    int j;
    int k;
    int l;
    int m;
    int n;
    double iter = 0;
    double t;
    double t1;
    double t2; 
    double b1;
    double b2;  
    double u11;
    double u12; 
    double u22; 
    double *T11 = new double [K];
    double *T12 = new double [K]; 
    double *T22 = new double [K]; 
    double s;
    double sold;
    double c;
    double r;
    double *delta1 = new double [K];     
    double *delta2 = new double [K];
    double *coef = new double [K];
    bool cc = 1;
    double **A = prepmat(X,K*p,p);  
    double **F = prepmat(X,K*p,p);  
    double **B = prepmat(b,p,p);
    double **H = prepmat(b,p,p);

     while(cc){ 
     iter++;
     if(iter>maxiter) break;
  
     cc = 0;
     
     for(k=0; k<K; k++){      
      for(i=0; i<p; i++){
       for(j=0; j<p; j++){
         H[i][j] = A[k*p+i][j];
       }
      }
      H = mult(H,B,p);
      for(l=0; l<p; l++){
       for(m=0; m<p; m++){
         F[k*p+l][m] = H[l][m];
       }
      }
     }

     for(l=0; l<p; l++){
      for(m=0; m<p; m++){
        H[l][m] = B[l][m];
      }
     }


     for(i=0; i<(p-1); i++){
      for(j=(i+1); j<p; j++){
       for(k=0; k<K; k++){
         T11[k] = F[k*p+i][i];
         T12[k] = F[k*p+i][j];  
         T22[k] = F[k*p+j][j];     
       }
     
       c = 1;
       s = 0;
       for(n=0; n<5; n++){
        sold = s;
        for(k=0; k<K; k++){
         delta1[k] = c*c*T11[k]+s*s*T22[k]+2*c*s*T12[k]; 
         delta2[k] = c*c*T22[k]+s*s*T11[k]-2*c*s*T12[k]; 
         coef[k] = w[k]*(delta1[k]-delta2[k])/(delta1[k]*delta2[k]);
        }
        u11 = 0;
        u12 = 0;
        u22 = 0;
        for(k=0; k<K; k++){
         u11 += coef[k]*T11[k]; 
         u12 += coef[k]*T12[k]; 
         u22 += coef[k]*T22[k]; 
        }
       

        if(abs(u12)>0.0000000001){ 
         t1 = ((u22-u11)/u12+sqrt(((u22-u11)/u12)*(u22-u11)/u12+4))/2;
         t2 = ((u22-u11)/u12-sqrt(((u22-u11)/u12)*(u22-u11)/u12+4))/2;
         t = t1;
         if(abs(t1) > abs(t2)) t = t2;
         c = 1/sqrt(1+t*t);
         s = t*c;  
        }
        if(abs(s-sold) < 0.0001) break;
       }
      
       r = s/(1+c);
       t1 = s/c;
       t2 = t1*t1;
       for(k=0; k<K; k++){
        F[k*p+i][i] = c*c*(T11[k]+2*t1*T12[k]+t2*T22[k]);
        F[k*p+j][j] = c*c*(t2*T11[k]-2*t1*T12[k]+T22[k]);
        F[k*p+i][j] = c*c*(t1*(T22[k]-T11[k])+(1-t2)*T12[k]);
        F[k*p+j][i] = F[k*p+i][j];
        for(l=0; l<p; l++){
         if((l!=i)&&(l!=j)){
          b1 = F[k*p+l][i]; 
          b2 = F[k*p+l][j];
          F[k*p+l][i] = b1+s*(b2-r*b1);
          F[k*p+l][j] = b2-s*(b1+r*b2);
          F[k*p+i][l] = F[k*p+l][i];
          F[k*p+j][l] = F[k*p+l][j];
         }
        } 
       }  
       for(l=0; l<p; l++){
        b1 = B[l][i]; 
        b2 = B[l][j];
        B[l][i] = b1+s*(b2-r*b1);
        B[l][j] = b2-s*(b1+r*b2);
       } 
      }
     } 
    
     t=0;
     for(i=0; i<p; i++){
      for(j=0; j<p; j++){
       t1 = abs(B[i][j]-H[i][j]);
       if(t1>t) t = t1;   
      }
     }
     if(t>eps[0]) cc = 1;
    }

    l=0;
    for(i=0;i<p;i++){   
     for(j=0;j<p;j++){
      result[l] = B[j][i];
      l++;
     }
    }
    result[p*p] = iter;
   
    for(i=0;i<(K*p);i++)   
     delete [] A[i]; 
    delete [] A;

    for(i=0;i<(K*p);i++)   
     delete [] F[i]; 
    delete [] F;
  
    for(i=0;i<p;i++)   
     delete [] B[i]; 
    delete [] B;
 
    for(i=0;i<p;i++)   
     delete [] H[i]; 
    delete [] H;
    
    delete [] T11;
    delete [] T12;
    delete [] T22;
    delete [] delta1;
    delete [] delta2;
    delete [] coef;
 }

}
