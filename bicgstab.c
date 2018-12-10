/* BiCGSTAB.f -- translated by f2c (version of 20 August 1993  13:15:44).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)

	
*/

#include "preproc.h"


/*  -- Iterative template routine --
*     Univ. of Tennessee and Oak Ridge National Laboratory
*     October 1, 1993
*     Details of this algorithm are described in "Templates for the
*     Solution of Linear Systems: Building Blocks for Iterative
*     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra,
*     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
*     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*	
*
*  Modification
*  ============
*  Polytechnique Montreal
*  August 8, 2014
*  The algorithm was modified to use Intel MKL functions in C. The f2c library 
*  is no more used. The A matrix is passed in the function call as a diagonal
*  sparse matrix. To use a general matrix, uncomment the cblas_dgemv lines 
*  and comment the mkl_ddiamv lines. The same must be done with the cblas_dtbsv
*  lines and mkl_ddiasv lines. The function prototype should also be updated
*
*	
*  Purpose
*  =======
*
*  BICGSTAB solves the linear system A*x = b using the
*  BiConjugate Gradient Stabilized iterative method with
*  preconditioning.
*
*  Convergence test: ( norm( b - A*x ) / norm( b ) ) < TOL.
*  For other measures, see the above reference.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER.
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
*
*  B       (input) DOUBLE PRECISION array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (input/output) DOUBLE PRECISION array, dimension N.
*          On input, the initial guess. This is commonly set to
*          the zero vector.
*          On exit, if INFO = 0, the iterated approximate solution.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LDW,7)  !!! DEPRECATED: Allocated internaly !!!
*          Workspace for residual, direction vector, etc.
*          Note that vectors R and S shared the same workspace.
*
*  LDW     (input) INTEGER   !!! DEPRECATED: WORK allocated internaly !!!
*          The leading dimension of the array WORK. LDW >= max(1,N).
*
*  ITER    (input/output) INTEGER
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  RESID   (input/output) DOUBLE PRECISION
*          On input, the allowable convergence measure for
*          norm( b - A*x ) / norm( b ).
*          On output, the final value of this measure.
*
*  MATVEC  (external subroutine)  !!! DEPRECATED: MKL FCT USED !!!
*          The user must provide a subroutine to perform the
*          matrix-vector product
*
*               y := alpha*A*x + beta*y,
*
*          where alpha and beta are scalars, x and y are vectors,
*          and A is a matrix. Vector x must remain unchanged.
*          The solution is over-written on vector y.
*
*          The call is:
*
*             CALL MATVEC( ALPHA, X, BETA, Y )
*
*          The matrix is passed into the routine in a common block.
*
*  PSOLVE  (external subroutine)  !!! DEPRECATED: MKL FCT USED !!!
*          The user must provide a subroutine to perform the
*          preconditioner solve routine for the linear system
*
*               M*x = b,
*
*          where x and b are vectors, and M a matrix. Vector b must
*          remain unchanged.
*          The solution is over-written on vector b.
*
*          The call is:
*
*             CALL PSOLVE( X, B )
*
*          The preconditioner is passed into the routine in a common block.
*
*  INFO    (output) INTEGER
*
*          =  0: Successful exit. Iterated approximate solution returned.
*
*          >  0: Convergence to tolerance not achieved. This will be
*                set to the number of iterations performed.
*
*          <  0: Illegal input parameter, or breakdown occurred
*                during iteration.
*
*                Illegal parameter:
*
*                   -1: matrix dimension N < 0
*                   -2: LDW < N      !!! DEPRECATED !!!
*                   -3: Maximum number of iterations ITER <= 0.
*                   -4: lval < N
*
*                BREAKDOWN: If parameters RHO or OMEGA become smaller
*                   than some tolerance, the program will terminate.
*                   Here we check against tolerance BREAKTOL.
*
*                  -10: RHO < BREAKTOL: RHO and RTLD have become
*                                       orthogonal.
*                  -11: OMEGA < BREAKTOL: S and T have become
*                                         orthogonal relative to T'*T.
*
*                  BREAKTOL is set in function GETBREAK.
*
*  BLAS CALLS: DAXPY, DCOPY, DDOT, DNRM2, DSCAL
*  ==============================================================
*/
 
int bicgstab(int *iter,_FLOAT *resid,int n,char *matdescra,_FLOAT *A,int lval,int ndiag,int *dist,char *predescra,_FLOAT *b,_FLOAT *x)
{
   /* System generated locals */
   int work_dim1;
   _FLOAT d__1;
   char transa,up,lo;

   /* Output variable */
   int info;
    
   /* Table of constant values */
   int inc=1;
   _FLOAT c_b5 = -1.;
   _FLOAT c_b6 = 1.;
   _FLOAT c_b25 = 0.;

   /* Local variables */
   double eps = 1.0E-12;
   int p,r,s,t,v,rtld,phat, shat,tmp;
   int maxit; 
   int job[7];
   int *IA,*JA,*ipar;
   _FLOAT alpha,beta,omega,rho,rho1;
   _FLOAT rhotol,omegatol, tol,bnrm2;
   _FLOAT *work,*temp,*P,*dpar;
   
   /* Executable Statements */
   work_dim1=n;
   info = 0;
   transa='N'; up='U'; lo='L';
   
/*     Test the input parameters. */
   if (n < 0)
   {
      info = -1;
   }
   else if (*iter <= 0) 
   {
      info = -3;
   }   
   else if (lval < n)
   {
      info = -4;
   }
   if (info != 0) 
   {
      return info;
   }

   maxit = *iter;
   tol = *resid;

/*     Alias workspace columns. */
   r = 0;
   rtld = 1;
   p = 2;
   v = 3;
   t = 4;
   phat = 5;
   shat = 6;
   s = 0;
   tmp=7;
   
 /*     Workspace allocation. */  
   work  = allocate_1d_array(8*work_dim1,"work");
   
 /*     Preconditionner allocation   */
   temp= allocate_1d_array(lval*ndiag,"temp");
   P  = allocate_1d_array(lval*ndiag,"P");
   JA = (int *) malloc((lval*ndiag)*sizeof(int));
   IA = (int *) malloc((lval+1)*sizeof(int));
   ipar=(int *) malloc(128*sizeof(int));
   dpar=allocate_1d_array(128,"dpar");
   
   job[0]=1; job[1]=1; job[2]=0; job[3]=0; job[4]=0; job[5]=1;
   ipar[1]=6; ipar[5]=1; ipar[30]=1; dpar[30]= eps; dpar[31]=1.0E-12;
   
   mkl_dcsrdia(job,&lval,temp,JA,IA,A,&lval,dist,&ndiag,temp,JA,IA,&info);
          
   dcsrilu0(&lval, temp, IA, JA, P, ipar, dpar, &info);

   free(temp);

   if (info!=0)
   {
      printf("Preconditionner ILU0 was not evaluated correctly: %d\n",info);
      exit(0);   
   }

/*     Set parameter tolerances. */
   rhotol = eps;
   omegatol = eps;

/*     Set initial residual. */
   cblas_dcopy(n, b, inc, &(work[r * work_dim1]), inc);

   if (cblas_dnrm2(n, x,inc)!= 0.) 
   {
    
      //cblas_dgemv(CblasRowMajor,CblasNoTrans, n, n, c_b5, A, n, x, inc, c_b6,&(work[r * work_dim1]), inc);
	   
	   mkl_ddiamv(&transa, &n, &n, &c_b5, matdescra, A, &lval, dist, &ndiag, x, &c_b6, &(work[r * work_dim1]));
	    
      if (cblas_dnrm2(n, &(work[r * work_dim1]),inc) <= tol) 
      {
         goto L30;
	   }
   }
   cblas_dcopy(n,&(work[r * work_dim1]), inc, &(work[rtld * work_dim1]), inc);

   bnrm2 = cblas_dnrm2(n, b,inc);//_SQRT(cblas_ddot(n, b, inc, b, inc));//

   if (bnrm2 == 0.) 
   {
	   bnrm2 = 1.;
   }

   *iter = 0;

L10:

/*     Perform BiConjugate Gradient Stabilized iteration. */
   ++(*iter);
   rho = cblas_ddot(n, &(work[rtld * work_dim1]), inc, &(work[r * work_dim1]), inc);

   if (_FABS(rho) < rhotol) 
   {
      goto L25;
   }

/*        Compute vector P. */
   if (*iter > 1) 
   {
      beta = rho / rho1 * (alpha / omega);
      d__1 = -omega;
      cblas_daxpy(n, d__1, &(work[v * work_dim1]), inc, &(work[p * work_dim1]), inc);
      cblas_dscal (n, beta, &(work[p * work_dim1]), inc);
      cblas_daxpy(n, c_b6, &(work[r * work_dim1]), inc, &(work[p * work_dim1]), inc);
   } 
   else 
   {
      cblas_dcopy(n,&(work[r * work_dim1]), inc, &(work[p * work_dim1]), inc);
   }  

/*        Compute direction adjusting vector PHAT and scalar ALPHA. */
  /* if (predescra[0]=='3')
   {
      cblas_dcopy(n,&(work[p * work_dim1]), inc, &(work[phat * work_dim1]), inc);
      cblas_dcopy(3*n,&(A[work_dim1]), inc, &(M[0]), inc);
      
      info=LAPACKE_dgtsv( LAPACK_COL_MAJOR, n, 1, &(M[0]), &(M[work_dim1]), &(M[2*work_dim1]), &(work[phat * work_dim1]), n );
    
      //tridiagonal(0,n,&(M[0]),&(M[work_dim1]),&(M[2*work_dim1]),&(work[phat * work_dim1]));
   }
   else
   {
      mkl_ddiasv(&transa, &n, &c_b6, predescra, A, &lval, dist, &ndiag, &(work[p * work_dim1]), &(work[phat * work_dim1]));
   }*/
   
   mkl_dcsrtrsv(&lo,&transa,&up,&lval,P,IA,JA,&(work[p * work_dim1]),&(work[tmp * work_dim1]));
   mkl_dcsrtrsv(&up,&transa,&transa,&lval,P,IA,JA,&(work[tmp * work_dim1]),&(work[phat * work_dim1]));

   mkl_ddiamv(&transa, &n, &n, &c_b6, matdescra, A, &lval, dist, &ndiag, &(work[phat * work_dim1]), &c_b25, &(work[v * work_dim1]));

   alpha = rho / cblas_ddot(n, &(work[rtld * work_dim1]), inc, &(work[v * work_dim1]), inc);

/*        Early check for tolerance. */
   d__1 = -alpha;
   cblas_daxpy(n, d__1, &(work[v * work_dim1]), inc, &(work[r * work_dim1]), inc);
   cblas_dcopy(n,&(work[r * work_dim1]), inc, &(work[s * work_dim1]), inc);

   if (cblas_dnrm2(n, &(work[s * work_dim1]),inc) <= tol) 
   {
      cblas_daxpy(n, alpha, &(work[phat * work_dim1]), inc,  x, inc);
      *resid = cblas_dnrm2(n, &(work[s * work_dim1]),inc) / bnrm2;
      goto L30;
   } 
   else 
   {
/*           Compute stabilizer vector SHAT and scalar OMEGA. */
    /*  if (predescra[0]=='3')
      {
         cblas_dcopy(n,&(work[s * work_dim1]), inc, &(work[shat * work_dim1]), inc);
         cblas_dcopy(3*n,&(A[work_dim1]), inc, &(M[0]), inc);
         
         info=LAPACKE_dgtsv( LAPACK_COL_MAJOR, n, 1, &(M[0]), &(M[work_dim1]), &(M[2*work_dim1]), &(work[shat * work_dim1]), n );
         
         //tridiagonal(1,n,&(M[0]),&(M[work_dim1]),&(M[2*work_dim1]),&(work[shat * work_dim1]));
      }
      else
      {
         mkl_ddiasv(&transa, &n, &c_b6, predescra, A, &lval, dist, &ndiag, &(work[s * work_dim1]), &(work[shat * work_dim1]));
      }*/
      
     // cblas_dcopy(n,&(work[s * work_dim1]), inc, &(work[shat * work_dim1]), inc);  // Added line to copy work[p...] to work[phat...]
     // cblas_dtbsv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit, n, 0, M, 1, &(work[shat * work_dim1]), inc);
      //cblas_dgemv(CblasRowMajor,CblasNoTrans, n, n, c_b6, A, n, &(work[shat * work_dim1]), inc, c_b25,&(work[t * work_dim1]), inc);
      
      mkl_dcsrtrsv(&lo,&transa,&up,&lval,P,IA,JA,&(work[s * work_dim1]),&(work[tmp * work_dim1]));
      mkl_dcsrtrsv(&up,&transa,&transa,&lval,P,IA,JA,&(work[tmp * work_dim1]),&(work[shat * work_dim1]));
      
      mkl_ddiamv(&transa, &n, &n, &c_b6, matdescra, A, &lval, dist, &ndiag, &(work[shat * work_dim1]), &c_b25, &(work[t * work_dim1]));

      omega = cblas_ddot(n, &(work[t * work_dim1]), inc, &(work[s * work_dim1]), inc)/cblas_ddot(n, &(work[t * work_dim1]), inc, &(work[t * work_dim1]), inc);

/*           Compute new solution approximation vector X. */
      cblas_daxpy(n, alpha, &(work[phat * work_dim1]), inc,  x, inc);

      cblas_daxpy(n, omega, &(work[shat * work_dim1]), inc,  x, inc);

/*           Compute residual R, check for tolerance. */
      d__1 = -omega;
      cblas_daxpy(n, d__1, &(work[t * work_dim1]), inc, &(work[r * work_dim1]), inc);

      *resid = cblas_dnrm2(n, &(work[r * work_dim1]),inc) / bnrm2;

      if (*resid <= tol) 
      {
         goto L30;
      }
      if (*iter == maxit) 
      {
         goto L20;
      }
   }

   if (_FABS(omega) < omegatol) 
   {
      goto L25;
   } 
   else 
   {
      rho1 = rho;
      goto L10;
   }

L20:

/*     Iteration fails. */
   info = 1;
   return info;

L25:

/*     Set breakdown flag. */

   if (_FABS(rho) < rhotol) 
   {
      info = -10;
   } 
   else if (_FABS(omega) < omegatol) 
   {
      info = -11;
   }
   return info;

L30:

/*     Iteration successful; return. */
   free(work); free(P); free(dpar); free(IA); free(JA); free(ipar);

   return info;

/*     End of BICGSTAB */

}
