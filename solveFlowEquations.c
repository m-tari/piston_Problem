#include "main_fsi.h"
#include "preproc.h"


void solveFlowEquations(int ndiag, double **A, double *B)
{
    int *dist;
	dist = (int *)malloc(2*ndiag*sizeof(int));
	if(dist==NULL)
    {
        printf("can not allocate dist");
    }

    double tolBG    = pow(10,-15);
	int iterMaxbg   = 500;          // maximum iteration for the convergence of residual in smbicgstab
	int info;                       // output of bicgstab
	char matdescra[6],predescra[6];

	matdescra[0]='G'; matdescra[1]='L'; matdescra[2]='N'; matdescra[3]='C'; 
	predescra[0]='D'; predescra[1]='L'; predescra[2]='N'; predescra[3]='C';
	dist[0]=-5; dist[1]=-4; dist[2]=-3; dist[3]=-2; dist[4]=-1; dist[5]=0; dist[6]=1;// vector with the different diagonal number using in the global matrix
	dist[7]=2;  dist[8]=3;  dist[9]=4;  dist[10]=5;

    info = bicgstab(&iterMaxbg,&tolBG,3*bagF->N,matdescra,&(A[0][0]),3*bagF->N,ndiag,dist,predescra,B,bagF->x); //inversion of the matrix using bicgstab algorithm 
	if (info<0 || tolBG!=tolBG)
	{
  		printf("BiCGStab could not solve the implicit system! info=%d iterMaxbg=%d tolBG=%le\n",info,iterMaxbg,tolBG); 
	}

}