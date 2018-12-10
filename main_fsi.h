#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define _FLOAT double
#define  gam   1.4

typedef struct S_bagS S_bagS;
typedef struct S_bagF S_bagF;

struct S_bagS
{
	double mp;  // piston mass
	double dt;  // structure timestep size
	double kp;  // spring stifness
	double A;   // piston area

	double u;
	double u_dot;
	double u_ddot;

	double u_p;
	double u_dot_p;
	double u_ddot_p;	
};

struct S_bagF
{
	double p0;  // pressure at chamber at rest

	double dt;
	double total_time;	
	double dx;
	double L;

	double q1;
	double q2;
	double q3;

	double *x;
	double *rho;
	double *u;
	double *e;
	double *p;
	double *H;	
	int N;

	double *w1;
	double *w2;
	double *w3;
	double *w1_p;
	double *w2_p;
	double *w3_p;

	double *c,*r,*s;
	double *sr,*sl,*rr,*rl;

	double *d1r,*d2r,*d3r; 
	double *d1l,*d2l,*d3l;  
  
	double *eps2r,*eps4r,*eps2l,*eps4l;

	double *dw1r,*dw2r,*dw3r; 
	double *dw1l,*dw2l,*dw3l;

	double *h1r,*h2r,*h3r; 
	double *h1l,*h2l,*h3l; 

	double *f1,*f2,*f3;	  		
};

extern S_bagS *bagS;
extern S_bagF *bagF;


extern void     initializeFlow(void);
extern void     initializeStructure(void);
extern void     structure(void);
extern void     updateStructure(void);
extern void     writeStructure(double time, double u, double u_dot, double u_ddot);
extern void     readInputStructure(void);
extern void     readInputFlow(void);
extern _FLOAT  *allocate_1d_array(int ntot,char *str);
extern _FLOAT **allocate_2d_arrays(int imax,int jmax,char *str);
extern void	    deallocate_2d_arrays(_FLOAT **tmp);
extern int      bicgstab(int *iter,_FLOAT *resid,int n,char *matdescra,_FLOAT *A,int lval,int ndiag,int *dist,char *predescra,_FLOAT *b,_FLOAT *x);
extern void 	printVectorD(double* vect1,int dim);
extern void 	printMatrixD(double** mat1,int nl,int nc);
extern void 	createAmatrix(double **A);
extern void 	createBmatrix(double *B);
extern void 	writeFlow(double time, double *rho, double *u, double *e, double *p);
extern void 	flow_explicit(void);
extern void 	updateFlow(void);
extern void     printVectorD(double* vect1,int dim);
extern void 	printMatrixD(double** mat1,int nl,int nc);
extern void     writeMatrixD(double** mat1,int nl,int nc, char *name);
extern void solveFlowEquations(int ndiag, double **A, double *B);