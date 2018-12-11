#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define _FLOAT double 
#define  gam   1.4

typedef struct S_bagS S_bagS;
typedef struct S_bagF S_bagF;
typedef struct S_com S_com;

struct S_bagS
{
	_FLOAT mp;  // piston mass
	_FLOAT dt;  // structure timestep size
	_FLOAT kp;  // spring stifness
	_FLOAT A;   // piston area

	_FLOAT u;
	_FLOAT u_dot;
	_FLOAT u_ddot;

	_FLOAT u_p;
	_FLOAT u_dot_p;
	_FLOAT u_ddot_p;	
};

struct S_bagF
{
	_FLOAT p0;  // pressure at chamber at rest

	_FLOAT dt;
	_FLOAT dx;
	_FLOAT L;

	_FLOAT q1;
	_FLOAT q2;
	_FLOAT q3;

	_FLOAT *x;
	_FLOAT *rho;
	_FLOAT *u;
	_FLOAT *e;
	_FLOAT *p;
	_FLOAT *H;	
	int N;

	_FLOAT *w1;
	_FLOAT *w2;
	_FLOAT *w3;
	_FLOAT *w1_p;
	_FLOAT *w2_p;
	_FLOAT *w3_p;

	_FLOAT *c,*r,*s;
	_FLOAT *sr,*sl,*rr,*rl;

	_FLOAT *d1r,*d2r,*d3r; 
	_FLOAT *d1l,*d2l,*d3l;  
  
	_FLOAT *eps2r,*eps4r,*eps2l,*eps4l;

	_FLOAT *dw1r,*dw2r,*dw3r; 
	_FLOAT *dw1l,*dw2l,*dw3l;

	_FLOAT *h1r,*h2r,*h3r; 
	_FLOAT *h1l,*h2l,*h3l; 

	_FLOAT *f1,*f2,*f3;	  		
};

struct S_com
{
	int Ssolve;  // switch for structure solver (0-off,1-on)
	int Fsolve;  // switch for flow solver (0-oof,1-on)
	_FLOAT total_time;	

};

  
extern S_bagS *bagS;
extern S_bagF *bagF;
extern S_com  *com;

extern void     initializeFlow(void);
extern void     initializeStructure(void);
extern void     structure(void);
extern void     updateStructure(void);
extern void     writeStructure(_FLOAT time, _FLOAT u, _FLOAT u_dot, _FLOAT u_ddot);
extern void     readInputFile(void);
extern _FLOAT  *allocate_1d_array(int ntot,char *str);
extern _FLOAT **allocate_2d_arrays(int imax,int jmax,char *str);
extern void	    deallocate_2d_arrays(_FLOAT **tmp);
extern int      bicgstab(int *iter,_FLOAT *resid,int n,char *matdescra,_FLOAT *A,int lval,int ndiag,int *dist,char *predescra,_FLOAT *b,_FLOAT *x);
extern void 	printVectorD(_FLOAT* vect1,int dim);
extern void 	printMatrixD(_FLOAT** mat1,int nl,int nc);
extern void 	createAmatrix(_FLOAT **A);
extern void 	createBmatrix(_FLOAT *B);
extern void 	writeFlow(_FLOAT time, _FLOAT *rho, _FLOAT *u, _FLOAT *e, _FLOAT *p);
extern void 	flow_explicit(void);
extern void 	updateFlow(void);
extern void     printVectorD(_FLOAT* vect1,int dim);
extern void 	printMatrixD(_FLOAT** mat1,int nl,int nc);
extern void     writeMatrixD(_FLOAT** mat1,int nl,int nc, char *name);
extern void     solveFlowEquations(int ndiag, _FLOAT **A, _FLOAT *B);