#include "main_fsi.h"

void createBmatrix(double *B)
{
	int N = bagF->N;	
	/* Boundary condition for i=0 and i=N nodes */
	/* i=1 */
	B[0] = 0.0; 
	B[1] = 0.0; 
	B[2] = 0.0; 
	/* i=N */
	B[3*N-3] = 0.0; 
	B[3*N-2] = 0.0; 
	B[3*N-1] = 0.0; 

	for (int i = 2; i <= bagF->N-1; ++i)
	{
		B[3*(i-1)]   = -(bagF->dt/(2.0*bagF->dx))*(bagF->rho[i+1]*bagF->u[i+1]-bagF->rho[i-1]*bagF->u[i-1]);
		B[3*(i-1)+1] = -(bagF->dt/(2.0*bagF->dx))*((bagF->rho[i+1]*pow(bagF->u[i+1],2)+bagF->p[i+1])-(bagF->rho[i-1]*pow(bagF->u[i-1],2)+bagF->p[i-1]));
		B[3*(i-1)+2] = -(bagF->dt/(2.0*bagF->dx))*(((bagF->e[i+1]+bagF->p[i+1])*bagF->u[i+1])-((bagF->e[i-1]+bagF->p[i-1])*bagF->u[i-1]));
	}
}