#include "main_fsi.h"

void updateFlow(void) 
{
	/* for implicit. commented for now
	double q1[bagF->N+1], q2[bagF->N+1], q3[bagF->N+1];

	for (int i = 0; i < bagF->N; ++i)
	{
		q1[i+1] = bagF->rho[i+1];
		q1[i+1] = q1[i+1] + bagF->x[3*i];

		q2[i+1] = bagF->rho[i+1]*bagF->u[i+1];
		q2[i+1] = q2[i+1] + bagF->x[3*i+1];

		q3[i+1] = bagF->rho[i+1]*bagF->e[i+1];
		q3[i+1] = q3[i+1] + bagF->x[3*i+2];
	}

	for (int i = 1; i <= bagF->N; ++i)
	{
		bagF->rho[i] = q1[i];
		bagF->u[i]   = q2[i]/q1[i];
		bagF->e[i]   = q3[i]/q1[i];
	}

	//  calculating pressure on the fluid nodes, equation of state  
	for (int i = 1; i <= bagF->N; ++i)
	{
		bagF->p[i] = (gam - 1.0)*(bagF->rho[i]*bagF->e[i] - 0.5*bagF->rho[i]*pow(bagF->u[i],2));
	}

	*/	
	for (int i = 0; i <= bagF->N; ++i)
	{
		bagF->w1_p[i] = bagF->w1[i];
		bagF->w2_p[i] = bagF->w1[i];
		bagF->w3_p[i] = bagF->w1[i];
	}		
}