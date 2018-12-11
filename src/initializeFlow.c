#include "main_fsi.h"

void initializeFlow(void)
{

	int numAlloc = bagF->N + 1;
	bagF->rho       = allocate_1d_array(numAlloc, "rho");
	bagF->u       	= allocate_1d_array(numAlloc, "u");
	bagF->e         = allocate_1d_array(numAlloc, "e");
	bagF->p         = allocate_1d_array(numAlloc, "p");
	
	bagF->x         = allocate_1d_array(3*bagF->N, "x");


	bagF->H         = allocate_1d_array(numAlloc, "H");
	bagF->w1         = allocate_1d_array(numAlloc, "w1");
	bagF->w2         = allocate_1d_array(numAlloc, "w2");
	bagF->w3         = allocate_1d_array(numAlloc, "w3");
	bagF->w1_p       = allocate_1d_array(numAlloc, "w1_p");
	bagF->w2_p       = allocate_1d_array(numAlloc, "w2_p");
	bagF->w3_p       = allocate_1d_array(numAlloc, "w3_p");

	bagF->c         = allocate_1d_array(numAlloc, "c");
	bagF->r         = allocate_1d_array(numAlloc, "r");
	bagF->s         = allocate_1d_array(numAlloc, "s");

	bagF->rr         = allocate_1d_array(numAlloc, "rr");
	bagF->rl         = allocate_1d_array(numAlloc, "rl");
	bagF->sr         = allocate_1d_array(numAlloc, "sr");
	bagF->sl         = allocate_1d_array(numAlloc, "sl");

	bagF->d1r         = allocate_1d_array(numAlloc, "d1r");
	bagF->d2r         = allocate_1d_array(numAlloc, "d2r");
	bagF->d3r         = allocate_1d_array(numAlloc, "d3r");
	bagF->d1l         = allocate_1d_array(numAlloc, "d1l");
	bagF->d2l         = allocate_1d_array(numAlloc, "d2l");
	bagF->d3l         = allocate_1d_array(numAlloc, "d3l");
	bagF->dw1r         = allocate_1d_array(numAlloc, "dw1r");
	bagF->dw2r         = allocate_1d_array(numAlloc, "dw2r");
	bagF->dw3r         = allocate_1d_array(numAlloc, "dw3r");
	bagF->dw1l         = allocate_1d_array(numAlloc, "dw1l");
	bagF->dw2l         = allocate_1d_array(numAlloc, "dw2l");
	bagF->dw3l         = allocate_1d_array(numAlloc, "dw3l");

	bagF->eps2r         = allocate_1d_array(numAlloc, "eps2r");
	bagF->eps2l         = allocate_1d_array(numAlloc, "eps2r");
	bagF->eps4r         = allocate_1d_array(numAlloc, "eps2r");
	bagF->eps4l         = allocate_1d_array(numAlloc, "eps2r");

	bagF->h1r         = allocate_1d_array(numAlloc, "h1r");
	bagF->h2r         = allocate_1d_array(numAlloc, "h2r");
	bagF->h3r         = allocate_1d_array(numAlloc, "h3r");
	bagF->h1l         = allocate_1d_array(numAlloc, "h1l");
	bagF->h2l         = allocate_1d_array(numAlloc, "h2l");
	bagF->h3l         = allocate_1d_array(numAlloc, "h3l");

	bagF->f1         = allocate_1d_array(numAlloc, "f1");
	bagF->f2         = allocate_1d_array(numAlloc, "f2");
	bagF->f3         = allocate_1d_array(numAlloc, "f3");

	bagF->dx = bagF->L/(bagF->N-1);
	_FLOAT R = 1.0;//287.0; // gas constant for the air
	_FLOAT cv = R/(gam-1.0);
	_FLOAT T = 1/gam; //25 + 273.15;

	/* nodes on the fluid solver starts from 0 to N, where 0 and N are BCs */
	for (int i = 0; i <= bagF->N; ++i)
	{
		bagF->u[i] = 0.0;
		bagF->e[i] = cv*T; // without T for testing
	}

	for (int i = 0; i <= bagF->N; ++i)
	{
		if (((_FLOAT)i/bagF->N)*bagF->dx < bagF->dx/2.0)
		{
			bagF->p[i] = 10.0/gam;
			bagF->rho[i] = 10.0;
			bagF->H[i] = bagF->e[i]+bagF->p[i]/bagF->rho[i];

		}
		else
		{
			bagF->p[i] = 1.0/gam;
			bagF->rho[i] = 1.0;
			bagF->H[i] = bagF->e[i]+bagF->p[i]/bagF->rho[i];			
		}
	}

	// initialize w and w_p	
	for (int i = 0; i <= bagF->N; ++i)
	{
		bagF->w1[i] = bagF->rho[i];
		bagF->w2[i] = bagF->rho[i]*bagF->u[i];
		bagF->w3[i] = bagF->rho[i]*bagF->e[i];

		bagF->w1_p[i] = bagF->w1[i];
		bagF->w2_p[i] = bagF->w2[i];
		bagF->w3_p[i] = bagF->w3[i];	
	}
}