#include "main_fsi.h"
#include "preproc.h"

void flow_explicit(void)
{
	double k2 = 1.; 
	double k4= 1./32.;
	double c4=2.;

	double *d1r   = bagF->d1r;
	double *d2r   = bagF->d2r;
	double *d3r   = bagF->d3r;
	double *d1l   = bagF->d1l;
	double *d2l   = bagF->d2l;
	double *d3l   = bagF->d3l;

	double *eps2r = bagF->eps2r;
	double *eps4r = bagF->eps4r;
	double *eps2l = bagF->eps2l;
	double *eps4l = bagF->eps4l;

	double *dw1r  = bagF->dw1r;
	double *dw2r  = bagF->dw2r;
	double *dw3r  = bagF->dw3r;
	double *dw1l  = bagF->dw1l;
	double *dw2l  = bagF->dw2l;
	double *dw3l  = bagF->dw3l;


	// r = abs(u)+c, s=f(p), c=sqrt(g*p/rho)
	for (int i = 1; i <= bagF->N-1; ++i)
	{
		bagF->c[i] = sqrt(gam*bagF->p[i]/bagF->rho[i]);
		bagF->r[i] = fabs(bagF->u[i]) + bagF->c[i];
		bagF->s[i] = fabs((bagF->p[i+1]-2.*bagF->p[i]+bagF->p[i-1])/(bagF->p[i+1]+2.*bagF->p[i]+bagF->p[i-1]));
	}

	//sr, sl, rr, rl (right, left) of s and r in JST
	for (int i = 2; i < bagF->N-1; ++i)
	{	
		bagF->sr[i] = max(bagF->s[i],  bagF->s[i+1]);
		bagF->sl[i] = max(bagF->s[i-1],bagF->s[i]);
		bagF->rr[i] = max(bagF->r[i+1],bagF->r[i]);
		bagF->rl[i] = max(bagF->r[i],  bagF->r[i-1]);

		bagF->eps2r[i]  = k2*bagF->sr[i]*bagF->rr[i];
		bagF->eps4r[i]  = max(0, k4*bagF->rr[i]-c4*bagF->eps2r[i]);

		bagF->eps2l[i]  = k2*bagF->sl[i]*bagF->rl[i];
		bagF->eps4l[i]  = max(0, k4*bagF->rl[i]-c4*bagF->eps2l[i]);				
	}

	// eps2 = f(k2,s,r) , eps4=f(k4,r,c4,eps2)
	for (int i = 1; i <= bagF->N-1; ++i)
	{	
		// right
		bagF->dw1r[i]   = bagF->w1[i+1]-bagF->w1[i];
		bagF->dw2r[i]   = bagF->w2[i+1]-bagF->w2[i];
		bagF->dw3r[i]   = bagF->w3[i+1]-bagF->w3[i];
		// left
		bagF->dw1l[i]   = bagF->w1[i]-bagF->w1[i-1];
		bagF->dw2l[i]   = bagF->w2[i]-bagF->w2[i-1];
		bagF->dw3l[i]   = bagF->w3[i]-bagF->w3[i-1];		
	}

	// compute d by eps 2 and 4 and dw	
	for (int i = 2; i < bagF->N-1; ++i)
	{	
		// right	
		d1r[i] = 0.0;//eps2r[i]*dw1r[i]-eps4r[i]*(dw1r[i+1]-2.*dw1r[i]+dw1r[i-1]);
		d2r[i] = 0.0;//eps2r[i]*dw2r[i]-eps4r[i]*(dw2r[i+1]-2.*dw2r[i]+dw2r[i-1]);
		d3r[i] = 0.0;//eps2r[i]*dw3r[i]-eps4r[i]*(dw3r[i+1]-2.*dw3r[i]+dw3r[i-1]);
		// left
		d1l[i] = 0.0;//eps2l[i]*dw1l[i]-eps4l[i]*(dw1l[i+1]-2.*dw1l[i]+dw1l[i-1]);
		d2l[i] = 0.0;//eps2l[i]*dw2l[i]-eps4l[i]*(dw2l[i+1]-2.*dw2l[i]+dw2l[i-1]);
		d3l[i] = 0.0;//eps2l[i]*dw3l[i]-eps4l[i]*(dw3l[i+1]-2.*dw3l[i]+dw3l[i-1]);
	}

	// convective fluxes
	for (int i = 2; i < bagF->N-1; ++i)
	{
		bagF->f1[i] = bagF->rho[i]*bagF->u[i];
		bagF->f2[i] = bagF->rho[i]*pow(bagF->u[i],2)+bagF->p[i];
		bagF->f3[i] = bagF->rho[i]*bagF->u[i]*bagF->H[i];
	}

	// total flux (fluid + artificial)
	for (int i = 2; i < bagF->N-1; ++i)
	{	
		// right	
		bagF->h1r[i] = 0.5*(bagF->f1[i+1]+bagF->f1[i])-d1r[i];
		bagF->h2r[i] = 0.5*(bagF->f2[i+1]+bagF->f2[i])-d2r[i];
		bagF->h3r[i] = 0.5*(bagF->f3[i+1]+bagF->f3[i])-d3r[i];
		// if (i==bagF->N-2) // flux at right wall bc
		// {
		// 	bagF->h1r[i] = 0.0        - d1r[i];
		// 	bagF->h2r[i] = bagF->p[i] - d2r[i];
		// 	bagF->h3r[i] = 0.0        - d3r[i];
		// }		
		// left
		bagF->h1l[i] = 0.5*(bagF->f1[i]+bagF->f1[i-1])-d1l[i];
		bagF->h2l[i] = 0.5*(bagF->f2[i]+bagF->f2[i-1])-d2l[i];
		bagF->h3l[i] = 0.5*(bagF->f3[i]+bagF->f3[i-1])-d3l[i];
		if (i==2)         // flux at left wall bc
		{
			bagF->h1l[i] = 0.0        - d1l[i];
			bagF->h2l[i] = bagF->p[i] - d2l[i];
			bagF->h3l[i] = 0.0        - d3l[i];
		}				
	}

	for (int i = 2; i < bagF->N-1; ++i)
	{
		// compute wn+1 by wn(known) and h(flux+dissipative terms) and dx(known) - 1st order euler
		bagF->w1[i] = -(bagF->dt/bagF->dx)*(bagF->h1r[i] - bagF->h1l[i]) + bagF->w1_p[i];
		bagF->w2[i] = -(bagF->dt/bagF->dx)*(bagF->h2r[i] - bagF->h2l[i]) + bagF->w2_p[i];
		bagF->w3[i] = -(bagF->dt/bagF->dx)*(bagF->h3r[i] - bagF->h3l[i]) + bagF->w3_p[i];
	}

	// update halos (left)
	bagF->w1[1] = 2.*bagF->w1[2] - bagF->w1[3]; 
	bagF->w2[1] = 2.*bagF->w2[2] - bagF->w2[3]; 
	bagF->w3[1] = 2.*bagF->w3[2] - bagF->w3[3]; 

	bagF->w1[0] = 3.*bagF->w1[2] - 2.*bagF->w1[3]; 
	bagF->w2[0] = 3.*bagF->w2[2] - 2.*bagF->w2[3]; 
	bagF->w3[0] = 3.*bagF->w3[2] - 2.*bagF->w3[3];

	// update halos (right)
	bagF->w1[bagF->N-1] = 2.*bagF->w1[bagF->N-2] - bagF->w1[bagF->N-3]; 
	bagF->w2[bagF->N-1] = 2.*bagF->w2[bagF->N-2] - bagF->w2[bagF->N-3]; 
	bagF->w3[bagF->N-1] = 2.*bagF->w3[bagF->N-2] - bagF->w3[bagF->N-3]; 

	bagF->w1[bagF->N]   = 3.*bagF->w1[bagF->N-2] - 2.*bagF->w1[bagF->N-3]; 
	bagF->w2[bagF->N]   = 3.*bagF->w2[bagF->N-2] - 2.*bagF->w2[bagF->N-3]; 
	bagF->w3[bagF->N]   = 3.*bagF->w3[bagF->N-2] - 2.*bagF->w3[bagF->N-3];

	// convert w1,w2,w3 to rho,u,p and e
	for (int i = 0; i <= bagF->N; ++i)
	{
		bagF->rho[i] = bagF->w1[i];
		bagF->u[i]   = bagF->w2[i]/bagF->w1[i];
		bagF->e[i]   = bagF->w3[i]/bagF->w1[i];

		bagF->p[i]   = (gam-1.)*bagF->rho[i]*(bagF->e[i]-pow(bagF->u[i],2)/2.);
		bagF->H[i]   = bagF->e[i]+bagF->p[i]/bagF->rho[i];
	}			
}