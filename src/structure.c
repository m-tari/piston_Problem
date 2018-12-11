#include "main_fsi.h"


void structure(void)
{
	_FLOAT du;
	_FLOAT denom; 
	/* compute du */
	denom        = (4.0*bagS->mp/pow(bagS->dt,2) + bagS->kp);
	du           = (bagS->A*(bagF->p[bagF->N-1]-bagF->p0) - bagS->kp*bagS->u_p + bagS->mp*((4.0/bagS->dt)*bagS->u_dot_p + bagS->u_ddot_p))/denom;

	/* update u */
	bagS->u      = bagS->u_p + du;

	/* compute u_ddot */
	bagS->u_ddot = (4.0/pow(bagS->dt,2))*du - (4.0/bagS->dt)*bagS->u_dot_p - bagS->u_ddot_p;

	/* compute u_dot */
	bagS->u_dot  = bagS->u_dot_p + (bagS->dt/2.0)*(bagS->u_ddot_p + bagS->u_ddot);
}

