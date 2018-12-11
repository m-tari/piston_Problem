#include "main_fsi.h"

void initializeStructure(void)
{
	/* u_ddot is calculated base on u and u_dot initial condition from equation of motion */
	bagS->u_ddot_p = (1.0/bagS->mp)*(-bagS->kp*bagS->u_p + bagS->A*(bagF->p[bagF->N-1] - bagF->p0));
}