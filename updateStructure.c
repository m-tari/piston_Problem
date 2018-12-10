#include "main_fsi.h"

void updateStructure(void)
{
	bagS->u_p      = bagS->u;
	bagS->u_dot_p  = bagS->u_dot;
	bagS->u_ddot_p = bagS->u_ddot;	
}