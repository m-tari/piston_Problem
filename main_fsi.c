#include "main_fsi.h"


struct S_bagS *bagS;
struct S_bagF *bagF;

int main(int argc, char const *argv[])
{
	double time       = 0.0;

	initializeFlow();
	writeFlow(time, bagF->rho, bagF->u, bagF->e, bagF->p);

	// initializeStructure();
	// writeStructure(time, bagS->u_p, bagS->u_dot_p, bagS->u_ddot_p);


	while (time < bagF->total_time)
	{

		printf("time = %f\n", time);	
		/* structue solver */
		// structure();

		/* flow solver */
		flow_explicit();
		// time = time + bagS->dt;
		time = time + bagF->dt;
		// writeStructure(time, bagS->u, bagS->u_dot, bagS->u_ddot);

		// updateStructure();
		updateFlow();

		writeFlow(time, bagF->rho, bagF->u, bagF->e, bagF->p);		

	}	
	return 0;
}