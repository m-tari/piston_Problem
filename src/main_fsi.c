#include "main_fsi.h"


struct S_bagS *bagS;
struct S_bagF *bagF;
struct S_com  *com;
 
int main(int argc, char const *argv[])
{
	_FLOAT time       = 0.0;

	bagS = malloc(sizeof(S_bagS));
	bagF = malloc(sizeof(S_bagF));	
	com  = malloc(sizeof(S_com));

	readInputFile();	

	if (com->Ssolve == 1)
	{
		initializeStructure();
		writeStructure(time, bagS->u_p, bagS->u_dot_p, bagS->u_ddot_p);
	}
	if (com->Fsolve == 1)
	{
		initializeFlow();
		writeFlow(time, bagF->rho, bagF->u, bagF->e, bagF->p);
	}

	printf("\nSolving...\n");
	while (time < com->total_time)
	{

		printf("time = %f\n", time);
			
		/* structue solver */
		if (com->Ssolve == 1)
		{
			structure();
			writeStructure(time, bagS->u, bagS->u_dot, bagS->u_ddot);
			updateStructure();
			time = time + bagS->dt;
		}
		/* flow solver */
		if (com->Fsolve == 1)
		{
			flow_explicit();
			writeFlow(time, bagF->rho, bagF->u, bagF->e, bagF->p);
			updateFlow();
			time = time + bagF->dt;
		}
	}	
	return 0;
}