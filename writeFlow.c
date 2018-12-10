#include "main_fsi.h"

void writeFlow(double time, double *rho, double *u, double *e, double *p)
{
	FILE *fout;

	if (time == 0.0)
	{
		fout = fopen("flow.dat", "w+");
		fprintf(fout, "\"time = %f\"\n", time);
		// fprintf(fout, "variables = rho, u, e, p\n");		
		for (int i = 0; i <= bagF->N; ++i)
		{
			fprintf(fout, "%d %.6f %.6f %.6f %.6f\n", i, bagF->rho[i], bagF->u[i], bagF->e[i], bagF->p[i]);
		}
		fprintf(fout, "\n\n");		
		fclose(fout);
	}
	else 
	{
		fout = fopen("flow.dat", "a+");
		fprintf(fout, "\"time = %f\"\n", time);
		// fprintf(fout, "variables = rho, u, e, p\n");						
		for (int i = 0; i <= bagF->N; ++i)
		{
			fprintf(fout, "%d %.6f %.6f %.6f %.6f\n", i, bagF->rho[i], bagF->u[i], bagF->e[i], bagF->p[i]);
		}
		fprintf(fout, "\n\n");		
		fclose(fout);
	}

	
}