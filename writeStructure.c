#include "main_fsi.h"

void writeStructure(double time, double u, double u_dot, double u_ddot)
{
	FILE *fout;

	if (time == 0.0)
	{
		fout = fopen("structure.dat", "w+");
		fprintf(fout, "TITLE = \"Structure solver displceament, velocity and acceleration\"\n");
		fprintf(fout, "Variables = t, u, u_dot, u_ddot\n");
		fclose(fout);
	}

	fout = fopen("structure.dat", "a+");
	fprintf(fout, "%f %.10f %.10f %.10f\n", time, u, u_dot, u_ddot);
	fclose(fout);
	
}