#include "main_fsi.h"

void readInputFlow(void)
{
	char line[100];
	char *filename = "InputFlow.in";
	FILE *Input;
	Input = fopen (filename, "r");
    if(Input==NULL)
    {
    	printf("Error: Input file (%s) not found\n ",filename);
    }
	// reading input parameters from input file
	printf("\n ---------------- Flow solver parameters ----------------\n");
	fgets(line,100, Input);	
	sscanf(line,"%lf ",&bagF->p0);
	printf(" p0         = %lf\n", bagF->p0);   

	fgets(line,100, Input);	
	sscanf(line,"%d ",&bagF->N);
	printf(" N          = %d\n", bagF->N);

	bagF->p  =  allocate_1d_array(bagF->N, "bagF->p");
	fgets(line,100, Input);	
	sscanf(line,"%lf ",&bagF->p[bagF->N-1]);
	printf(" p          = %lf\n", bagF->p[bagF->N-1]);	

	fgets(line,100, Input);	
	sscanf(line,"%lf ",&bagF->L);
	printf(" L          = %lf\n", bagF->L);

	fgets(line,100, Input);	
	sscanf(line,"%lf ",&bagF->dt);
	printf(" dt         = %lf\n", bagF->dt);

	fgets(line,100, Input);	
	sscanf(line,"%lf ",&bagF->total_time);
	printf(" total_time = %lf\n", bagF->total_time);

	fclose(Input);
			 	 		 
}