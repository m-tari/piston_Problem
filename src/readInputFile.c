#include "main_fsi.h"

void readInputFile(void)
{

	char line[100];
	char *filename = "Input.in";
	FILE *Input;
	Input = fopen (filename, "r");
    if(Input==NULL)
    {
    	printf("Error: Input file (%s) not found\n ",filename);
    }
	// reading input parameters from input file

    // common
	fgets(line,100, Input);	
	sscanf(line,"%d %d",&com->Ssolve, &com->Fsolve);
	printf(" Ssolve, Fsolve  = %d %d\n", com->Ssolve, com->Fsolve); 


	fgets(line,100, Input);	
	sscanf(line,"%lf ",&com->total_time);
	printf(" total_time      = %lf\n", com->total_time);

	// structure
	printf("\n ---------------- Structure solver parameters ----------------\n");

	fgets(line,100, Input);	
	sscanf(line,"%lf ",&bagS->mp);
	printf(" mp       = %lf\n", bagS->mp);   

	fgets(line,100, Input);	
	sscanf(line,"%lf ",&bagS->dt);
	printf(" dt       = %lf\n", bagS->dt);  

	fgets(line,100, Input);	
	sscanf(line,"%lf ",&bagS->kp);
	printf(" kp       = %lf\n", bagS->kp); 

	fgets(line,100, Input);	
	sscanf(line,"%lf ",&bagS->A);
	printf(" A        = %lf\n", bagS->A); 

	fgets(line,100, Input);	
	sscanf(line,"%lf ",&bagS->u_p);
	printf(" u_p      = %lf\n", bagS->u_p); 

	fgets(line,100, Input);	
	sscanf(line,"%lf ",&bagS->u_dot_p);
	printf(" u_dot_p  = %lf\n", bagS->u_dot_p);

	// flow
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


	fclose(Input);	
}