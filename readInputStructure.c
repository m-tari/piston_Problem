#include "main_fsi.h"

void readInputStructure(void)
{
	char line[100];
	char *filename = "InputStructure.in";
	FILE *Input;
	Input = fopen (filename, "r");
    if(Input==NULL)
    {
    	printf("Error: Input file (%s) not found\n ",filename);
    }
	// reading input parameters from input file
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

	fclose(Input); 			 	 		 
}