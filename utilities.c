#include "preproc.h"
#include "main_fsi.h"

/* Function that prints double Vector */
void printVectorD(double* vect1,int dim)
{
    int i;
    printf(" \n");
    for (i = 0 ; i < dim ; i++) {
        printf(" %.9f \n", vect1[i]);
    }
}

/* Function that prints double Matrix */
void printMatrixD(double** mat1,int nl,int nc)
{
    int i,j;
    printf("\n");
    for(i=0;i<nl;i++)
    {
        for(j=0;j<nc;j++) 
        {
            printf(" %.4f",mat1[i][j]);         
        }
        printf("\n");
    }
    printf("\n");
}

/* Function that writes double Matrix */
void writeMatrixD(double** mat1,int nl,int nc, char *name)
{
    FILE *filename;
    filename = fopen(name, "w+");
    int i,j;
    for(i=0;i<nl;i++)
    {
        for(j=0;j<nc;j++) 
        {
            fprintf(filename, " %+.4f",mat1[i][j]);         
        }
        fprintf(filename, "\n");
    }
    fclose(filename);
}
