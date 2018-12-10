#include "preproc.h"

/******************************************************/
/*	This function allocate the memory for a 1-D        */
/*	array of a given size                              */
/******************************************************/

_FLOAT *allocate_1d_array(int nbpt,char *str)
{
   // int i;
   char mess[STLEN];
   _FLOAT *tmp;  

   // tmp=(_FLOAT *)malloc(nbpt*sizeof(_FLOAT));
   tmp=(_FLOAT *)calloc(nbpt,sizeof(_FLOAT));

   if(tmp==NULL)
   {
      sprintf(mess,"%s can not be allocated",str);
   }

   // for (i=0;i<nbpt;i++) tmp[i]=0.;

   return tmp;
}

