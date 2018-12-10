#include "preproc.h"

/******************************************************/
/*	This function allocate the memory for a 2-D _FLOAT */
/*	array of a given size, all contiguous              */
/******************************************************/

_FLOAT **allocate_2d_arrays(int imax,int jmax,char *str)
{
  int i;
  char mess[STLEN];
  _FLOAT *tmpij,**tmp;
  
  
  /*printf("allocate_2d_arrays of size %d %d for %s\n",imax,jmax,str);*/
  tmpij=(_FLOAT *)calloc((imax*jmax),sizeof(_FLOAT));
  if(tmpij==NULL)
  {
    sprintf(mess,"%s can not be allocated",str);
  }
  //for (i=0;i<imax*jmax;i++) tmpij[i]=0.;
  
  tmp=(_FLOAT **)calloc(imax,sizeof(_FLOAT *));
  if(tmp==NULL)
  {
    sprintf(mess,"%s can not be allocated",str);
  }

  for (i=0;i<imax;i++) tmp[i]=&(tmpij[i*jmax]);

  return tmp;
}
