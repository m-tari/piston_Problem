#include "preproc.h"

/******************************************************/
/*	This function deallocate the memory for a 2-D      */
/*	array of a given size, all contiguous              */
/******************************************************/

void deallocate_2d_arrays(_FLOAT **tmp)
{
  free(&tmp[0][0]); 
  free(&tmp[0]); tmp = NULL;
  return;
}
