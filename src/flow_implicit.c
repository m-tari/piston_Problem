#include "main_fsi.h"
#include "preproc.h"

void flow_implicit(void)
{
	int ndiag       = 11;
	_FLOAT **A      = allocate_2d_arrays(ndiag,3*bagF->N, "A");
	_FLOAT *B       = allocate_1d_array(3*bagF->N, "B");


	/* create coefficients matrix A */
	createAmatrix(A);

	/* create history term matrix B */
	createBmatrix(B);

	/* solve set of equations Ax=B */
	solveFlowEquations(ndiag, A, B);

	deallocate_2d_arrays(A);
	free(B);
	// free(dist);

}