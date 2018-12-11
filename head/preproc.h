#ifndef HEAD_PREPROC
#define HEAD_PREPROC

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <unistd.h>
#include <time.h>
#include "mkl.h"

#define STLEN 128
#define _FLOAT double

#define _SQRT sqrt
#define _FABS fabs
#define _POW pow
#define max(a,b) ((a)>(b) ? (a):(b))
#define min(a,b) ((a)<(b) ? (a):(b))
#define sign(a)  ((a)>=(0) ? (1):(-1))
#define trig(a)  ((a)>(0.0) ? (1):(0))
#define fsign(a) ((a)>=(0.0) ? (1.0):(-1.0))
#define ftrig(a) ((a)>(0.0) ? (1.0):(0.0))

#define PI   (4.*atan(1.))
#define INF  1.0E300 /*Highest number available with double precision*/

#endif

