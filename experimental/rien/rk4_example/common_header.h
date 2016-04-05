
#ifndef COMMON_HEADER
#define COMMON_HEADER

#include <math.h>

typedef double real_t;

/** Number of differential variables. */
#define NX  2
/** Number of control variables. */
#define NU  1

/** Fixed number of integration steps of fixed step size. */
#define FIXED_STEP_SIZE 0

#if FIXED_STEP_SIZE == 1
/** Integration step size. */
real_t H_INT = 1.0/100;
/** Number of integration steps. */
#define NSTEPS  10
#endif

typedef struct sim_in_s {
	real_t x[NX];
	real_t u[NU];

#if FIXED_STEP_SIZE == 0
	real_t step;
	unsigned int nSteps;
#endif
} sim_in;

typedef struct sim_out_s {
	real_t xn[NX];
	real_t Sx[NX*NX];
	real_t Su[NX*NU];
	real_t cpuTime;
} sim_out;

typedef unsigned int uint;

#endif /* COMMON_HEADER */
