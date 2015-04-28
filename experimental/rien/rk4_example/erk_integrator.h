
#ifndef ERK_INTEGRATOR_HEADER
#define ERK_INTEGRATOR_HEADER

#include "common_header.h"
#include "timing_functions.h"
extern void printMatrix( const char* name, const real_t* mat, uint nRows, uint nCols );


/** Fixed number of stages for the Explicit Runge-Kutta method. */
#define NUM_STAGES 4

/** Matrix of size: 4 x 4 (row major format) */
static const real_t A_MAT[ 16 ] =
{ 0,  0,   0, 0,
 0.5, 0,   0, 0,
 0,   0.5, 0, 0,
 0,   0,   1, 0 };

/** Vector of size: 1 x 4 */
static const real_t B_VEC[ 4 ] =
{ 1.0/6, 2.0/6, 2.0/6, 1.0/6 };

/** Vector of size: 1 x 4 */
static const real_t C_VEC[ 4 ] =
{ 0, 0.5, 0.5, 1.0 };

real_t vec[NX*(1+NX+NU)+NU];
real_t K_tmp[NUM_STAGES*NX*(1+NX+NU)];
real_t tmp[NX*(1+NX+NU)];

#endif /* ERK_INTEGRATOR_HEADER */
