
#include "common_header.h"
#include "timing_functions.h"
extern void printMatrix( const char* name, const real_t* mat, uint nRows, uint nCols );

real_t vec[NX*(1+NX+NU)+NU];
real_t k1[NX*(1+NX+NU)];
real_t k2[NX*(1+NX+NU)];
real_t k3[NX*(1+NX+NU)];
real_t k4[NX*(1+NX+NU)];
real_t tmp[NX*(1+NX+NU)];

int integrate( const sim_in* in, sim_out* out ) {

	unsigned int i, istep;
	timer tmr;
#if FIXED_STEP_SIZE == 0
	real_t H_INT = in->step;
	unsigned int NSTEPS = in->nSteps;
#endif

	for( i = 0; i < NX; i++ ) tmp[i] = in->x[i];
	for( i = 0; i < NX*(NX+NU); i++ ) tmp[NX+i] = 0.0; // sensitivities
	for( i = 0; i < NX; i++ ) tmp[NX+i*NX+i] = 1.0;    // sensitivities wrt x
	
	for( i = 0; i < NU; i++ ) vec[NX*(1+NX+NU)+i] = in->u[i];
	
	tic(&tmr);
	for( istep = 0; istep < NSTEPS; istep++ ) {
//		printMatrix("cur_x", tmp, 1, NX);
		for( i = 0; i < NX*(1+NX+NU); i++ ) vec[i] = tmp[i];
		VDE_fun( vec, k1 );  // k1
		
		for( i = 0; i < NX*(1+NX+NU); i++ ) vec[i] = tmp[i] + H_INT*0.5*k1[i];
		VDE_fun( vec, k2 );  // k2
		
		for( i = 0; i < NX*(1+NX+NU); i++ ) vec[i] = tmp[i] + H_INT*0.5*k2[i];
		VDE_fun( vec, k3 );  // k3
		
		for( i = 0; i < NX*(1+NX+NU); i++ ) vec[i] = tmp[i] + H_INT*k3[i];
		VDE_fun( vec, k4 );  // k4
		
		for( i = 0; i < NX*(1+NX+NU); i++ ) tmp[i] += H_INT/6*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]); // RK4 step
	}
	out->cpuTime = toc(&tmr);
	for( i = 0; i < NX; i++ ) out->xn[i] = tmp[i];
	for( i = 0; i < NX*NX; i++ ) out->Sx[i] = tmp[NX+i];
	for( i = 0; i < NX*NU; i++ ) out->Su[i] = tmp[NX+NX*NX+i];
	
	return 0;
}
