
#include "erk_integrator.h"

int integrate( const sim_in* in, sim_out* out ) {

	unsigned int i, s, j, istep;
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
		for( s = 0; s < NUM_STAGES; s++ ) {
			for( i = 0; i < NX*(1+NX+NU); i++ ) vec[i] = tmp[i];
			for( j = 0; j < s; j++ ) {
				if( A_MAT[s*NUM_STAGES+j] != 0 ) {
//					printMatrix("A_MAT[s*NUM_STAGES+j]", &A_MAT[s*NUM_STAGES+j], 1, 1);
					for( i = 0; i < NX*(1+NX+NU); i++ ) vec[i] += H_INT*A_MAT[s*NUM_STAGES+j]*K_tmp[j*NX*(1+NX+NU)+i];
				}
			}
			VDE_fun( vec, &(K_tmp[s*NX*(1+NX+NU)]) );  // k evaluation
		}
		for( s = 0; s < NUM_STAGES; s++ ) {
			for( i = 0; i < NX*(1+NX+NU); i++ ) tmp[i] += H_INT*B_VEC[s]*K_tmp[s*NX*(1+NX+NU)+i]; // ERK step
		}
	}
	out->cpuTime = toc(&tmr);
	for( i = 0; i < NX; i++ ) out->xn[i] = tmp[i];
	for( i = 0; i < NX*NX; i++ ) out->Sx[i] = tmp[NX+i];
	for( i = 0; i < NX*NU; i++ ) out->Su[i] = tmp[NX+NX*NX+i];
	
	return 0;
}
