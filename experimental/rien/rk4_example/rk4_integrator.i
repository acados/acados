%module rk4_integrator
%include typemaps.i

%{
 /* Put header files here or function declarations like below */
#include "common_header.h"
%}

extern int integrate( const sim_in* in, sim_out* out );
extern void printMatrix( const char* name, real_t* mat, uint nRows, uint nCols );

%include "common_header.h"		// Just grab original C header file
%extend sim_in_s {
	sim_in_s() {
		uint i;
		sim_in *s;
		s = (sim_in*) malloc(sizeof(sim_in));
		for( i = 0; i < NX; i++ ) s->x[i] = 0.0;
		for( i = 0; i < NU; i++ ) s->u[i] = 0.0;
		
		#if FIXED_STEP_SIZE == 0
			s->step = 0.1;
			s->nSteps = 1;
		#endif
		return s;
	}
};


