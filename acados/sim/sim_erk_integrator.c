#include "acados/sim/sim_erk_integrator.h"
#include "acados/sim/model.h"

// Fixed number of stages for the Explicit Runge-Kutta method.
#define NUM_STAGES 4

#if FIXED_STEP_SIZE == 1
// Integration step size.
real_t H_INT = 1.0/100;
// Number of integration steps.
#define NSTEPS  10
#endif

// Butcher tableau (row major format)
static const real_t A_MAT[ 16 ] = {0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 1, 0};
static const real_t B_VEC[ 4 ] = {1.0/6, 2.0/6, 2.0/6, 1.0/6};

real_t vec[NX*(1+NX+NU)+NU];
real_t K_tmp[NUM_STAGES*NX*(1+NX+NU)];
real_t tmp[NX*(1+NX+NU)];

void integrate(const sim_in *in, sim_out *out) {
    unsigned int i, s, j, istep;
#if FIXED_STEP_SIZE == 0
    real_t H_INT = in->step;
    unsigned int NSTEPS = in->nSteps;
#endif

    for (i = 0; i < NX; i++) tmp[i] = in->x[i];
    for (i = 0; i < NX*(NX+NU); i++) tmp[NX+i] = 0.0;  // sensitivities
    for (i = 0; i < NX; i++) tmp[NX+i*NX+i] = 1.0;     // sensitivities wrt x

    for (i = 0; i < NU; i++) vec[NX*(1+NX+NU)+i] = in->u[i];

    for (istep = 0; istep < NSTEPS; istep++) {
        for (s = 0; s < NUM_STAGES; s++) {
            for (i = 0; i < NX*(1+NX+NU); i++) {
                vec[i] = tmp[i];
            }
            for (j = 0; j < s; j++) {
                if (A_MAT[s*NUM_STAGES+j] != 0) {
                    for (i = 0; i < NX*(1+NX+NU); i++) {
                        vec[i] += H_INT*A_MAT[s*NUM_STAGES+j]*K_tmp[j*NX*(1+NX+NU)+i];
                    }
                }
            }
            VDE_fun(vec, &(K_tmp[s*NX*(1+NX+NU)]));  // k evaluation
        }
        for (s = 0; s < NUM_STAGES; s++) {
            for (i = 0; i < NX*(1+NX+NU); i++) {
                tmp[i] += H_INT*B_VEC[s]*K_tmp[s*NX*(1+NX+NU)+i];  // ERK step
            }
        }
    }
    for (i = 0; i < NX; i++)      out->xn[i] = tmp[i];
    for (i = 0; i < NX*NX; i++) out->Sx[i] = tmp[NX+i];
    for (i = 0; i < NX*NU; i++) out->Su[i] = tmp[NX+NX*NX+i];
}
