#ifndef ACADOS_SIM_SIM_ERK_INTEGRATOR_H_
#define ACADOS_SIM_SIM_ERK_INTEGRATOR_H_

#include "acados/utils/types.h"

#define FIXED_STEP_SIZE 0

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
} sim_out;

void integrate(const sim_in *in, sim_out *out);

#endif  // ACADOS_SIM_SIM_ERK_INTEGRATOR_H_
