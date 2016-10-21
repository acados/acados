#ifndef ACADOS_SIM_SIM_COMMON_H_
#define ACADOS_SIM_SIM_COMMON_H_

#include "acados/utils/types.h"

#define FIXED_STEP_SIZE 0

typedef struct sim_in_ {
    int_t nx;   // NX
    int_t nu;   // NU
    real_t *x;  // x[NX]
    real_t *u;  // u[NU]

    void (*VDE_fun)(const real_t*, real_t*);
    void (*jac_fun)(const real_t*, real_t*);

#if FIXED_STEP_SIZE == 0
    real_t step;
    unsigned int nSteps;
#endif
} sim_in;

typedef struct sim_out_ {
    real_t *xn;     // xn[NX]
    real_t *Sx;     // Sx[NX*NX]
    real_t *Su;     // Su[NX*NU]
} sim_out;


#endif  // ACADOS_SIM_SIM_COMMON_H_
