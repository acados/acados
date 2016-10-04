#ifndef ACADOS_SIM_SIM_ERK_INTEGRATOR_H_
#define ACADOS_SIM_SIM_ERK_INTEGRATOR_H_

#include "acados/utils/types.h"

#define FIXED_STEP_SIZE 0

typedef struct sim_in_ {
    int_t nx;   // NX
    int_t nu;   // NU
    real_t *x;  // x[NX]
    real_t *u;  // u[NU]

    void (*VDE_fun)(const real_t*, real_t*);

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

typedef struct sim_RK_opts_ {
    int_t num_stages;
    real_t *A_mat;
    real_t *c_vec;
    real_t *b_vec;
} sim_RK_opts;

typedef struct sim_erk_workspace_ {
    real_t *rhs_in;
    real_t *K_tmp;
    real_t *out_tmp;
} sim_erk_workspace;


void integrate(const sim_in *in, sim_out *out, const sim_RK_opts *opts, sim_erk_workspace *work);

void sim_erk_create_workspace(const sim_in *in, sim_RK_opts *opts, sim_erk_workspace *work);

void sim_erk_create_opts(int_t num_stages, sim_RK_opts *opts);

#endif  // ACADOS_SIM_SIM_ERK_INTEGRATOR_H_
