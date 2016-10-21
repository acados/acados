#include <stdlib.h>
#include <string.h>

#include "hpmpc/include/aux_d.h"
#include "acados/sim/sim_erk_integrator.h"
#include "acados/utils/print.h"

#if FIXED_STEP_SIZE == 1
// Integration step size.
real_t H_INT = 1.0/100;
// Number of integration steps.
#define NSTEPS  10
#endif

void sim_erk(const sim_in *in, sim_out *out, const sim_RK_opts *opts, sim_erk_workspace *work) {
    int_t nx = in->nx;
    int_t nu = in->nu;
    int_t num_stages = opts->num_stages;
    int_t i, s, j, istep;
#if FIXED_STEP_SIZE == 0
    real_t H_INT = in->step;
    int_t NSTEPS = in->nSteps;
#endif

    real_t *A_mat = opts->A_mat;
    real_t *b_vec = opts->b_vec;
//    real_t *c_vec = opts->c_vec;

//    print_matrix("stdout", A_mat, num_stages, num_stages);

    real_t *K_tmp = work->K_tmp;
    real_t *out_tmp = work->out_tmp;
    real_t *rhs_in = work->rhs_in;

    for (i = 0; i < nx; i++) out_tmp[i] = in->x[i];
    for (i = 0; i < nx*(nx+nu); i++) out_tmp[nx+i] = 0.0;  // sensitivities
    for (i = 0; i < nx; i++) out_tmp[nx+i*nx+i] = 1.0;     // sensitivities wrt x

    for (i = 0; i < nu; i++) rhs_in[nx*(1+nx+nu)+i] = in->u[i];

    for (istep = 0; istep < NSTEPS; istep++) {
        for (s = 0; s < num_stages; s++) {
            for (i = 0; i < nx*(1+nx+nu); i++) {
                rhs_in[i] = out_tmp[i];
            }
            for (j = 0; j < s; j++) {
                if (A_mat[j*num_stages+s] != 0) {
                    for (i = 0; i < nx*(1+nx+nu); i++) {
                        rhs_in[i] += H_INT*A_mat[j*num_stages+s]*K_tmp[j*nx*(1+nx+nu)+i];
                    }
                }
            }
            in->VDE_fun(rhs_in, &(K_tmp[s*nx*(1+nx+nu)]));  // k evaluation
        }
        for (s = 0; s < num_stages; s++) {
            for (i = 0; i < nx*(1+nx+nu); i++) {
                out_tmp[i] += H_INT*b_vec[s]*K_tmp[s*nx*(1+nx+nu)+i];  // ERK step
            }
        }
    }
    for (i = 0; i < nx; i++)      out->xn[i] = out_tmp[i];
    for (i = 0; i < nx*nx; i++) out->Sx[i] = out_tmp[nx+i];
    for (i = 0; i < nx*nu; i++) out->Su[i] = out_tmp[nx+nx*nx+i];
}


void sim_erk_create_workspace(const sim_in *in, sim_RK_opts *opts, sim_erk_workspace *work) {
    int_t nx = in->nx;
    int_t nu = in->nu;
    int_t num_stages = opts->num_stages;

    work->rhs_in = malloc(sizeof(*work->rhs_in) * (nx*(1+nx+nu)+nu));
    work->K_tmp = malloc(sizeof(*work->K_tmp) * (num_stages*nx*(1+nx+nu)));
    work->out_tmp = malloc(sizeof(*work->out_tmp) * (nx*(1+nx+nu)));
}


void sim_erk_create_opts(const int_t num_stages, sim_RK_opts *opts) {
    if ( num_stages == 1 ) {
        opts->num_stages = 1;       // explicit Euler
        opts->A_mat = malloc(sizeof(*opts->A_mat) * (num_stages*num_stages));
        opts->b_vec = malloc(sizeof(*opts->b_vec) * (num_stages));
        opts->c_vec = malloc(sizeof(*opts->c_vec) * (num_stages));
        opts->A_mat[0] = 0;
        opts->b_vec[0] = 1.0;
        opts->c_vec[0] = 0;
    } else if ( num_stages == 4 ) {
        opts->num_stages = 4;       // 4-stage, 4th order ERK method
        opts->A_mat = malloc(sizeof(*opts->A_mat) * (num_stages*num_stages));
        opts->b_vec = malloc(sizeof(*opts->b_vec) * (num_stages));
        opts->c_vec = malloc(sizeof(*opts->c_vec) * (num_stages));

        memcpy(opts->A_mat,
                ((real_t[]) {0, 0.5, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 1, 0, 0, 0, 0}),
                sizeof(*opts->A_mat) * (num_stages*num_stages));
        memcpy(opts->b_vec,
                ((real_t[]) {1.0/6, 2.0/6, 2.0/6, 1.0/6}), sizeof(*opts->b_vec) * (num_stages));
        memcpy(opts->c_vec,
                ((real_t[]) {0.0, 0.5, 0.5, 1.0}), sizeof(*opts->c_vec) * (num_stages));
    } else {
        // throw error somehow?
    }
}
