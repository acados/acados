/*
 *    This file is part of acados.
 *
 *    acados is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    acados is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with acados; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#ifndef ACADOS_SIM_SIM_COMMON_H_
#define ACADOS_SIM_SIM_COMMON_H_

#include <stdbool.h>

#include "acados/sim/sim_collocation_utils.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"

#include "acados/utils/external_function_generic.h"

// maximum number of integration stages
#define NS_MAX 15



typedef enum
{
    // ERK and LIFTED_ERK
    EXPL_ODE_FUN,
    EXPL_ODE_HES,  // wrt x and u ???
    EXPL_VDE_FOR,
    EXPL_VDE_ADJ,
    // IRK
    IMPL_ODE_FUN,
    IMPL_ODE_FUN_JAC_X_XDOT,
    IMPL_ODE_JAC_X_XDOT_U,
    IMPL_ODE_FUN_JAC_X_XDOT_U,
    IMPL_ODE_HESS,
    // gnsf
    PHI_FUN,
    PHI_FUN_JAC_Y,
    PHI_JAC_Y_UHAT,
    LO_FUN,
    GET_GNSF_MATRICES
} sim_function_t;



typedef struct
{
    void *dims;

    double *x;  // x[NX] - initial state value for simulation
    double *u;  // u[NU] - control - constant over simulation time

    double *xdot;  // xdot[NX] - initialization for state derivatives k within the integrator
    double *z;     // z[NZ] - initialization for algebraic variables z

    double *S_forw;  // forward seed [Sx, Su]
    double *S_adj;   // backward seed

    void *model;

    double T;  // simulation time

} sim_in;



typedef struct
{
    double CPUtime;  // in seconds
    double LAtime;   // in seconds
    double ADtime;   // in seconds
} sim_info;



typedef struct
{
    double *xn;      // xn[NX]
    double *S_forw;  // S_forw[NX*(NX+NU)]
    double *S_adj;   //
    double *S_hess;  //

    double *zn;           // z - algebraic variables - reported at start of simulation interval
    double *S_algebraic;  // sensitivities of reported value of algebraic variables w.r.t.
                          // initial stat & control (x_n,u)

    double *grad;  // gradient correction

    sim_info *info;
} sim_out;



typedef struct
{
    int ns;  // number of integration stages

    int num_steps;
    int num_forw_sens;

    int tableau_size;  // check that is consistent with ns
            // only update when butcher tableau is changed
            // kind of private -> no setter!
    double *A_mat;
    double *c_vec;
    double *b_vec;

    bool sens_forw;
    bool sens_adj;
    bool sens_hess;

    bool output_z;        // 1 -- if zn should be computed
    bool sens_algebraic;  // 1 -- if S_algebraic should be computed

    // for explicit integrators: newton_iter == 0 && scheme == NULL
    // && jac_reuse=false
    int newton_iter;
    bool jac_reuse;
    Newton_scheme *scheme;

    // workspace
    void *work;

} sim_opts;



typedef struct
{
    int (*evaluate)(void *config_, sim_in *in, sim_out *out, void *opts, void *mem, void *work);
    int (*precompute)(void *config_, sim_in *in, sim_out *out, void *opts, void *mem, void *work);
    int (*opts_calculate_size)(void *config_, void *dims);
    void *(*opts_assign)(void *config_, void *dims, void *raw_memory);
    void (*opts_initialize_default)(void *config_, void *dims, void *opts);
    void (*opts_update)(void *config_, void *dims, void *opts);
    int (*opts_set)(void *config_, void *opts_, const char *field, void *value);
    int (*memory_calculate_size)(void *config, void *dims, void *opts);
    void *(*memory_assign)(void *config, void *dims, void *opts, void *raw_memory);
    int (*workspace_calculate_size)(void *config, void *dims, void *opts);
    int (*model_calculate_size)(void *config, void *dims);
    void *(*model_assign)(void *config, void *dims, void *raw_memory);
    int (*model_set)(void *model, const char *field, void *value);
    void (*config_initialize_default)(void *config);
//    int (*dims_calculate_size)(void *config);
    int (*dims_calculate_size)();
    void *(*dims_assign)(void *config, void *raw_memory);
    void (*dims_set)(void *config, void *dims, const char *field, const int *value);
    void (*dims_get)(void *config, void *dims, const char *field, int *value);

} sim_config;



/* config */
//
int sim_config_calculate_size();
//
sim_config *sim_config_assign(void *raw_memory);

/* in */
//
int sim_in_calculate_size(void *config, void *dims);
//
sim_in *sim_in_assign(void *config, void *dims, void *raw_memory);
//
int sim_in_set_(void *config_, void *dims_, sim_in *in, const char *field, void *value);

/* out */
//
int sim_out_calculate_size(void *config, void *dims);
//
sim_out *sim_out_assign(void *config, void *dims, void *raw_memory);
//
int sim_out_get_(void *config, void *dims, sim_out *out, const char *field, void *value);

/* opts */
//
int sim_opts_set_(sim_opts *opts, const char *field, void *value);


#endif  // ACADOS_SIM_SIM_COMMON_H_
