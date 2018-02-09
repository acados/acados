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

#ifndef ACADOS_SIM_SIM_GNSF_H_
#define ACADOS_SIM_SIM_GNSF_H_

#include <stdbool.h>

#include "acados/sim/sim_collocation.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"


typedef struct {
    int num_stages;
    int nx;
    int nu;
    int nz;
    int nx1;
    int nx2;
    int num_steps;
    int n_out;
    int nff;
} gnsf_dims;

typedef struct {
    double *ff;
    double *x0_1;
    double *u_0;
} gnsf_res_in;

// TODO

typedef struct {
    double *xf; //TODO
} gnsf_out;

typedef struct {
    double *A_dt;
    double *b_dt;
    double *c_butcher;

    double *x_0;
    double *u_0;

    double *KKf; //etc...

    double *S_forw;  // forward seed
    double *S_adj;   // backward seed

        // casadi_functions
    casadi_function_t jac_res_ffx1u;
    void (*jac_res_ffx1u_wrapped)(const int, const int, const double *, double *, casadi_function_t);

    casadi_function_t res_inc_Jff;
    void (*res_inc_Jff_wrapped)(const int, const int, const double *, double *, casadi_function_t);

} gnsf_in;

typedef struct {
    double *KKf;
    double *KKx;
    double *KKu;
} gnsf_fixed;
/*
typedef struct {
    double CPUtime;
    double LAtime;
    double ADtime;
} sim_info; */


/*typedef struct {
    double *xn;      // xn[NX]
    double *S_forw;  // S_forw[NX*(NX+NU)]
    double *S_adj;   //
  //  double *S_hess;  //

//    double *grad;  // gradient correction

//    sim_info *info;
} sim_out; */


typedef struct {
//    double interval;
//    int num_stages;
    int num_steps;
//    int num_forw_sens;
    bool sens_forw;
    bool sens_adj;
    bool sens_hess;
    int newton_iter;
   // Newton_scheme *scheme;
} sim_gnsf_opts;

void print_gnsf_dims(gnsf_dims* dims);
void print_gnsf_res_in( gnsf_dims dims, double* res_in );
void print_gnsf_res_out( gnsf_dims dims, double* res_out );
void gnsf_get_dims( gnsf_dims* dims, casadi_function_t get_ints_fun);
// void gnsf_get_KK_mat( gnsf_in *in, casadi_function_t KK_mat_fun);
void gnsf_get_KK_mat(gnsf_dims* dims, gnsf_fixed *fix, casadi_function_t KK_mat_fun);
void gnsf_simulate( gnsf_dims* dims, gnsf_in* in, gnsf_out out);
void gnsf_allocate_fixed( gnsf_dims *dims, gnsf_fixed *fix);
/*
typedef struct {
    int (*fun)(sim_in *in, sim_out *out, void *args, void *mem, void *work);
    int (*calculate_args_size)(sim_dims *dims);
    void *(*assign_args)(sim_dims *dims, void *raw_memory);
    void (*initialize_default_args)(sim_dims *dims, void *args);
    int (*calculate_memory_size)(sim_dims *dims, void *args);
    void *(*assign_memory)(sim_dims *dims, void *args, void *raw_memory);
    int (*calculate_workspace_size)(sim_dims *dims, void *args);
} sim_solver_fcn_ptrs;


int sim_dims_calculate_size();

sim_dims *assign_sim_dims(void *raw_memory);

int sim_in_calculate_size(sim_dims *dims);

sim_in *assign_sim_in(sim_dims *dims, void *raw_memory);

int sim_out_calculate_size(sim_dims *dims);

sim_out *assign_sim_out(sim_dims *dims, void *raw_memory);
*/
#endif  // ACADOS_SIM_SIM_COMMON_H_
