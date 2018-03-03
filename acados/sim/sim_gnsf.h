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

#include "acados/utils/timing.h"
#include "acados/utils/types.h"

#include "acados/sim/sim_common.h"

#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_blas.h"
#include "blasfeo/include/blasfeo_d_kernel.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_aux.h"


typedef struct {
    int num_stages;
    int nx;
    int nu;
    int nz;
    int nx1;
    int nx2;
    int num_steps;
    int n_out;
} gnsf_dims;

typedef struct {
    double *ff;
    double *x0_1;
    double *u_0;
} gnsf_res_in;


typedef struct {
    double *x;
    double *u;
    double *S_forw;  // forward seed
    double *S_adj;   // backward seed

} gnsf_in;

typedef struct {
    struct blasfeo_dmat KKf;
    struct blasfeo_dmat KKx;
    struct blasfeo_dmat KKu;

    struct blasfeo_dmat ZZf;
    struct blasfeo_dmat ZZx;
    struct blasfeo_dmat ZZu;

    struct blasfeo_dmat ALO;
    struct blasfeo_dmat M2inv;
    struct blasfeo_dmat dK2_dx2;

    double* A_dt;
    double* b_dt;
    double* c;
    double dt;

    // external functions
    external_function_generic *res_inc_Jff;
    external_function_generic *jac_res_ffx1u;
    external_function_generic *f_LO_inc_J_x1k1uz;


} gnsf_fixed;

typedef struct { //workspace
    double *res_in;
    double *res_out;
    double *f_LO_in;
    double *f_LO_out;
    double *Z_out;

    int *ipiv; // index of pivot vector

    struct blasfeo_dmat J_r_ff;
    struct blasfeo_dmat J_r_x1u; 
    
    struct blasfeo_dmat dK1_dx1;
    struct blasfeo_dmat dK1_du;
    struct blasfeo_dmat dZ_dx1;
    struct blasfeo_dmat dZ_du;
    struct blasfeo_dmat aux_G2_x1;
    struct blasfeo_dmat aux_G2_u;
    struct blasfeo_dmat J_G2_K1;
    
    struct blasfeo_dmat dK2_dx1;
    struct blasfeo_dmat dK2_du;
    struct blasfeo_dmat dK2_dff;
    struct blasfeo_dmat dxf_dwn;
    struct blasfeo_dmat S_forw_new;
    struct blasfeo_dmat S_forw;

    struct blasfeo_dvec K2_val;
    struct blasfeo_dvec x0_traj;
    struct blasfeo_dvec res_val;
    struct blasfeo_dvec u0;
    struct blasfeo_dvec lambda;
    struct blasfeo_dvec lambda_old;

    struct blasfeo_dmat aux_G2_ff;
    struct blasfeo_dmat dPsi_dff;
    struct blasfeo_dmat dPsi_dx;
    struct blasfeo_dmat dPsi_du;

    struct blasfeo_dvec *K1_val;
    struct blasfeo_dvec *x1_val;
    struct blasfeo_dvec *ff_val;
    struct blasfeo_dvec *Z_val;
    struct blasfeo_dvec *f_LO_val;

    struct blasfeo_dmat *f_LO_jac;

} gnsf_workspace;
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
    int  newton_max;
    bool jac_reuse;
   // Newton_scheme *scheme;
} gnsf_opts;

void print_gnsf_dims(gnsf_dims* dims);
void print_gnsf_res_in( gnsf_dims* dims, double* res_in );
void print_gnsf_res_out( gnsf_dims* dims, double* res_out );
void gnsf_get_dims( gnsf_dims* dims, casadi_function_t get_ints_fun);
void gnsf_import(gnsf_dims* dims, gnsf_fixed *fix, casadi_function_t But_KK_ZZ_LO_fun);

void *gnsf_cast_workspace(gnsf_dims* dims, void *raw_memory);
int gnsf_calculate_workspace_size(gnsf_dims *dims, gnsf_opts* opts);

int gnsf_fixed_calculate_size(gnsf_dims *dims, gnsf_opts* opts);
void gnsf_allocate_fixed( gnsf_dims *dims, gnsf_fixed *fix);
gnsf_fixed *gnsf_fixed_assign(gnsf_dims *dims, void *raw_memory, int memsize);

void gnsf_simulate( gnsf_dims* dims, gnsf_fixed* fix, gnsf_in* in, sim_out* out, gnsf_opts* opts, void *work_);

void gnsf_neville(double *out, double xx, int n, double *x, double *Q);
double minimum_of_doubles(double *x, int n);

int gnsf_dims_calculate_size();
gnsf_dims *gnsf_dims_assign(void *raw_memory);

int gnsf_in_calculate_size(gnsf_dims *dims);
gnsf_in *gnsf_in_assign(gnsf_dims *dims, void *raw_memory);

int gnsf_opts_calculate_size(gnsf_dims *dims);
gnsf_opts *gnsf_opts_assign(gnsf_dims *dims, void *raw_memory);

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
