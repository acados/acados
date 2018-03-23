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

typedef struct
{
    int nx;
    int nu;
    int num_stages;
    int nz;
    int nx1;
    int nx2;
    int num_steps;
    int n_out;
    int n_in;
    
} gnsf2_dims;


typedef struct {
    double *ff;
    double *x0_1;
    double *u_0;
} gnsf_res_in;


typedef struct
{
    // external functions
    external_function_generic *Phi_inc_dy;
    external_function_generic *jac_Phi_y;
    external_function_generic *f_LO_inc_J_x1k1uz;

    // precomputed matrices
    struct blasfeo_dmat KKf;
    struct blasfeo_dmat KKx;
    struct blasfeo_dmat KKu;

    struct blasfeo_dmat YYf;
    struct blasfeo_dmat YYx;
    struct blasfeo_dmat YYu;

    struct blasfeo_dmat ZZf;
    struct blasfeo_dmat ZZx;
    struct blasfeo_dmat ZZu;

    struct blasfeo_dmat ALO;
    struct blasfeo_dmat M2inv;
    struct blasfeo_dmat dK2_dx2;

    // model defining matrices
    double *A;
    double *B;
    double *C;
    double *E;
    
    double *L_x;
    double *L_xdot;
    double *L_z;
    double *L_u;
    double *A_LO;

    // butcher table maybe remove
    double* A_dt;
    double* b_dt;
    double* c;
    double dt;
} gnsf2_model;

typedef struct {
    struct blasfeo_dmat E11;
    struct blasfeo_dmat E12;
    struct blasfeo_dmat E21;
    struct blasfeo_dmat E22;

    struct blasfeo_dmat A1;
    struct blasfeo_dmat A2;
    struct blasfeo_dmat B1;
    struct blasfeo_dmat B2;
    struct blasfeo_dmat C1;
    struct blasfeo_dmat C2;

    struct blasfeo_dmat AA1;
    struct blasfeo_dmat AA2;
    struct blasfeo_dmat BB1;
    struct blasfeo_dmat BB2;
    struct blasfeo_dmat CC1;
    struct blasfeo_dmat CC2;
    struct blasfeo_dmat DD1;
    struct blasfeo_dmat DD2;
    struct blasfeo_dmat EE1;
    struct blasfeo_dmat EE2;

    struct blasfeo_dmat QQ1;
    struct blasfeo_dmat PP1;
    // struct blasfeo_dmat PP3;

    struct blasfeo_dmat LLZ;
    struct blasfeo_dmat LLx;
    struct blasfeo_dmat LLK;
    struct blasfeo_dmat LLu;

    struct blasfeo_dmat M2;
    struct blasfeo_dmat dK2_dx2_work;

    int *ipivEE1; // index of pivot vector
    int *ipivEE2; // index of pivot vector
    int *ipivQQ1; // index of pivot vector
    int *ipivPP1; // index of pivot vector
    int* ipivM2;
} gnsf2_pre_workspace;

typedef struct { //workspace
    double *phi_in;
    double *phi_out;
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

    struct blasfeo_dvec yyu;
    struct blasfeo_dvec yyss;

    struct blasfeo_dmat aux_G2_ff;
    struct blasfeo_dmat dPsi_dff;
    struct blasfeo_dmat dPsi_dx;
    struct blasfeo_dmat dPsi_du;

    struct blasfeo_dmat dPHI_dy;

    struct blasfeo_dvec *K1_val;
    struct blasfeo_dvec *x1_val;
    struct blasfeo_dvec *ff_val;
    struct blasfeo_dvec *yy_val;
    struct blasfeo_dvec *Z_val;
    struct blasfeo_dvec *f_LO_val;

    struct blasfeo_dmat *f_LO_jac;

} gnsf2_workspace;

int sim_gnsf2_model_calculate_size(void *config, sim_dims *dims);
//
void *sim_gnsf2_model_assign(void *config, sim_dims *dim_in, void *raw_memory);

void *gnsf2_cast_workspace(gnsf2_dims* dims, void *raw_memory);
int gnsf2_workspace_calculate_size(void *config, sim_dims *dim_in, void *args);

int gnsf2_pre_workspace_calculate_size(gnsf2_dims *dims);
void *gnsf2_cast_pre_workspace(gnsf2_dims* dims, void *raw_memory);

int gnsf2_dims_calculate_size();
gnsf2_dims *gnsf2_dims_assign(void *raw_memory);

int sim_gnsf2_memory_calculate_size(void *config, sim_dims *dims, void *opts_);
void *sim_gnsf2_memory_assign(void *config, sim_dims *dims, void *opts_, void *raw_memory);

void gnsf2_get_dims( gnsf2_dims* dims, casadi_function_t get_ints_fun);
void gnsf2_import_matrices(gnsf2_dims* dims, gnsf2_model *model, casadi_function_t get_matrices_fun);
void gnsf2_import_precomputed(gnsf2_dims* dims, gnsf2_model *model, casadi_function_t But_KK_YY_ZZ_LO_fun);


int sim_gnsf2_opts_calculate_size(void *config, sim_dims *dims);
void *sim_gnsf2_opts_assign(void *config, sim_dims *dims, void *raw_memory);

void gnsf2_precompute(gnsf2_dims* dims, gnsf2_model *model, sim_rk_opts *opts, sim_in *in);

void sim_gnsf2_config_initialize_default(void *config_);

void sim_gnsf2_opts_initialize_default(void *config, sim_dims *dims, void *opts_);


int gnsf2_simulate(void *config, sim_in *in, sim_out *out, void *opts, void *mem_, void *work_);
double minimum_of_doubles(double *x, int n);
void gnsf2_neville(double *out, double xx, int n, double *x, double *Q);



// double minimum_of_doubles(double *x, int n);
// void gnsf_neville(double *out, double xx, int n, double *x, double *Q);

#endif  // ACADOS_SIM_SIM_COMMON_H_
