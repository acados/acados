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
    int num_stages;
    int nx;
    int nu;
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

typedef struct {
    double *x;
    double *u;
    double *S_forw;  // forward seed
    double *S_adj;   // backward seed

} gnsf2_in;

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

typedef struct
{
	/* external functions */
	// nonlinearity functions
	external_function_generic *Phi_inc_dy;
	// Linear output function
	external_function_generic *f_LO_inc_J_x1k1uz;

    // model defining matrices
    double *A;
    double *B;
    double *C;
    double *E;
    double *L_x;
    double *L_xdot;
    double *L_z;
    double *L_u;
    double *ALO;
} gnsf2_model;

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
} gnsf2_opts;

typedef struct {
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

    double* A_dt;
    double* b_dt;
    double* c;
    double dt;

    // external functions
    external_function_generic *res_inc_Jff;
    external_function_generic *jac_res_ffx1u;
    external_function_generic *f_LO_inc_J_x1k1uz;

} gnsf2_fixed;

int gnsf2_dims_calculate_size();
gnsf2_dims *gnsf2_dims_assign(void *raw_memory);

int gnsf2_in_calculate_size(gnsf2_dims *dims);
gnsf2_in *gnsf2_in_assign(gnsf2_dims *dims, void *raw_memory);

void gnsf2_get_dims( gnsf2_dims* dims, casadi_function_t get_ints_fun);

int gnsf2_opts_calculate_size(gnsf2_dims *dims);
gnsf2_opts *gnsf2_opts_assign(gnsf2_dims *dims, void *raw_memory);

int gnsf2_fixed_calculate_size(gnsf2_dims *dims, gnsf2_opts* opts);
gnsf2_fixed *gnsf2_fixed_assign(gnsf2_dims *dims, void *raw_memory, int memsize);
void gnsf2_import(gnsf2_dims* dims, gnsf2_fixed *fix, casadi_function_t But_KK_YY_ZZ_LO_fun);


#endif  // ACADOS_SIM_SIM_COMMON_H_
