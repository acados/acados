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

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// #include <xmmintrin.h>

#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"

#include "acados_c/external_function_interface.h"
#include "acados_c/ocp_nlp_interface.h"

// TODO(dimitris): use only the strictly necessary includes here

#include "acados/utils/mem.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"

#include "acados/ocp_nlp/ocp_nlp_sqp.h"
#include "acados/ocp_nlp/ocp_nlp_cost_common.h"
#include "acados/ocp_nlp/ocp_nlp_cost_ls.h"
#include "acados/ocp_nlp/ocp_nlp_cost_nls.h"
#include "acados/ocp_nlp/ocp_nlp_cost_external.h"
#include "acados/ocp_nlp/ocp_nlp_dynamics_cont.h"

#include "examples/c/chain_model/chain_model.h"
#include "examples/c/implicit_chain_model/chain_model_impl.h"

// x0
#include "examples/c/chain_model/x0_nm2.c"
#include "examples/c/chain_model/x0_nm3.c"
#include "examples/c/chain_model/x0_nm4.c"
#include "examples/c/chain_model/x0_nm5.c"
#include "examples/c/chain_model/x0_nm6.c"

// xN
#include "examples/c/chain_model/xN_nm2.c"
#include "examples/c/chain_model/xN_nm3.c"
#include "examples/c/chain_model/xN_nm4.c"
#include "examples/c/chain_model/xN_nm5.c"
#include "examples/c/chain_model/xN_nm6.c"

#define NN 15
#define TF 3.75
#define MAX_SQP_ITERS 10
#define NREP 1


static void shift_states(ocp_nlp_dims *dims, ocp_nlp_out *out, double *x_end)
{
	int N = dims->N;

    for (int i = 0; i < N; i++)
 		blasfeo_dveccp(dims->nx[i], &out->ux[i], dims->nu[i], &out->ux[i+1], dims->nu[i+1]);
 	blasfeo_pack_dvec(dims->nx[N], x_end, &out->ux[N], dims->nu[N]);
}



static void shift_controls(ocp_nlp_dims *dims, ocp_nlp_out *out, double *u_end)
{
	int N = dims->N;

    for (int i = 0; i < N-1; i++)
 		blasfeo_dveccp(dims->nu[i], &out->ux[i], 0, &out->ux[i+1], 0);
 	blasfeo_pack_dvec(dims->nu[N-1], u_end, &out->ux[N-1], 0);
}



static void select_dynamics_casadi(int N, int num_free_masses, external_function_casadi *forw_vde, external_function_casadi *jac_ode, external_function_casadi *impl_ode, external_function_casadi *impl_jac_x, external_function_casadi *impl_jac_xdot, external_function_casadi *impl_jac_u)
{
	switch (num_free_masses)
	{
		case 1:
			for (int ii = 0; ii < N; ii++)
			{
				forw_vde[ii].casadi_fun = &vde_chain_nm2;
				forw_vde[ii].casadi_work = &vde_chain_nm2_work;
				forw_vde[ii].casadi_sparsity_in = &vde_chain_nm2_sparsity_in;
				forw_vde[ii].casadi_sparsity_out = &vde_chain_nm2_sparsity_out;
				forw_vde[ii].casadi_n_in = &vde_chain_nm2_n_in;
				forw_vde[ii].casadi_n_out = &vde_chain_nm2_n_out;
				jac_ode[ii].casadi_fun = &jac_chain_nm2;
				jac_ode[ii].casadi_work = &jac_chain_nm2_work;
				jac_ode[ii].casadi_sparsity_in = &jac_chain_nm2_sparsity_in;
				jac_ode[ii].casadi_sparsity_out = &jac_chain_nm2_sparsity_out;
				jac_ode[ii].casadi_n_in = &jac_chain_nm2_n_in;
				jac_ode[ii].casadi_n_out = &jac_chain_nm2_n_out;

				impl_ode[ii].casadi_fun = &impl_odeFun_chain_nm2;
				impl_ode[ii].casadi_work = &impl_odeFun_chain_nm2_work;
				impl_ode[ii].casadi_sparsity_in = &impl_odeFun_chain_nm2_sparsity_in;
				impl_ode[ii].casadi_sparsity_out = &impl_odeFun_chain_nm2_sparsity_out;
				impl_ode[ii].casadi_n_in = &impl_odeFun_chain_nm2_n_in;
				impl_ode[ii].casadi_n_out = &impl_odeFun_chain_nm2_n_out;
				impl_jac_x[ii].casadi_fun = &impl_jacFun_x_chain_nm2;
				impl_jac_x[ii].casadi_work = &impl_jacFun_x_chain_nm2_work;
				impl_jac_x[ii].casadi_sparsity_in = &impl_jacFun_x_chain_nm2_sparsity_in;
				impl_jac_x[ii].casadi_sparsity_out = &impl_jacFun_x_chain_nm2_sparsity_out;
				impl_jac_x[ii].casadi_n_in = &impl_jacFun_x_chain_nm2_n_in;
				impl_jac_x[ii].casadi_n_out = &impl_jacFun_x_chain_nm2_n_out;
				impl_jac_xdot[ii].casadi_fun = &impl_jacFun_xdot_chain_nm2;
				impl_jac_xdot[ii].casadi_work = &impl_jacFun_xdot_chain_nm2_work;
				impl_jac_xdot[ii].casadi_sparsity_in = &impl_jacFun_xdot_chain_nm2_sparsity_in;
				impl_jac_xdot[ii].casadi_sparsity_out = &impl_jacFun_xdot_chain_nm2_sparsity_out;
				impl_jac_xdot[ii].casadi_n_in = &impl_jacFun_xdot_chain_nm2_n_in;
				impl_jac_xdot[ii].casadi_n_out = &impl_jacFun_xdot_chain_nm2_n_out;
				impl_jac_u[ii].casadi_fun = &impl_jacFun_u_chain_nm2;
				impl_jac_u[ii].casadi_work = &impl_jacFun_u_chain_nm2_work;
				impl_jac_u[ii].casadi_sparsity_in = &impl_jacFun_u_chain_nm2_sparsity_in;
				impl_jac_u[ii].casadi_sparsity_out = &impl_jacFun_u_chain_nm2_sparsity_out;
				impl_jac_u[ii].casadi_n_in = &impl_jacFun_u_chain_nm2_n_in;
				impl_jac_u[ii].casadi_n_out = &impl_jacFun_u_chain_nm2_n_out;
			}
			break;
		case 2:
			for (int ii = 0; ii < N; ii++)
			{
				forw_vde[ii].casadi_fun = &vde_chain_nm3;
				forw_vde[ii].casadi_work = &vde_chain_nm3_work;
				forw_vde[ii].casadi_sparsity_in = &vde_chain_nm3_sparsity_in;
				forw_vde[ii].casadi_sparsity_out = &vde_chain_nm3_sparsity_out;
				forw_vde[ii].casadi_n_in = &vde_chain_nm3_n_in;
				forw_vde[ii].casadi_n_out = &vde_chain_nm3_n_out;
				jac_ode[ii].casadi_fun = &jac_chain_nm3;
				jac_ode[ii].casadi_work = &jac_chain_nm3_work;
				jac_ode[ii].casadi_sparsity_in = &jac_chain_nm3_sparsity_in;
				jac_ode[ii].casadi_sparsity_out = &jac_chain_nm3_sparsity_out;
				jac_ode[ii].casadi_n_in = &jac_chain_nm3_n_in;
				jac_ode[ii].casadi_n_out = &jac_chain_nm3_n_out;

				impl_ode[ii].casadi_fun = &impl_odeFun_chain_nm3;
				impl_ode[ii].casadi_work = &impl_odeFun_chain_nm3_work;
				impl_ode[ii].casadi_sparsity_in = &impl_odeFun_chain_nm3_sparsity_in;
				impl_ode[ii].casadi_sparsity_out = &impl_odeFun_chain_nm3_sparsity_out;
				impl_ode[ii].casadi_n_in = &impl_odeFun_chain_nm3_n_in;
				impl_ode[ii].casadi_n_out = &impl_odeFun_chain_nm3_n_out;
				impl_jac_x[ii].casadi_fun = &impl_jacFun_x_chain_nm3;
				impl_jac_x[ii].casadi_work = &impl_jacFun_x_chain_nm3_work;
				impl_jac_x[ii].casadi_sparsity_in = &impl_jacFun_x_chain_nm3_sparsity_in;
				impl_jac_x[ii].casadi_sparsity_out = &impl_jacFun_x_chain_nm3_sparsity_out;
				impl_jac_x[ii].casadi_n_in = &impl_jacFun_x_chain_nm3_n_in;
				impl_jac_x[ii].casadi_n_out = &impl_jacFun_x_chain_nm3_n_out;
				impl_jac_xdot[ii].casadi_fun = &impl_jacFun_xdot_chain_nm3;
				impl_jac_xdot[ii].casadi_work = &impl_jacFun_xdot_chain_nm3_work;
				impl_jac_xdot[ii].casadi_sparsity_in = &impl_jacFun_xdot_chain_nm3_sparsity_in;
				impl_jac_xdot[ii].casadi_sparsity_out = &impl_jacFun_xdot_chain_nm3_sparsity_out;
				impl_jac_xdot[ii].casadi_n_in = &impl_jacFun_xdot_chain_nm3_n_in;
				impl_jac_xdot[ii].casadi_n_out = &impl_jacFun_xdot_chain_nm3_n_out;
				impl_jac_u[ii].casadi_fun = &impl_jacFun_u_chain_nm3;
				impl_jac_u[ii].casadi_work = &impl_jacFun_u_chain_nm3_work;
				impl_jac_u[ii].casadi_sparsity_in = &impl_jacFun_u_chain_nm3_sparsity_in;
				impl_jac_u[ii].casadi_sparsity_out = &impl_jacFun_u_chain_nm3_sparsity_out;
				impl_jac_u[ii].casadi_n_in = &impl_jacFun_u_chain_nm3_n_in;
				impl_jac_u[ii].casadi_n_out = &impl_jacFun_u_chain_nm3_n_out;
			}
			break;
		case 3:
			for (int ii = 0; ii < N; ii++)
			{
				forw_vde[ii].casadi_fun = &vde_chain_nm4;
				forw_vde[ii].casadi_work = &vde_chain_nm4_work;
				forw_vde[ii].casadi_sparsity_in = &vde_chain_nm4_sparsity_in;
				forw_vde[ii].casadi_sparsity_out = &vde_chain_nm4_sparsity_out;
				forw_vde[ii].casadi_n_in = &vde_chain_nm4_n_in;
				forw_vde[ii].casadi_n_out = &vde_chain_nm4_n_out;
				jac_ode[ii].casadi_fun = &jac_chain_nm4;
				jac_ode[ii].casadi_work = &jac_chain_nm4_work;
				jac_ode[ii].casadi_sparsity_in = &jac_chain_nm4_sparsity_in;
				jac_ode[ii].casadi_sparsity_out = &jac_chain_nm4_sparsity_out;
				jac_ode[ii].casadi_n_in = &jac_chain_nm4_n_in;
				jac_ode[ii].casadi_n_out = &jac_chain_nm4_n_out;

				impl_ode[ii].casadi_fun = &impl_odeFun_chain_nm4;
				impl_ode[ii].casadi_work = &impl_odeFun_chain_nm4_work;
				impl_ode[ii].casadi_sparsity_in = &impl_odeFun_chain_nm4_sparsity_in;
				impl_ode[ii].casadi_sparsity_out = &impl_odeFun_chain_nm4_sparsity_out;
				impl_ode[ii].casadi_n_in = &impl_odeFun_chain_nm4_n_in;
				impl_ode[ii].casadi_n_out = &impl_odeFun_chain_nm4_n_out;
				impl_jac_x[ii].casadi_fun = &impl_jacFun_x_chain_nm4;
				impl_jac_x[ii].casadi_work = &impl_jacFun_x_chain_nm4_work;
				impl_jac_x[ii].casadi_sparsity_in = &impl_jacFun_x_chain_nm4_sparsity_in;
				impl_jac_x[ii].casadi_sparsity_out = &impl_jacFun_x_chain_nm4_sparsity_out;
				impl_jac_x[ii].casadi_n_in = &impl_jacFun_x_chain_nm4_n_in;
				impl_jac_x[ii].casadi_n_out = &impl_jacFun_x_chain_nm4_n_out;
				impl_jac_xdot[ii].casadi_fun = &impl_jacFun_xdot_chain_nm4;
				impl_jac_xdot[ii].casadi_work = &impl_jacFun_xdot_chain_nm4_work;
				impl_jac_xdot[ii].casadi_sparsity_in = &impl_jacFun_xdot_chain_nm4_sparsity_in;
				impl_jac_xdot[ii].casadi_sparsity_out = &impl_jacFun_xdot_chain_nm4_sparsity_out;
				impl_jac_xdot[ii].casadi_n_in = &impl_jacFun_xdot_chain_nm4_n_in;
				impl_jac_xdot[ii].casadi_n_out = &impl_jacFun_xdot_chain_nm4_n_out;
				impl_jac_u[ii].casadi_fun = &impl_jacFun_u_chain_nm4;
				impl_jac_u[ii].casadi_work = &impl_jacFun_u_chain_nm4_work;
				impl_jac_u[ii].casadi_sparsity_in = &impl_jacFun_u_chain_nm4_sparsity_in;
				impl_jac_u[ii].casadi_sparsity_out = &impl_jacFun_u_chain_nm4_sparsity_out;
				impl_jac_u[ii].casadi_n_in = &impl_jacFun_u_chain_nm4_n_in;
				impl_jac_u[ii].casadi_n_out = &impl_jacFun_u_chain_nm4_n_out;
			}
			break;
		case 4:
			for (int ii = 0; ii < N; ii++)
			{
				forw_vde[ii].casadi_fun = &vde_chain_nm5;
				forw_vde[ii].casadi_work = &vde_chain_nm5_work;
				forw_vde[ii].casadi_sparsity_in = &vde_chain_nm5_sparsity_in;
				forw_vde[ii].casadi_sparsity_out = &vde_chain_nm5_sparsity_out;
				forw_vde[ii].casadi_n_in = &vde_chain_nm5_n_in;
				forw_vde[ii].casadi_n_out = &vde_chain_nm5_n_out;
				jac_ode[ii].casadi_fun = &jac_chain_nm5;
				jac_ode[ii].casadi_work = &jac_chain_nm5_work;
				jac_ode[ii].casadi_sparsity_in = &jac_chain_nm5_sparsity_in;
				jac_ode[ii].casadi_sparsity_out = &jac_chain_nm5_sparsity_out;
				jac_ode[ii].casadi_n_in = &jac_chain_nm5_n_in;
				jac_ode[ii].casadi_n_out = &jac_chain_nm5_n_out;

				impl_ode[ii].casadi_fun = &impl_odeFun_chain_nm5;
				impl_ode[ii].casadi_work = &impl_odeFun_chain_nm5_work;
				impl_ode[ii].casadi_sparsity_in = &impl_odeFun_chain_nm5_sparsity_in;
				impl_ode[ii].casadi_sparsity_out = &impl_odeFun_chain_nm5_sparsity_out;
				impl_ode[ii].casadi_n_in = &impl_odeFun_chain_nm5_n_in;
				impl_ode[ii].casadi_n_out = &impl_odeFun_chain_nm5_n_out;
				impl_jac_x[ii].casadi_fun = &impl_jacFun_x_chain_nm5;
				impl_jac_x[ii].casadi_work = &impl_jacFun_x_chain_nm5_work;
				impl_jac_x[ii].casadi_sparsity_in = &impl_jacFun_x_chain_nm5_sparsity_in;
				impl_jac_x[ii].casadi_sparsity_out = &impl_jacFun_x_chain_nm5_sparsity_out;
				impl_jac_x[ii].casadi_n_in = &impl_jacFun_x_chain_nm5_n_in;
				impl_jac_x[ii].casadi_n_out = &impl_jacFun_x_chain_nm5_n_out;
				impl_jac_xdot[ii].casadi_fun = &impl_jacFun_xdot_chain_nm5;
				impl_jac_xdot[ii].casadi_work = &impl_jacFun_xdot_chain_nm5_work;
				impl_jac_xdot[ii].casadi_sparsity_in = &impl_jacFun_xdot_chain_nm5_sparsity_in;
				impl_jac_xdot[ii].casadi_sparsity_out = &impl_jacFun_xdot_chain_nm5_sparsity_out;
				impl_jac_xdot[ii].casadi_n_in = &impl_jacFun_xdot_chain_nm5_n_in;
				impl_jac_xdot[ii].casadi_n_out = &impl_jacFun_xdot_chain_nm5_n_out;
				impl_jac_u[ii].casadi_fun = &impl_jacFun_u_chain_nm5;
				impl_jac_u[ii].casadi_work = &impl_jacFun_u_chain_nm5_work;
				impl_jac_u[ii].casadi_sparsity_in = &impl_jacFun_u_chain_nm5_sparsity_in;
				impl_jac_u[ii].casadi_sparsity_out = &impl_jacFun_u_chain_nm5_sparsity_out;
				impl_jac_u[ii].casadi_n_in = &impl_jacFun_u_chain_nm5_n_in;
				impl_jac_u[ii].casadi_n_out = &impl_jacFun_u_chain_nm5_n_out;
			}
			break;
		case 5:
			for (int ii = 0; ii < N; ii++)
			{
				forw_vde[ii].casadi_fun = &vde_chain_nm6;
				forw_vde[ii].casadi_work = &vde_chain_nm6_work;
				forw_vde[ii].casadi_sparsity_in = &vde_chain_nm6_sparsity_in;
				forw_vde[ii].casadi_sparsity_out = &vde_chain_nm6_sparsity_out;
				forw_vde[ii].casadi_n_in = &vde_chain_nm6_n_in;
				forw_vde[ii].casadi_n_out = &vde_chain_nm6_n_out;
				jac_ode[ii].casadi_fun = &jac_chain_nm6;
				jac_ode[ii].casadi_work = &jac_chain_nm6_work;
				jac_ode[ii].casadi_sparsity_in = &jac_chain_nm6_sparsity_in;
				jac_ode[ii].casadi_sparsity_out = &jac_chain_nm6_sparsity_out;
				jac_ode[ii].casadi_n_in = &jac_chain_nm6_n_in;
				jac_ode[ii].casadi_n_out = &jac_chain_nm6_n_out;

				impl_ode[ii].casadi_fun = &impl_odeFun_chain_nm6;
				impl_ode[ii].casadi_work = &impl_odeFun_chain_nm6_work;
				impl_ode[ii].casadi_sparsity_in = &impl_odeFun_chain_nm6_sparsity_in;
				impl_ode[ii].casadi_sparsity_out = &impl_odeFun_chain_nm6_sparsity_out;
				impl_ode[ii].casadi_n_in = &impl_odeFun_chain_nm6_n_in;
				impl_ode[ii].casadi_n_out = &impl_odeFun_chain_nm6_n_out;
				impl_jac_x[ii].casadi_fun = &impl_jacFun_x_chain_nm6;
				impl_jac_x[ii].casadi_work = &impl_jacFun_x_chain_nm6_work;
				impl_jac_x[ii].casadi_sparsity_in = &impl_jacFun_x_chain_nm6_sparsity_in;
				impl_jac_x[ii].casadi_sparsity_out = &impl_jacFun_x_chain_nm6_sparsity_out;
				impl_jac_x[ii].casadi_n_in = &impl_jacFun_x_chain_nm6_n_in;
				impl_jac_x[ii].casadi_n_out = &impl_jacFun_x_chain_nm6_n_out;
				impl_jac_xdot[ii].casadi_fun = &impl_jacFun_xdot_chain_nm6;
				impl_jac_xdot[ii].casadi_work = &impl_jacFun_xdot_chain_nm6_work;
				impl_jac_xdot[ii].casadi_sparsity_in = &impl_jacFun_xdot_chain_nm6_sparsity_in;
				impl_jac_xdot[ii].casadi_sparsity_out = &impl_jacFun_xdot_chain_nm6_sparsity_out;
				impl_jac_xdot[ii].casadi_n_in = &impl_jacFun_xdot_chain_nm6_n_in;
				impl_jac_xdot[ii].casadi_n_out = &impl_jacFun_xdot_chain_nm6_n_out;
				impl_jac_u[ii].casadi_fun = &impl_jacFun_u_chain_nm6;
				impl_jac_u[ii].casadi_work = &impl_jacFun_u_chain_nm6_work;
				impl_jac_u[ii].casadi_sparsity_in = &impl_jacFun_u_chain_nm6_sparsity_in;
				impl_jac_u[ii].casadi_sparsity_out = &impl_jacFun_u_chain_nm6_sparsity_out;
				impl_jac_u[ii].casadi_n_in = &impl_jacFun_u_chain_nm6_n_in;
				impl_jac_u[ii].casadi_n_out = &impl_jacFun_u_chain_nm6_n_out;
			}
			break;
		default:
			printf("Problem size not available\n");
			exit(1);
			break;
	}
	return;
}



static void select_ls_stage_cost_jac_casadi(int indx, int N, int num_free_masses, external_function_casadi *ls_cost_jac)
{
	switch (num_free_masses)
	{
		case 1:
			if (indx < N)
			{
				ls_cost_jac->casadi_fun = &ls_cost_nm2;
				ls_cost_jac->casadi_work = &ls_cost_nm2_work;
				ls_cost_jac->casadi_sparsity_in = &ls_cost_nm2_sparsity_in;
				ls_cost_jac->casadi_sparsity_out = &ls_cost_nm2_sparsity_out;
				ls_cost_jac->casadi_n_in = &ls_cost_nm2_n_in;
				ls_cost_jac->casadi_n_out = &ls_cost_nm2_n_out;
			}
			else
			{
				ls_cost_jac->casadi_fun = &ls_costN_nm2;
				ls_cost_jac->casadi_work = &ls_costN_nm2_work;
				ls_cost_jac->casadi_sparsity_in = &ls_costN_nm2_sparsity_in;
				ls_cost_jac->casadi_sparsity_out = &ls_costN_nm2_sparsity_out;
				ls_cost_jac->casadi_n_in = &ls_costN_nm2_n_in;
				ls_cost_jac->casadi_n_out = &ls_costN_nm2_n_out;
			}
			break;
		case 2:
			if (indx < N)
			{
				ls_cost_jac->casadi_fun = &ls_cost_nm3;
				ls_cost_jac->casadi_work = &ls_cost_nm3_work;
				ls_cost_jac->casadi_sparsity_in = &ls_cost_nm3_sparsity_in;
				ls_cost_jac->casadi_sparsity_out = &ls_cost_nm3_sparsity_out;
				ls_cost_jac->casadi_n_in = &ls_cost_nm3_n_in;
				ls_cost_jac->casadi_n_out = &ls_cost_nm3_n_out;
			}
			else
			{
				ls_cost_jac->casadi_fun = &ls_costN_nm3;
				ls_cost_jac->casadi_work = &ls_costN_nm3_work;
				ls_cost_jac->casadi_sparsity_in = &ls_costN_nm3_sparsity_in;
				ls_cost_jac->casadi_sparsity_out = &ls_costN_nm3_sparsity_out;
				ls_cost_jac->casadi_n_in = &ls_costN_nm3_n_in;
				ls_cost_jac->casadi_n_out = &ls_costN_nm3_n_out;
			}
			break;
		case 3:
			if (indx < N)
			{
				ls_cost_jac->casadi_fun = &ls_cost_nm4;
				ls_cost_jac->casadi_work = &ls_cost_nm4_work;
				ls_cost_jac->casadi_sparsity_in = &ls_cost_nm4_sparsity_in;
				ls_cost_jac->casadi_sparsity_out = &ls_cost_nm4_sparsity_out;
				ls_cost_jac->casadi_n_in = &ls_cost_nm4_n_in;
				ls_cost_jac->casadi_n_out = &ls_cost_nm4_n_out;
			}
			else
			{
				ls_cost_jac->casadi_fun = &ls_costN_nm4;
				ls_cost_jac->casadi_work = &ls_costN_nm4_work;
				ls_cost_jac->casadi_sparsity_in = &ls_costN_nm4_sparsity_in;
				ls_cost_jac->casadi_sparsity_out = &ls_costN_nm4_sparsity_out;
				ls_cost_jac->casadi_n_in = &ls_costN_nm4_n_in;
				ls_cost_jac->casadi_n_out = &ls_costN_nm4_n_out;
			}
			break;
		case 4:
			if (indx < N)
			{
				ls_cost_jac->casadi_fun = &ls_cost_nm5;
				ls_cost_jac->casadi_work = &ls_cost_nm5_work;
				ls_cost_jac->casadi_sparsity_in = &ls_cost_nm5_sparsity_in;
				ls_cost_jac->casadi_sparsity_out = &ls_cost_nm5_sparsity_out;
				ls_cost_jac->casadi_n_in = &ls_cost_nm5_n_in;
				ls_cost_jac->casadi_n_out = &ls_cost_nm5_n_out;
			}
			else
			{
				ls_cost_jac->casadi_fun = &ls_costN_nm5;
				ls_cost_jac->casadi_work = &ls_costN_nm5_work;
				ls_cost_jac->casadi_sparsity_in = &ls_costN_nm5_sparsity_in;
				ls_cost_jac->casadi_sparsity_out = &ls_costN_nm5_sparsity_out;
				ls_cost_jac->casadi_n_in = &ls_costN_nm5_n_in;
				ls_cost_jac->casadi_n_out = &ls_costN_nm5_n_out;
			}
			break;
		case 5:
			if (indx < N)
			{
				ls_cost_jac->casadi_fun = &ls_cost_nm6;
				ls_cost_jac->casadi_work = &ls_cost_nm6_work;
				ls_cost_jac->casadi_sparsity_in = &ls_cost_nm6_sparsity_in;
				ls_cost_jac->casadi_sparsity_out = &ls_cost_nm6_sparsity_out;
				ls_cost_jac->casadi_n_in = &ls_cost_nm6_n_in;
				ls_cost_jac->casadi_n_out = &ls_cost_nm6_n_out;
			}
			else
			{
				ls_cost_jac->casadi_fun = &ls_costN_nm6;
				ls_cost_jac->casadi_work = &ls_costN_nm6_work;
				ls_cost_jac->casadi_sparsity_in = &ls_costN_nm6_sparsity_in;
				ls_cost_jac->casadi_sparsity_out = &ls_costN_nm6_sparsity_out;
				ls_cost_jac->casadi_n_in = &ls_costN_nm6_n_in;
				ls_cost_jac->casadi_n_out = &ls_costN_nm6_n_out;
			}
			break;
		default:
			printf("Problem size not available\n");
			exit(1);
			break;
	}

	return;
}



void read_initial_state(const int nx, const int num_free_masses, double *x0)
{
	double *ptr;
    switch (num_free_masses)
    {
        case 1:
            ptr = x0_nm2;
            break;
        case 2:
            ptr = x0_nm3;
            break;
        case 3:
            ptr = x0_nm4;
            break;
        case 4:
            ptr = x0_nm5;
            break;
        case 5:
            ptr = x0_nm6;
            break;
        default:
            printf("\nwrong number of free masses\n");
			exit(1);
            break;
    }
    for (int i = 0; i < nx; i++)
		x0[i] = ptr[i];
}



void read_final_state(const int nx, const int num_free_masses, double *xN)
{
	double *ptr;
    switch (num_free_masses)
    {
        case 1:
            ptr = xN_nm2;
            break;
        case 2:
            ptr = xN_nm3;
            break;
        case 3:
            ptr = xN_nm4;
            break;
        case 4:
            ptr = xN_nm5;
            break;
        case 5:
            ptr = xN_nm6;
            break;
        default:
            printf("\nwrong number of free masses\n");
			exit(1);
            break;
    }
    for (int i = 0; i < nx; i++)
		xN[i] = ptr[i];
}



// hand-generated external function for externally provided hessian and gradient
void ext_cost_nm2(void *fun, double *in, double *out)
{

	int ii;

	int nu = 3;
	int nx = 6;

	int nv = nu+nx;

	// ref
	double ref[nu+nx];
	for (ii=0; ii<nx; ii++)
		ref[ii] = xN_nm2[ii];
	for (ii=0; ii<nu; ii++)
		ref[nx+ii] = 0.0;

	// Hessian
	double *hess = out+nv;
	for (ii=0; ii<nv*nv; ii++)
		hess[ii] = 0.0;
	for (ii=0; ii<nx; ii++)
		hess[ii*(nv+1)] = 1e-2;
	for (; ii<nx+nu; ii++)
		hess[ii*(nv+1)] = 1.0;

	// gradient
	double *xu= in;
	double *grad = out;
	for (ii=0; ii<nv; ii++)
		grad[ii] = 0.0;
	for (ii=0; ii<nv; ii++)
		grad[ii] = hess[ii*(nv+1)] * (xu[ii] - ref[ii]);

	return;

}

void ext_cost_nm3(void *fun, double *in, double *out)
{

	int ii;

	int nu = 3;
	int nx = 12;

	int nv = nu+nx;

	// ref
	double ref[nu+nx];
	for (ii=0; ii<nx; ii++)
		ref[ii] = xN_nm3[ii];
	for (ii=0; ii<nu; ii++)
		ref[nx+ii] = 0.0;

	// Hessian
	double *hess = out+nv;
	for (ii=0; ii<nv*nv; ii++)
		hess[ii] = 0.0;
	for (ii=0; ii<nx; ii++)
		hess[ii*(nv+1)] = 1e-2;
	for (; ii<nx+nu; ii++)
		hess[ii*(nv+1)] = 1.0;

	// gradient
	double *xu= in;
	double *grad = out;
	for (ii=0; ii<nv; ii++)
		grad[ii] = 0.0;
	for (ii=0; ii<nv; ii++)
		grad[ii] = hess[ii*(nv+1)] * (xu[ii] - ref[ii]);

	return;

}

void ext_cost_nm4(void *fun, double *in, double *out)
{

	int ii;

	int nu = 3;
	int nx = 18;

	int nv = nu+nx;

	// ref
	double ref[nu+nx];
	for (ii=0; ii<nx; ii++)
		ref[ii] = xN_nm4[ii];
	for (ii=0; ii<nu; ii++)
		ref[nx+ii] = 0.0;

	// Hessian
	double *hess = out+nv;
	for (ii=0; ii<nv*nv; ii++)
		hess[ii] = 0.0;
	for (ii=0; ii<nx; ii++)
		hess[ii*(nv+1)] = 1e-2;
	for (; ii<nx+nu; ii++)
		hess[ii*(nv+1)] = 1.0;

	// gradient
	double *xu= in;
	double *grad = out;
	for (ii=0; ii<nv; ii++)
		grad[ii] = 0.0;
	for (ii=0; ii<nv; ii++)
		grad[ii] = hess[ii*(nv+1)] * (xu[ii] - ref[ii]);

	return;

}

void ext_cost_nm5(void *fun, double *in, double *out)
{

	int ii;

	int nu = 3;
	int nx = 24;

	int nv = nu+nx;

	// ref
	double ref[nu+nx];
	for (ii=0; ii<nx; ii++)
		ref[ii] = xN_nm5[ii];
	for (ii=0; ii<nu; ii++)
		ref[nx+ii] = 0.0;

	// Hessian
	double *hess = out+nv;
	for (ii=0; ii<nv*nv; ii++)
		hess[ii] = 0.0;
	for (ii=0; ii<nx; ii++)
		hess[ii*(nv+1)] = 1e-2;
	for (; ii<nx+nu; ii++)
		hess[ii*(nv+1)] = 1.0;

	// gradient
	double *xu= in;
	double *grad = out;
	for (ii=0; ii<nv; ii++)
		grad[ii] = 0.0;
	for (ii=0; ii<nv; ii++)
		grad[ii] = hess[ii*(nv+1)] * (xu[ii] - ref[ii]);

	return;

}

void ext_cost_nm6(void *fun, double *in, double *out)
{

	int ii;

	int nu = 3;
	int nx = 30;

	int nv = nu+nx;

	// ref
	double ref[nu+nx];
	for (ii=0; ii<nx; ii++)
		ref[ii] = xN_nm6[ii];
	for (ii=0; ii<nu; ii++)
		ref[nx+ii] = 0.0;

	// Hessian
	double *hess = out+nv;
	for (ii=0; ii<nv*nv; ii++)
		hess[ii] = 0.0;
	for (ii=0; ii<nx; ii++)
		hess[ii*(nv+1)] = 1e-2;
	for (; ii<nx+nu; ii++)
		hess[ii*(nv+1)] = 1.0;

	// gradient
	double *xu= in;
	double *grad = out;
	for (ii=0; ii<nv; ii++)
		grad[ii] = 0.0;
	for (ii=0; ii<nv; ii++)
		grad[ii] = hess[ii*(nv+1)] * (xu[ii] - ref[ii]);

	return;

}



// hand-wirtten box constraints on states as nonlinear constraints
void nonlin_constr_nm2(void *evaluate, double *in, double *out)
{

	int ii;

	int nu = 3;
	int nx = 6;

	int nv = nu+nx;
	int nh = nx;

	// fun
	double *fun = out;
	for (ii=0; ii<nx; ii++)
		fun[ii] = in[ii]; // x

	// jacobian
	double *jac = out+nh;
	for (ii=0; ii<nv*nh; ii++)
		jac[ii] = 0.0;
	for (ii=0; ii<nh; ii++)
		jac[ii*(nh+1)] = 1.0;

	return;

}

void nonlin_constr_nm3(void *evaluate, double *in, double *out)
{

	int ii;

	int nu = 3;
	int nx = 12;

	int nv = nu+nx;
	int nh = nx;

	// fun
	double *fun = out;
	for (ii=0; ii<nx; ii++)
		fun[ii] = in[ii]; // x

	// jacobian
	double *jac = out+nh;
	for (ii=0; ii<nv*nh; ii++)
		jac[ii] = 0.0;
	for (ii=0; ii<nh; ii++)
		jac[ii*(nh+1)] = 1.0;

	return;

}

void nonlin_constr_nm4(void *evaluate, double *in, double *out)
{

	int ii;

	int nu = 3;
	int nx = 18;

	int nv = nu+nx;
	int nh = nx;

	// fun
	double *fun = out;
	for (ii=0; ii<nx; ii++)
		fun[ii] = in[ii]; // x

	// jacobian
	double *jac = out+nh;
	for (ii=0; ii<nv*nh; ii++)
		jac[ii] = 0.0;
	for (ii=0; ii<nh; ii++)
		jac[ii*(nh+1)] = 1.0;

	return;

}

void nonlin_constr_nm5(void *evaluate, double *in, double *out)
{

	int ii;

	int nu = 3;
	int nx = 24;

	int nv = nu+nx;
	int nh = nx;

	// fun
	double *fun = out;
	for (ii=0; ii<nx; ii++)
		fun[ii] = in[ii]; // x

	// jacobian
	double *jac = out+nh;
	for (ii=0; ii<nv*nh; ii++)
		jac[ii] = 0.0;
	for (ii=0; ii<nh; ii++)
		jac[ii*(nh+1)] = 1.0;

	return;

}

void nonlin_constr_nm6(void *evaluate, double *in, double *out)
{

	int ii;

	int nu = 3;
	int nx = 30;

	int nv = nu+nx;
	int nh = nx;

	// fun
	double *fun = out;
	for (ii=0; ii<nx; ii++)
		fun[ii] = in[ii]; // x

	// jacobian
	double *jac = out+nh;
	for (ii=0; ii<nv*nh; ii++)
		jac[ii] = 0.0;
	for (ii=0; ii<nh; ii++)
		jac[ii*(nh+1)] = 1.0;

	return;

}



/************************************************
* main
************************************************/

// TODO(dimitris): compile on windows

int main()
{
    // _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);

	const int NMF = 3;  // number of free masses

    int NX = 6 * NMF;
    int NU = 3;

    /************************************************
    * problem dimensions
    ************************************************/

    int nx[NN + 1] = {0};
    int nu[NN + 1] = {0};
    int nbx[NN + 1] = {0};
    int nbu[NN + 1] = {0};
    int nb[NN + 1] = {0};
    int ng[NN + 1] = {0};
    int nh[NN + 1] = {0};
    int nq[NN + 1] = {0};
    int ns[NN+1] = {0};
	int ny[NN+1] = {0};

    nx[0] = NX;
    nu[0] = NU;
    nbx[0] = nx[0];
    nbu[0] = nu[0];
    nb[0] = nbu[0]+nbx[0];
	ng[0] = 0;
	nh[0] = 0;
	ny[0] = nx[0]+nu[0];

    for (int i = 1; i < NN; i++)
    {
        nx[i] = NX;
        nu[i] = NU;
        nbx[i] = NMF;
        nbu[i] = NU;
		nb[i] = nbu[i]+nbx[i];
		ng[i] = 0;
		nh[i] = 0;
		ny[i] = nx[i]+nu[i];
    }

    nx[NN] = NX;
    nu[NN] = 0;
    nbx[NN] = NX;
    nbu[NN] = 0;
    nb[NN] = nbu[NN]+nbx[NN];
	ng[NN] = 0;
	nh[NN] = 0;
	ny[NN] = nx[NN]+nu[NN];

    /************************************************
    * problem data
    ************************************************/

    double *x_end = malloc(sizeof(double)*NX);
    double *u_end = malloc(sizeof(double)*NU);

	for (int i = 0; i < NX; i++) x_end[i] = 0;
	for (int i = 0; i < NU; i++) u_end[i] = 0;

    double wall_pos = -0.01;
    double UMAX = 10;

	double x_pos_inf = +1e4;
	double x_neg_inf = -1e4;

    double *xref = malloc(NX*sizeof(double));
    read_final_state(NX, NMF, xref);

    double uref[3] = {0.0, 0.0, 0.0};

    double *diag_cost_x = malloc(NX*sizeof(double));

    for (int i = 0; i < NX; i++)
        diag_cost_x[i] = 1e-2;

    double diag_cost_u[3] = {1.0, 1.0, 1.0};


	// idxb0
    int idxb0[nb[0]];
    for (int i = 0; i < nb[0]; i++) idxb0[i] = i;

	// idxb1
	int idxb1[nb[1]];
    for (int i = 0; i < NU; i++) idxb1[i] = i;

    for (int i = 0; i < NMF; i++) idxb1[NU+i] = NU + 6*i + 1;

	// idxbN
	int idxbN[nb[NN]];
    for (int i = 0; i < nb[NN]; i++)
        idxbN[i] = i;

	// lb0, ub0
    double lb0[NX+NU], ub0[NX+NU];
    for (int i = 0; i < NU; i++)
	{
        lb0[i] = -UMAX;
        ub0[i] = +UMAX;
    }
    read_initial_state(NX, NMF, lb0+NU);
    read_initial_state(NX, NMF, ub0+NU);

	// lb1, ub1
    double lb1[NMF+NU], ub1[NMF+NU];
    for (int j = 0; j < NU; j++)
	{
        lb1[j] = -UMAX;  // umin
        ub1[j] = +UMAX;  // umax
    }
    for (int j = 0; j < NMF; j++)
	{
        lb1[NU+j] = wall_pos;  // wall position
        ub1[NU+j] = x_pos_inf;
    }

	// lbN, ubN
    double lbN[NX], ubN[NX];
    for (int i = 0; i < NX; i++)
	{
        lbN[i] = x_neg_inf;
        ubN[i] = x_pos_inf;
    }

    /************************************************
    * plan + config
    ************************************************/

	ocp_nlp_solver_plan *plan = ocp_nlp_plan_create(NN);

	// TODO(dimitris): not necessarily GN, depends on cost module
	plan->nlp_solver = SQP_GN;

	// NOTE(dimitris): switching between different objectives on each stage to test everything
	for (int i = 0; i <= NN; i++)
	{
		if (i < 3)
			plan->nlp_cost[i] = EXTERNALLY_PROVIDED;  // also implements linear LS for this example
		else if (i%2 == 0)
			plan->nlp_cost[i] = LINEAR_LS;
		else if (i%2 == 1)
			plan->nlp_cost[i] = NONLINEAR_LS;  // also implements linear LS for this example
	}

	plan->ocp_qp_solver_plan.qp_solver = PARTIAL_CONDENSING_HPIPM;

	// NOTE(dimitris): switching between different integrators on each stage to test everything
	for (int i = 0; i < NN; i++)
	{
		plan->nlp_dynamics[i] = CONTINUOUS_MODEL;
		if (i < 3)
			plan->sim_solver_plan[i].sim_solver = LIFTED_IRK;
		else if (i%2 == 0)
			plan->sim_solver_plan[i].sim_solver = IRK;
		else if (i%2 == 1)
			plan->sim_solver_plan[i].sim_solver = ERK;
	}

	// TODO(dimitris): fix minor memory leak here
	ocp_nlp_solver_config *config = ocp_nlp_config_create(*plan, NN);

    /************************************************
    * ocp_nlp_dims
    ************************************************/

	ocp_nlp_dims *dims = ocp_nlp_dims_create(config);
	ocp_nlp_dims_initialize(config, nx, nu, ny, nbx, nbu, ng, nh, ns, nq, dims);

    /************************************************
    * dynamics
    ************************************************/

	#if 0
	// NOTE(dimitris): temp code to test casadi integrator
	int integrator_nx = 12;
	int integrator_nu = 3;

	double *integrator_in = malloc(sizeof(double)*(integrator_nx+integrator_nu));
	double *integrator_out = malloc(sizeof(double)*(integrator_nx + integrator_nx*(integrator_nx+integrator_nu)));

	integrator_in[0] = 0.1;

	external_function_casadi casadi_integrator;
	casadi_integrator.casadi_fun = &casadi_erk4_chain_nm3;
	casadi_integrator.casadi_work = &casadi_erk4_chain_nm3_work;
	casadi_integrator.casadi_sparsity_in = &casadi_erk4_chain_nm3_sparsity_in;
	casadi_integrator.casadi_sparsity_out = &casadi_erk4_chain_nm3_sparsity_out;
	casadi_integrator.casadi_n_in = &casadi_erk4_chain_nm3_n_in;
	casadi_integrator.casadi_n_out = &casadi_erk4_chain_nm3_n_out;
	external_function_casadi_create(&casadi_integrator);

	d_print_mat(1, integrator_nx+integrator_nu, integrator_in, 1);

	casadi_integrator.evaluate(&casadi_integrator.evaluate, integrator_in, integrator_out);

	d_print_mat(1, integrator_nx, integrator_out, 1);

	d_print_mat(integrator_nx, integrator_nx + integrator_nu, integrator_out+integrator_nx, integrator_nx);

	free(integrator_in);
	free(integrator_out);
	exit(1);
	#endif

	// explicit
	external_function_casadi *forw_vde_casadi = malloc(NN*sizeof(external_function_casadi));
	external_function_casadi *jac_ode_casadi = malloc(NN*sizeof(external_function_casadi));

	// implicit
	external_function_casadi *impl_ode_casadi = malloc(NN*sizeof(external_function_casadi));
	external_function_casadi *impl_jac_x_casadi = malloc(NN*sizeof(external_function_casadi));
	external_function_casadi *impl_jac_xdot_casadi = malloc(NN*sizeof(external_function_casadi));
	external_function_casadi *impl_jac_u_casadi = malloc(NN*sizeof(external_function_casadi));

	select_dynamics_casadi(NN, NMF, forw_vde_casadi, jac_ode_casadi, impl_ode_casadi, impl_jac_x_casadi, impl_jac_xdot_casadi, impl_jac_u_casadi);

	// forw_vde
	external_function_casadi_create_array(NN, forw_vde_casadi);
	// jac_ode
	external_function_casadi_create_array(NN, jac_ode_casadi);

	// impl_ode
	external_function_casadi_create_array(NN, impl_ode_casadi);
	// jac_x
	external_function_casadi_create_array(NN, impl_jac_x_casadi);
	// jac_xdot
	external_function_casadi_create_array(NN, impl_jac_xdot_casadi);
	// jac_u
	external_function_casadi_create_array(NN, impl_jac_u_casadi);

    /************************************************
    * nonlinear least squares
    ************************************************/

	external_function_casadi *ls_cost_jac_casadi = malloc((NN+1)*sizeof(external_function_casadi));
	external_function_generic *ext_cost_generic = malloc(NN*sizeof(external_function_casadi));

	for (int i = 0; i <= NN; i++)
	{
		switch (plan->nlp_cost[i])
		{
			case LINEAR_LS:
				// do nothing
				break;

			case NONLINEAR_LS:
				select_ls_stage_cost_jac_casadi(i, NN, NMF, &ls_cost_jac_casadi[i]);
				external_function_casadi_create(&ls_cost_jac_casadi[i]);
				break;

			case EXTERNALLY_PROVIDED:
				// TODO(dimitris): move inside select_ls_stage_cost_jac_casadi?
				switch(NMF)
				{
					case 1:
						ext_cost_generic[i].evaluate = &ext_cost_nm2;
						break;
					case 2:
						ext_cost_generic[i].evaluate = &ext_cost_nm3;
						break;
					case 3:
						ext_cost_generic[i].evaluate = &ext_cost_nm4;
						break;
					case 4:
						ext_cost_generic[i].evaluate = &ext_cost_nm5;
						break;
					case 5:
						ext_cost_generic[i].evaluate = &ext_cost_nm6;
						break;
					default:
						printf("\next cost not implemented for this numer of masses\n\n");
						exit(1);
				}
				break;
		}
	}

    /************************************************
    * nlp_in
    ************************************************/

	ocp_nlp_in *nlp_in = ocp_nlp_in_create(config, dims);

	// sampling times
	for (int ii=0; ii<NN; ii++)
		nlp_in->Ts[ii] = TF/NN;

	// output definition: y = [x; u]

	/* cost */
	ocp_nlp_cost_ls_model *stage_cost_ls;
	ocp_nlp_cost_nls_model *stage_cost_nls;
	ocp_nlp_cost_external_model *stage_cost_external;

	for (int i = 0; i <= NN; i++)
	{
		switch (plan->nlp_cost[i])
		{
			case LINEAR_LS:

				stage_cost_ls = (ocp_nlp_cost_ls_model *) nlp_in->cost[i];

				// Cyt
				blasfeo_dgese(nu[i]+nx[i], ny[i], 0.0, &stage_cost_ls->Cyt, 0, 0);
					for (int j = 0; j < nu[i]; j++)
				BLASFEO_DMATEL(&stage_cost_ls->Cyt, j, nx[i]+j) = 1.0;
					for (int j = 0; j < nx[i]; j++)
				BLASFEO_DMATEL(&stage_cost_ls->Cyt, nu[i]+j, j) = 1.0;

				// W
				blasfeo_dgese(ny[i], ny[i], 0.0, &stage_cost_ls->W, 0, 0);
					for (int j = 0; j < nx[i]; j++)
				BLASFEO_DMATEL(&stage_cost_ls->W, j, j) = diag_cost_x[j];
					for (int j = 0; j < nu[i]; j++)
				BLASFEO_DMATEL(&stage_cost_ls->W, nx[i]+j, nx[i]+j) = diag_cost_u[j];

				// y_ref
				blasfeo_pack_dvec(nx[i], xref, &stage_cost_ls->y_ref, 0);
				blasfeo_pack_dvec(nu[i], uref, &stage_cost_ls->y_ref, nx[i]);
				break;

			case NONLINEAR_LS:

				stage_cost_nls = (ocp_nlp_cost_nls_model *) nlp_in->cost[i];

				// nls_jac
				stage_cost_nls->nls_jac = (external_function_generic *) &ls_cost_jac_casadi[i];

				// W
				blasfeo_dgese(ny[i], ny[i], 0.0, &stage_cost_nls->W, 0, 0);
				for (int j = 0; j < nx[i]; j++)
					BLASFEO_DMATEL(&stage_cost_nls->W, j, j) = diag_cost_x[j];
				for (int j = 0; j < nu[i]; j++)
					BLASFEO_DMATEL(&stage_cost_nls->W, nx[i]+j, nx[i]+j) = diag_cost_u[j];

				// y_ref
				blasfeo_pack_dvec(nx[i], xref, &stage_cost_nls->y_ref, 0);
				blasfeo_pack_dvec(nu[i], uref, &stage_cost_nls->y_ref, nx[i]);
				break;

			case EXTERNALLY_PROVIDED:

				stage_cost_external = (ocp_nlp_cost_external_model *) nlp_in->cost[i];

				stage_cost_external->ext_cost = &ext_cost_generic[i];

				assert(i < NN && "externally provided cost not implemented for last stage!");

				break;
		}
	}

	/* dynamics */
	int set_fun_status;

	for (int i=0; i<NN; i++)
	{
		switch (plan->nlp_dynamics[i])
		{
			case CONTINUOUS_MODEL:

				if (plan->sim_solver_plan[i].sim_solver == ERK)
				{
					set_fun_status = nlp_set_model_in_stage(config, nlp_in, i, "forward_vde", &forw_vde_casadi[i]);
					if (set_fun_status != 0) exit(1);
					set_fun_status = nlp_set_model_in_stage(config, nlp_in, i, "explicit_jacobian", &jac_ode_casadi[i]);
					if (set_fun_status != 0) exit(1);
				}
				else if (plan->sim_solver_plan[i].sim_solver == LIFTED_IRK)
				{
					set_fun_status = nlp_set_model_in_stage(config, nlp_in, i, "forward_vde", &forw_vde_casadi[i]);
					if (set_fun_status != 0) exit(1);
					set_fun_status = nlp_set_model_in_stage(config, nlp_in, i, "explicit_jacobian", &jac_ode_casadi[i]);
					if (set_fun_status != 0) exit(1);
				}
				else if (plan->sim_solver_plan[i].sim_solver == IRK)
				{
					set_fun_status = nlp_set_model_in_stage(config, nlp_in, i, "implicit_ode", &impl_ode_casadi[i]);
					if (set_fun_status != 0) exit(1);
					set_fun_status = nlp_set_model_in_stage(config, nlp_in, i, "implicit_jacobian_x", &impl_jac_x_casadi[i]);
					if (set_fun_status != 0) exit(1);
					set_fun_status = nlp_set_model_in_stage(config, nlp_in, i, "implicit_jacobian_xdot", &impl_jac_xdot_casadi[i]);
					if (set_fun_status != 0) exit(1);
					set_fun_status = nlp_set_model_in_stage(config, nlp_in, i, "implicit_jacobian_u", &impl_jac_u_casadi[i]);
					if (set_fun_status != 0) exit(1);
				}
				break;

			case DISCRETE_MODEL:
				// TODO(dimitris): add some discrete model instances to test
				break;
		}
	}

    nlp_in->freezeSens = false;
	// if (scheme > 2)
    //     nlp_in->freezeSens = true;

    /* constraints */
	ocp_nlp_constraints_model **constraints = (ocp_nlp_constraints_model **) nlp_in->constraints;

	// fist stage
	blasfeo_pack_dvec(nb[0], lb0, &constraints[0]->d, 0);
	blasfeo_pack_dvec(nb[0], ub0, &constraints[0]->d, nb[0]+ng[0]);
    constraints[0]->idxb = idxb0;

	// other stages
    for (int i = 1; i < NN; i++)
	{
		blasfeo_pack_dvec(nb[i], lb1, &constraints[i]->d, 0);
		blasfeo_pack_dvec(nb[i], ub1, &constraints[i]->d, nb[i]+ng[i]);
        constraints[i]->idxb = idxb1;
    }
	blasfeo_pack_dvec(nb[NN], lbN, &constraints[NN]->d, 0);
	blasfeo_pack_dvec(nb[NN], ubN, &constraints[NN]->d, nb[NN]+ng[NN]);
    constraints[NN]->idxb = idxbN;

#if 0
	for (int ii=0; ii<=NN; ii++)
	{
		blasfeo_print_dmat(nu[ii]+nx[ii], ng[ii], &constraints[ii]->DCt, 0, 0);
		blasfeo_print_tran_dvec(2*nb[ii]+2*ng[ii]+2*nh[ii], &constraints[ii]->d, 0);
	}
	exit(1);
#endif

    /************************************************
    * sqp opts
    ************************************************/

	void *nlp_opts = ocp_nlp_opts_create(config, dims);
	ocp_nlp_sqp_opts *sqp_opts = (ocp_nlp_sqp_opts *) nlp_opts;

    for (int i = 0; i < NN; ++i)
	{
		ocp_nlp_dynamics_cont_opts *dynamics_stage_opts = sqp_opts->dynamics[i];
        sim_rk_opts *sim_opts = dynamics_stage_opts->sim_solver;

		if (plan->sim_solver_plan[i].sim_solver == ERK)
		{
			sim_opts->ns = 4;
		}
		else if (plan->sim_solver_plan[i].sim_solver == LIFTED_IRK)
		{
			sim_opts->ns = 2;
		}
		else if (plan->sim_solver_plan[i].sim_solver == IRK)
		{
			sim_opts->ns = 2;
			sim_opts->jac_reuse = true;
		}
    }

    sqp_opts->maxIter = MAX_SQP_ITERS;
    sqp_opts->min_res_g = 1e-9;
    sqp_opts->min_res_b = 1e-9;
    sqp_opts->min_res_d = 1e-9;
    sqp_opts->min_res_m = 1e-9;

	// update after user-defined opts
	config->opts_update(config, dims, nlp_opts); // TODO ocp_nlp_opts_update

    /************************************************
    * ocp_nlp out
    ************************************************/

	ocp_nlp_out *nlp_out = ocp_nlp_out_create(config, dims);

	ocp_nlp_solver *solver = ocp_nlp_create(config, dims, nlp_opts);

    /************************************************
    * sqp solve
    ************************************************/

	int nmpc_problems = 100;

    int status;

    acados_timer timer;
    acados_tic(&timer);

    for (int rep = 0; rep < NREP; rep++)
    {
		// warm start output initial guess of solution
		for (int i=0; i<=NN; i++)
		{
			blasfeo_pack_dvec(nu[i], uref, nlp_out->ux+i, 0);
			blasfeo_pack_dvec(nx[i], xref, nlp_out->ux+i, nu[i]);
		}

   	 	for (int idx = 0; idx < nmpc_problems; idx++)
		{
        	status = ocp_nlp_solve(solver, nlp_in, nlp_out);

			// TODO(dimitris): why nbu[0] is 0??
			// printf("NU = %d %d\n", NU, dims->qp_solver->nbu[0]);

			// TODO(dimitris): simulate system instead of passing x[1] as next state
			blasfeo_unpack_dvec(dims->nx[1], &nlp_out->ux[1], dims->nu[1], lb0+NU);
			blasfeo_unpack_dvec(dims->nx[1], &nlp_out->ux[1], dims->nu[1], ub0+NU);

			blasfeo_pack_dvec(nb[0], lb0, &constraints[0]->d, 0);
			blasfeo_pack_dvec(nb[0], ub0, &constraints[0]->d, nb[0]+ng[0]);

			if (true)
			{
				shift_states(dims, nlp_out, x_end);
				shift_controls(dims, nlp_out, u_end);
			}
			// print info
			if (true)
			{
				printf("\nproblem #%d, status %d, iters %d\n", idx, status, ((ocp_nlp_sqp_memory *)solver->mem)->sqp_iter);
				printf("xsim = \n");
				blasfeo_print_tran_dvec(dims->nx[0], &nlp_out->ux[0], dims->nu[0]);
				printf("xref = \n");
				blasfeo_print_tran_dvec(dims->nx[0], &stage_cost_nls->y_ref, 0);
			}
		}
    }

    double time = acados_toc(&timer)/NREP;

    printf("\n\ntotal time = %f ms\n\n", time*1e3);

    /************************************************
    * free memory
    ************************************************/

	// TODO(dimitris): VALGRIND!
 	external_function_casadi_free(forw_vde_casadi);
	external_function_casadi_free(jac_ode_casadi);
	free(forw_vde_casadi);
	free(jac_ode_casadi);

	external_function_casadi_free(impl_ode_casadi);
	external_function_casadi_free(impl_jac_x_casadi);
	external_function_casadi_free(impl_jac_xdot_casadi);
	external_function_casadi_free(impl_jac_u_casadi);
	free(impl_ode_casadi);
	free(impl_jac_x_casadi);
	free(impl_jac_xdot_casadi);
	free(impl_jac_u_casadi);

	free(ls_cost_jac_casadi);
	free(ext_cost_generic);

	free(nlp_opts);
	free(nlp_in);
	free(nlp_out);
	free(solver);
	free(dims);
	free(config);
	free(plan);

	free(xref);
	free(diag_cost_x);

	free(x_end);
	free(u_end);

	for (int i = 0; i <= NN; i++)
	{
		switch (plan->nlp_cost[i])
		{
			case NONLINEAR_LS:
				external_function_casadi_free(&ls_cost_jac_casadi[i]);
				break;
			default:
				break;
		}
	}

	/************************************************
	* return
	************************************************/

	int check_sqp_iter = ((ocp_nlp_sqp_memory *)solver->mem)->sqp_iter;

	if (status == 0)
		printf("\nsuccess! (%d iter) \n\n", check_sqp_iter);
	else
		printf("\nfailure!\n\n");

	return 0;
}
