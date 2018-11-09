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
#include "acados/ocp_nlp/ocp_nlp_dynamics_disc.h"
#include "acados/ocp_nlp/ocp_nlp_constraints_bgh.h"

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
#define NREP 10

// constraints (at stage 0): 0 box, 1 general, 2 general+nonlinear
#define CONSTRAINTS 2


// TODO(dimitris): DOES THIS EVEN WORK ATM?
enum sensitivities_scheme {
    EXACT_NEWTON,
    INEXACT_NEWTON,
    INIS,
    FROZEN_INEXACT_NEWTON,
    FROZEN_INIS
};



static void print_problem_info(//enum sensitivities_scheme sensitivities_type,
                               const int num_free_masses) //, const int num_stages)
{
    // char scheme_name[MAX_STR_LEN];
    // switch (sensitivities_type) {
    //     case EXACT_NEWTON:
    //         snprintf(scheme_name, sizeof(scheme_name), "EXACT_NEWTON");
    //         break;
    //     case INEXACT_NEWTON:
    //         snprintf(scheme_name, sizeof(scheme_name), "INEXACT_NEWTON");
    //         break;
    //     case INIS:
    //         snprintf(scheme_name, sizeof(scheme_name), "INIS");
    //         break;
    //     case FROZEN_INEXACT_NEWTON:
    //         snprintf(scheme_name, sizeof(scheme_name), "FROZEN_INEXACT_NEWTON");
    //         break;
    //     case FROZEN_INIS:
    //         snprintf(scheme_name, sizeof(scheme_name), "FROZEN_INIS");
    //         break;
    //     default:
    //         printf("Chose sensitivities type not available");
    //         exit(1);
    // }
    // printf("\n----- NUMBER OF FREE MASSES = %d, stages = %d (%s) -----\n",
    //        num_free_masses, num_stages, scheme_name);
	printf("\n----- NUMBER OF FREE MASSES = %d -----\n", num_free_masses);
}



#if 0
// example of hand-generated external function
void ls_cost_jac_nm4(external_function_generic *fun, double *in, double *out)
{

	int ii;

	int nv = 21;

	double *d_ptr = out;

	for (ii=0; ii<nv; ii++)
		d_ptr[ii] = in[ii];
	d_ptr += nv;

	for (ii=0; ii<nv*nv; ii++)
		d_ptr[ii] = 0.0;
	for (ii=0; ii<nv; ii++)
		d_ptr[ii*(nv+1)] = 1.0;
	d_ptr += nv;

	return;

}
#endif



static void select_dynamics_casadi(int N, int num_free_masses,
	external_function_casadi *forw_vde,
	external_function_casadi *impl_ode_fun,
	external_function_casadi *impl_ode_fun_jac_x_xdot,
	external_function_casadi *impl_ode_fun_jac_x_xdot_u,
	external_function_casadi *impl_ode_jac_x_xdot_u,
	external_function_casadi *erk4_casadi)
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

				impl_ode_fun[ii].casadi_fun = &casadi_impl_ode_fun_chain_nm2;
				impl_ode_fun[ii].casadi_work = &casadi_impl_ode_fun_chain_nm2_work;
				impl_ode_fun[ii].casadi_sparsity_in = &casadi_impl_ode_fun_chain_nm2_sparsity_in;
				impl_ode_fun[ii].casadi_sparsity_out = &casadi_impl_ode_fun_chain_nm2_sparsity_out;
				impl_ode_fun[ii].casadi_n_in = &casadi_impl_ode_fun_chain_nm2_n_in;
				impl_ode_fun[ii].casadi_n_out = &casadi_impl_ode_fun_chain_nm2_n_out;

				impl_ode_fun_jac_x_xdot[ii].casadi_fun = &casadi_impl_ode_fun_jac_x_xdot_chain_nm2;
				impl_ode_fun_jac_x_xdot[ii].casadi_work = &casadi_impl_ode_fun_jac_x_xdot_chain_nm2_work;
				impl_ode_fun_jac_x_xdot[ii].casadi_sparsity_in = &casadi_impl_ode_fun_jac_x_xdot_chain_nm2_sparsity_in;
				impl_ode_fun_jac_x_xdot[ii].casadi_sparsity_out = &casadi_impl_ode_fun_jac_x_xdot_chain_nm2_sparsity_out;
				impl_ode_fun_jac_x_xdot[ii].casadi_n_in = &casadi_impl_ode_fun_jac_x_xdot_chain_nm2_n_in;
				impl_ode_fun_jac_x_xdot[ii].casadi_n_out = &casadi_impl_ode_fun_jac_x_xdot_chain_nm2_n_out;

				impl_ode_fun_jac_x_xdot_u[ii].casadi_fun = &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm2;
				impl_ode_fun_jac_x_xdot_u[ii].casadi_work = &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm2_work;
				impl_ode_fun_jac_x_xdot_u[ii].casadi_sparsity_in = &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm2_sparsity_in;
				impl_ode_fun_jac_x_xdot_u[ii].casadi_sparsity_out = &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm2_sparsity_out;
				impl_ode_fun_jac_x_xdot_u[ii].casadi_n_in = &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm2_n_in;
				impl_ode_fun_jac_x_xdot_u[ii].casadi_n_out = &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm2_n_out;

				impl_ode_jac_x_xdot_u[ii].casadi_fun = &casadi_impl_ode_jac_x_xdot_u_chain_nm2;
				impl_ode_jac_x_xdot_u[ii].casadi_work = &casadi_impl_ode_jac_x_xdot_u_chain_nm2_work;
				impl_ode_jac_x_xdot_u[ii].casadi_sparsity_in = &casadi_impl_ode_jac_x_xdot_u_chain_nm2_sparsity_in;
				impl_ode_jac_x_xdot_u[ii].casadi_sparsity_out = &casadi_impl_ode_jac_x_xdot_u_chain_nm2_sparsity_out;
				impl_ode_jac_x_xdot_u[ii].casadi_n_in = &casadi_impl_ode_jac_x_xdot_u_chain_nm2_n_in;
				impl_ode_jac_x_xdot_u[ii].casadi_n_out = &casadi_impl_ode_jac_x_xdot_u_chain_nm2_n_out;

				erk4_casadi[ii].casadi_fun = &casadi_erk4_chain_nm2;
				erk4_casadi[ii].casadi_work = &casadi_erk4_chain_nm2_work;
				erk4_casadi[ii].casadi_sparsity_in = &casadi_erk4_chain_nm2_sparsity_in;
				erk4_casadi[ii].casadi_sparsity_out = &casadi_erk4_chain_nm2_sparsity_out;
				erk4_casadi[ii].casadi_n_in = &casadi_erk4_chain_nm2_n_in;
				erk4_casadi[ii].casadi_n_out = &casadi_erk4_chain_nm2_n_out;
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

				impl_ode_fun[ii].casadi_fun = &casadi_impl_ode_fun_chain_nm3;
				impl_ode_fun[ii].casadi_work = &casadi_impl_ode_fun_chain_nm3_work;
				impl_ode_fun[ii].casadi_sparsity_in = &casadi_impl_ode_fun_chain_nm3_sparsity_in;
				impl_ode_fun[ii].casadi_sparsity_out = &casadi_impl_ode_fun_chain_nm3_sparsity_out;
				impl_ode_fun[ii].casadi_n_in = &casadi_impl_ode_fun_chain_nm3_n_in;
				impl_ode_fun[ii].casadi_n_out = &casadi_impl_ode_fun_chain_nm3_n_out;

				impl_ode_fun_jac_x_xdot[ii].casadi_fun = &casadi_impl_ode_fun_jac_x_xdot_chain_nm3;
				impl_ode_fun_jac_x_xdot[ii].casadi_work = &casadi_impl_ode_fun_jac_x_xdot_chain_nm3_work;
				impl_ode_fun_jac_x_xdot[ii].casadi_sparsity_in = &casadi_impl_ode_fun_jac_x_xdot_chain_nm3_sparsity_in;
				impl_ode_fun_jac_x_xdot[ii].casadi_sparsity_out = &casadi_impl_ode_fun_jac_x_xdot_chain_nm3_sparsity_out;
				impl_ode_fun_jac_x_xdot[ii].casadi_n_in = &casadi_impl_ode_fun_jac_x_xdot_chain_nm3_n_in;
				impl_ode_fun_jac_x_xdot[ii].casadi_n_out = &casadi_impl_ode_fun_jac_x_xdot_chain_nm3_n_out;

				impl_ode_fun_jac_x_xdot_u[ii].casadi_fun = &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm3;
				impl_ode_fun_jac_x_xdot_u[ii].casadi_work = &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm3_work;
				impl_ode_fun_jac_x_xdot_u[ii].casadi_sparsity_in = &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm3_sparsity_in;
				impl_ode_fun_jac_x_xdot_u[ii].casadi_sparsity_out = &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm3_sparsity_out;
				impl_ode_fun_jac_x_xdot_u[ii].casadi_n_in = &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm3_n_in;
				impl_ode_fun_jac_x_xdot_u[ii].casadi_n_out = &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm3_n_out;

				impl_ode_jac_x_xdot_u[ii].casadi_fun = &casadi_impl_ode_jac_x_xdot_u_chain_nm3;
				impl_ode_jac_x_xdot_u[ii].casadi_work = &casadi_impl_ode_jac_x_xdot_u_chain_nm3_work;
				impl_ode_jac_x_xdot_u[ii].casadi_sparsity_in = &casadi_impl_ode_jac_x_xdot_u_chain_nm3_sparsity_in;
				impl_ode_jac_x_xdot_u[ii].casadi_sparsity_out = &casadi_impl_ode_jac_x_xdot_u_chain_nm3_sparsity_out;
				impl_ode_jac_x_xdot_u[ii].casadi_n_in = &casadi_impl_ode_jac_x_xdot_u_chain_nm3_n_in;
				impl_ode_jac_x_xdot_u[ii].casadi_n_out = &casadi_impl_ode_jac_x_xdot_u_chain_nm3_n_out;

				erk4_casadi[ii].casadi_fun = &casadi_erk4_chain_nm3;
				erk4_casadi[ii].casadi_work = &casadi_erk4_chain_nm3_work;
				erk4_casadi[ii].casadi_sparsity_in = &casadi_erk4_chain_nm3_sparsity_in;
				erk4_casadi[ii].casadi_sparsity_out = &casadi_erk4_chain_nm3_sparsity_out;
				erk4_casadi[ii].casadi_n_in = &casadi_erk4_chain_nm3_n_in;
				erk4_casadi[ii].casadi_n_out = &casadi_erk4_chain_nm3_n_out;
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

				impl_ode_fun[ii].casadi_fun = &casadi_impl_ode_fun_chain_nm4;
				impl_ode_fun[ii].casadi_work = &casadi_impl_ode_fun_chain_nm4_work;
				impl_ode_fun[ii].casadi_sparsity_in = &casadi_impl_ode_fun_chain_nm4_sparsity_in;
				impl_ode_fun[ii].casadi_sparsity_out = &casadi_impl_ode_fun_chain_nm4_sparsity_out;
				impl_ode_fun[ii].casadi_n_in = &casadi_impl_ode_fun_chain_nm4_n_in;
				impl_ode_fun[ii].casadi_n_out = &casadi_impl_ode_fun_chain_nm4_n_out;

				impl_ode_fun_jac_x_xdot[ii].casadi_fun = &casadi_impl_ode_fun_jac_x_xdot_chain_nm4;
				impl_ode_fun_jac_x_xdot[ii].casadi_work = &casadi_impl_ode_fun_jac_x_xdot_chain_nm4_work;
				impl_ode_fun_jac_x_xdot[ii].casadi_sparsity_in = &casadi_impl_ode_fun_jac_x_xdot_chain_nm4_sparsity_in;
				impl_ode_fun_jac_x_xdot[ii].casadi_sparsity_out = &casadi_impl_ode_fun_jac_x_xdot_chain_nm4_sparsity_out;
				impl_ode_fun_jac_x_xdot[ii].casadi_n_in = &casadi_impl_ode_fun_jac_x_xdot_chain_nm4_n_in;
				impl_ode_fun_jac_x_xdot[ii].casadi_n_out = &casadi_impl_ode_fun_jac_x_xdot_chain_nm4_n_out;

				impl_ode_fun_jac_x_xdot_u[ii].casadi_fun = &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm4;
				impl_ode_fun_jac_x_xdot_u[ii].casadi_work = &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm4_work;
				impl_ode_fun_jac_x_xdot_u[ii].casadi_sparsity_in = &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm4_sparsity_in;
				impl_ode_fun_jac_x_xdot_u[ii].casadi_sparsity_out = &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm4_sparsity_out;
				impl_ode_fun_jac_x_xdot_u[ii].casadi_n_in = &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm4_n_in;
				impl_ode_fun_jac_x_xdot_u[ii].casadi_n_out = &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm4_n_out;

				impl_ode_jac_x_xdot_u[ii].casadi_fun = &casadi_impl_ode_jac_x_xdot_u_chain_nm4;
				impl_ode_jac_x_xdot_u[ii].casadi_work = &casadi_impl_ode_jac_x_xdot_u_chain_nm4_work;
				impl_ode_jac_x_xdot_u[ii].casadi_sparsity_in = &casadi_impl_ode_jac_x_xdot_u_chain_nm4_sparsity_in;
				impl_ode_jac_x_xdot_u[ii].casadi_sparsity_out = &casadi_impl_ode_jac_x_xdot_u_chain_nm4_sparsity_out;
				impl_ode_jac_x_xdot_u[ii].casadi_n_in = &casadi_impl_ode_jac_x_xdot_u_chain_nm4_n_in;
				impl_ode_jac_x_xdot_u[ii].casadi_n_out = &casadi_impl_ode_jac_x_xdot_u_chain_nm4_n_out;

				erk4_casadi[ii].casadi_fun = &casadi_erk4_chain_nm4;
				erk4_casadi[ii].casadi_work = &casadi_erk4_chain_nm4_work;
				erk4_casadi[ii].casadi_sparsity_in = &casadi_erk4_chain_nm4_sparsity_in;
				erk4_casadi[ii].casadi_sparsity_out = &casadi_erk4_chain_nm4_sparsity_out;
				erk4_casadi[ii].casadi_n_in = &casadi_erk4_chain_nm4_n_in;
				erk4_casadi[ii].casadi_n_out = &casadi_erk4_chain_nm4_n_out;
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

				impl_ode_fun[ii].casadi_fun = &casadi_impl_ode_fun_chain_nm5;
				impl_ode_fun[ii].casadi_work = &casadi_impl_ode_fun_chain_nm5_work;
				impl_ode_fun[ii].casadi_sparsity_in = &casadi_impl_ode_fun_chain_nm5_sparsity_in;
				impl_ode_fun[ii].casadi_sparsity_out = &casadi_impl_ode_fun_chain_nm5_sparsity_out;
				impl_ode_fun[ii].casadi_n_in = &casadi_impl_ode_fun_chain_nm5_n_in;
				impl_ode_fun[ii].casadi_n_out = &casadi_impl_ode_fun_chain_nm5_n_out;

				impl_ode_fun_jac_x_xdot[ii].casadi_fun = &casadi_impl_ode_fun_jac_x_xdot_chain_nm5;
				impl_ode_fun_jac_x_xdot[ii].casadi_work = &casadi_impl_ode_fun_jac_x_xdot_chain_nm5_work;
				impl_ode_fun_jac_x_xdot[ii].casadi_sparsity_in = &casadi_impl_ode_fun_jac_x_xdot_chain_nm5_sparsity_in;
				impl_ode_fun_jac_x_xdot[ii].casadi_sparsity_out = &casadi_impl_ode_fun_jac_x_xdot_chain_nm5_sparsity_out;
				impl_ode_fun_jac_x_xdot[ii].casadi_n_in = &casadi_impl_ode_fun_jac_x_xdot_chain_nm5_n_in;
				impl_ode_fun_jac_x_xdot[ii].casadi_n_out = &casadi_impl_ode_fun_jac_x_xdot_chain_nm5_n_out;

				impl_ode_fun_jac_x_xdot_u[ii].casadi_fun = &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm5;
				impl_ode_fun_jac_x_xdot_u[ii].casadi_work = &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm5_work;
				impl_ode_fun_jac_x_xdot_u[ii].casadi_sparsity_in = &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm5_sparsity_in;
				impl_ode_fun_jac_x_xdot_u[ii].casadi_sparsity_out = &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm5_sparsity_out;
				impl_ode_fun_jac_x_xdot_u[ii].casadi_n_in = &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm5_n_in;
				impl_ode_fun_jac_x_xdot_u[ii].casadi_n_out = &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm5_n_out;

				impl_ode_jac_x_xdot_u[ii].casadi_fun = &casadi_impl_ode_jac_x_xdot_u_chain_nm5;
				impl_ode_jac_x_xdot_u[ii].casadi_work = &casadi_impl_ode_jac_x_xdot_u_chain_nm5_work;
				impl_ode_jac_x_xdot_u[ii].casadi_sparsity_in = &casadi_impl_ode_jac_x_xdot_u_chain_nm5_sparsity_in;
				impl_ode_jac_x_xdot_u[ii].casadi_sparsity_out = &casadi_impl_ode_jac_x_xdot_u_chain_nm5_sparsity_out;
				impl_ode_jac_x_xdot_u[ii].casadi_n_in = &casadi_impl_ode_jac_x_xdot_u_chain_nm5_n_in;
				impl_ode_jac_x_xdot_u[ii].casadi_n_out = &casadi_impl_ode_jac_x_xdot_u_chain_nm5_n_out;

				// erk4_casadi[ii].casadi_fun = &casadi_erk4_chain_nm5;
				// erk4_casadi[ii].casadi_work = &casadi_erk4_chain_nm5_work;
				// erk4_casadi[ii].casadi_sparsity_in = &casadi_erk4_chain_nm5_sparsity_in;
				// erk4_casadi[ii].casadi_sparsity_out = &casadi_erk4_chain_nm5_sparsity_out;
				// erk4_casadi[ii].casadi_n_in = &casadi_erk4_chain_nm5_n_in;
				// erk4_casadi[ii].casadi_n_out = &casadi_erk4_chain_nm5_n_out;
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

				impl_ode_fun[ii].casadi_fun = &casadi_impl_ode_fun_chain_nm6;
				impl_ode_fun[ii].casadi_work = &casadi_impl_ode_fun_chain_nm6_work;
				impl_ode_fun[ii].casadi_sparsity_in = &casadi_impl_ode_fun_chain_nm6_sparsity_in;
				impl_ode_fun[ii].casadi_sparsity_out = &casadi_impl_ode_fun_chain_nm6_sparsity_out;
				impl_ode_fun[ii].casadi_n_in = &casadi_impl_ode_fun_chain_nm6_n_in;
				impl_ode_fun[ii].casadi_n_out = &casadi_impl_ode_fun_chain_nm6_n_out;

				impl_ode_fun_jac_x_xdot[ii].casadi_fun = &casadi_impl_ode_fun_jac_x_xdot_chain_nm6;
				impl_ode_fun_jac_x_xdot[ii].casadi_work = &casadi_impl_ode_fun_jac_x_xdot_chain_nm6_work;
				impl_ode_fun_jac_x_xdot[ii].casadi_sparsity_in = &casadi_impl_ode_fun_jac_x_xdot_chain_nm6_sparsity_in;
				impl_ode_fun_jac_x_xdot[ii].casadi_sparsity_out = &casadi_impl_ode_fun_jac_x_xdot_chain_nm6_sparsity_out;
				impl_ode_fun_jac_x_xdot[ii].casadi_n_in = &casadi_impl_ode_fun_jac_x_xdot_chain_nm6_n_in;
				impl_ode_fun_jac_x_xdot[ii].casadi_n_out = &casadi_impl_ode_fun_jac_x_xdot_chain_nm6_n_out;

				impl_ode_fun_jac_x_xdot_u[ii].casadi_fun = &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm6;
				impl_ode_fun_jac_x_xdot_u[ii].casadi_work = &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm6_work;
				impl_ode_fun_jac_x_xdot_u[ii].casadi_sparsity_in = &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm6_sparsity_in;
				impl_ode_fun_jac_x_xdot_u[ii].casadi_sparsity_out = &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm6_sparsity_out;
				impl_ode_fun_jac_x_xdot_u[ii].casadi_n_in = &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm6_n_in;
				impl_ode_fun_jac_x_xdot_u[ii].casadi_n_out = &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm6_n_out;

				impl_ode_jac_x_xdot_u[ii].casadi_fun = &casadi_impl_ode_jac_x_xdot_u_chain_nm6;
				impl_ode_jac_x_xdot_u[ii].casadi_work = &casadi_impl_ode_jac_x_xdot_u_chain_nm6_work;
				impl_ode_jac_x_xdot_u[ii].casadi_sparsity_in = &casadi_impl_ode_jac_x_xdot_u_chain_nm6_sparsity_in;
				impl_ode_jac_x_xdot_u[ii].casadi_sparsity_out = &casadi_impl_ode_jac_x_xdot_u_chain_nm6_sparsity_out;
				impl_ode_jac_x_xdot_u[ii].casadi_n_in = &casadi_impl_ode_jac_x_xdot_u_chain_nm6_n_in;
				impl_ode_jac_x_xdot_u[ii].casadi_n_out = &casadi_impl_ode_jac_x_xdot_u_chain_nm6_n_out;

				// erk4_casadi[ii].casadi_fun = &casadi_erk4_chain_nm6;
				// erk4_casadi[ii].casadi_work = &casadi_erk4_chain_nm6_work;
				// erk4_casadi[ii].casadi_sparsity_in = &casadi_erk4_chain_nm6_sparsity_in;
				// erk4_casadi[ii].casadi_sparsity_out = &casadi_erk4_chain_nm6_sparsity_out;
				// erk4_casadi[ii].casadi_n_in = &casadi_erk4_chain_nm6_n_in;
				// erk4_casadi[ii].casadi_n_out = &casadi_erk4_chain_nm6_n_out;
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


#if 0
static void select_ls_cost_jac_casadi(int N, int num_free_masses, external_function_casadi *ls_cost_jac)
{
	for (int ii = 0; ii <= N; ii++)
		select_ls_stage_cost_jac_casadi(ii, N, num_free_masses, &ls_cost_jac[ii]);
}
#endif



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
void ext_cost_nm2(void *fun, ext_fun_arg_t *type_in, void **in, ext_fun_arg_t *type_out, void **out)
{

	int ii;

	int nu = 3;
	int nx = 6;

	int nv = nu+nx;

	// ref
	double *ref = calloc(nx+nu, sizeof(double));
	for (ii=0; ii<nu; ii++)
		ref[ii] = 0.0;
	for (ii=0; ii<nx; ii++)
		ref[nu+ii] = xN_nm2[ii];

	// Hessian
	double *hess = out[1];
	for (ii=0; ii<nv*nv; ii++)
		hess[ii] = 0.0;
	for (ii=0; ii<nu; ii++)
		hess[ii*(nv+1)] = 1.0;
	for (; ii<nu+nx; ii++)
		hess[ii*(nv+1)] = 1e-2;

	// gradient
	double *ux = in[0];
	double *grad = out[0];
	for (ii=0; ii<nv; ii++)
		grad[ii] = 0.0;
	for (ii=0; ii<nv; ii++)
		grad[ii] = hess[ii*(nv+1)] * (ux[ii] - ref[ii]);

    free(ref);
	return;

}

void ext_cost_nm3(void *fun, ext_fun_arg_t *type_in, void **in, ext_fun_arg_t *type_out, void **out)
{

	int ii;

	int nu = 3;
	int nx = 12;

	int nv = nu+nx;

	// ref
    double *ref = calloc(nx+nu, sizeof(double));
	for (ii=0; ii<nu; ii++)
		ref[ii] = 0.0;
	for (ii=0; ii<nx; ii++)
		ref[nu+ii] = xN_nm3[ii];

	// Hessian
	double *hess = out[1];
	for (ii=0; ii<nv*nv; ii++)
		hess[ii] = 0.0;
	for (ii=0; ii<nu; ii++)
		hess[ii*(nv+1)] = 1.0;
	for (; ii<nu+nx; ii++)
		hess[ii*(nv+1)] = 1e-2;

	// gradient
	double *ux = in[0];
	double *grad = out[0];
	for (ii=0; ii<nv; ii++)
		grad[ii] = 0.0;
	for (ii=0; ii<nv; ii++)
		grad[ii] = hess[ii*(nv+1)] * (ux[ii] - ref[ii]);

    free(ref);
	return;

}

void ext_cost_nm4(void *fun, ext_fun_arg_t *type_in, void **in, ext_fun_arg_t *type_out, void **out)
{

	int ii;

	int nu = 3;
	int nx = 18;

	int nv = nu+nx;

	// ref
    double *ref = calloc(nx+nu, sizeof(double));
	for (ii=0; ii<nu; ii++)
		ref[ii] = 0.0;
	for (ii=0; ii<nx; ii++)
		ref[nu+ii] = xN_nm4[ii];

	// Hessian
	double *hess = out[1];
	for (ii=0; ii<nv*nv; ii++)
		hess[ii] = 0.0;
	for (ii=0; ii<nu; ii++)
		hess[ii*(nv+1)] = 1.0;
	for (; ii<nu+nx; ii++)
		hess[ii*(nv+1)] = 1e-2;

	// gradient
	double *ux = in[0];
	double *grad = out[0];
	for (ii=0; ii<nv; ii++)
		grad[ii] = 0.0;
	for (ii=0; ii<nv; ii++)
		grad[ii] = hess[ii*(nv+1)] * (ux[ii] - ref[ii]);

    free(ref);
	return;

}

void ext_cost_nm5(void *fun, ext_fun_arg_t *type_in, void **in, ext_fun_arg_t *type_out, void **out)
{

	int ii;

	int nu = 3;
	int nx = 24;

	int nv = nu+nx;

	// ref
    double *ref = calloc(nx+nu, sizeof(double));
	for (ii=0; ii<nu; ii++)
		ref[ii] = 0.0;
	for (ii=0; ii<nx; ii++)
		ref[nu+ii] = xN_nm5[ii];

	// Hessian
	double *hess = out[1];
	for (ii=0; ii<nv*nv; ii++)
		hess[ii] = 0.0;
	for (ii=0; ii<nu; ii++)
		hess[ii*(nv+1)] = 1.0;
	for (; ii<nu+nx; ii++)
		hess[ii*(nv+1)] = 1e-2;

	// gradient
	double *ux = in[0];
	double *grad = out[0];
	for (ii=0; ii<nv; ii++)
		grad[ii] = 0.0;
	for (ii=0; ii<nv; ii++)
		grad[ii] = hess[ii*(nv+1)] * (ux[ii] - ref[ii]);

    free(ref);
	return;

}

void ext_cost_nm6(void *fun, ext_fun_arg_t *type_in, void **in, ext_fun_arg_t *type_out, void **out)
{

	int ii;

	int nu = 3;
	int nx = 30;

	int nv = nu+nx;

	// ref
    double *ref = calloc(nx+nu, sizeof(double));
	for (ii=0; ii<nu; ii++)
		ref[ii] = 0.0;
	for (ii=0; ii<nx; ii++)
		ref[nu+ii] = xN_nm6[ii];

	// Hessian
	double *hess = out[1];
	for (ii=0; ii<nv*nv; ii++)
		hess[ii] = 0.0;
	for (ii=0; ii<nu; ii++)
		hess[ii*(nv+1)] = 1.0;
	for (; ii<nu+nx; ii++)
		hess[ii*(nv+1)] = 1e-2;

	// gradient
	double *ux = in[0];
	double *grad = out[0];
	for (ii=0; ii<nv; ii++)
		grad[ii] = 0.0;
	for (ii=0; ii<nv; ii++)
		grad[ii] = hess[ii*(nv+1)] * (ux[ii] - ref[ii]);

    free(ref);
	return;

}



// hand-wirtten box constraints on states as nonlinear constraints
void nonlin_constr_nm2(void *evaluate, ext_fun_arg_t *type_in, void **in, ext_fun_arg_t *type_out, void **out)
{

	int ii;

	int nu = 3;
	int nx = 6;

	int nh = nx;

	// fun
	struct blasfeo_dvec_args *fun_args = out[0];
	struct blasfeo_dvec *fun = fun_args->x;
	int xi = fun_args->xi;
	struct blasfeo_dvec *ux = in[0];
	blasfeo_dveccp(nx, ux, nu, fun, xi);

	// jacobian
	struct blasfeo_dmat_args *jac_args = out[1];
	struct blasfeo_dmat *jac = jac_args->A;
	int ai = jac_args->ai;
	int aj = jac_args->aj;
	blasfeo_dgese(nu+nx, nh, 0.0, jac, ai, aj);
	for (ii=0; ii<nh; ii++)
		BLASFEO_DMATEL(jac, ai+nu+ii, aj+ii) = 1.0;

	return;

}

void nonlin_constr_nm3(void *evaluate, ext_fun_arg_t *type_in, void **in, ext_fun_arg_t *type_out, void **out)
{

	int ii;

	int nu = 3;
	int nx = 12;

	int nh = nx;

	// fun
	struct blasfeo_dvec_args *fun_args = out[0];
	struct blasfeo_dvec *fun = fun_args->x;
	int xi = fun_args->xi;
	struct blasfeo_dvec *ux = in[0];
	blasfeo_dveccp(nx, ux, nu, fun, xi);

	// jacobian
	struct blasfeo_dmat_args *jac_args = out[1];
	struct blasfeo_dmat *jac = jac_args->A;
	int ai = jac_args->ai;
	int aj = jac_args->aj;
	blasfeo_dgese(nu+nx, nh, 0.0, jac, ai, aj);
	for (ii=0; ii<nh; ii++)
		BLASFEO_DMATEL(jac, ai+nu+ii, aj+ii) = 1.0;

	return;

}

void nonlin_constr_nm4(void *evaluate, ext_fun_arg_t *type_in, void **in, ext_fun_arg_t *type_out, void **out)
{

	int ii;

	int nu = 3;
	int nx = 18;

	int nh = nx;

	// fun
	struct blasfeo_dvec_args *fun_args = out[0];
	struct blasfeo_dvec *fun = fun_args->x;
	int xi = fun_args->xi;
	struct blasfeo_dvec *ux = in[0];
	blasfeo_dveccp(nx, ux, nu, fun, xi);

	// jacobian
	struct blasfeo_dmat_args *jac_args = out[1];
	struct blasfeo_dmat *jac = jac_args->A;
	int ai = jac_args->ai;
	int aj = jac_args->aj;
	blasfeo_dgese(nu+nx, nh, 0.0, jac, ai, aj);
	for (ii=0; ii<nh; ii++)
		BLASFEO_DMATEL(jac, ai+nu+ii, aj+ii) = 1.0;

	return;

}

void nonlin_constr_nm5(void *evaluate, ext_fun_arg_t *type_in, void **in, ext_fun_arg_t *type_out, void **out)
{

	int ii;

	int nu = 3;
	int nx = 24;

	int nh = nx;

	// fun
	struct blasfeo_dvec_args *fun_args = out[0];
	struct blasfeo_dvec *fun = fun_args->x;
	int xi = fun_args->xi;
	struct blasfeo_dvec *ux = in[0];
	blasfeo_dveccp(nx, ux, nu, fun, xi);

	// jacobian
	struct blasfeo_dmat_args *jac_args = out[1];
	struct blasfeo_dmat *jac = jac_args->A;
	int ai = jac_args->ai;
	int aj = jac_args->aj;
	blasfeo_dgese(nu+nx, nh, 0.0, jac, ai, aj);
	for (ii=0; ii<nh; ii++)
		BLASFEO_DMATEL(jac, ai+nu+ii, aj+ii) = 1.0;

	return;

}

void nonlin_constr_nm6(void *evaluate, ext_fun_arg_t *type_in, void **in, ext_fun_arg_t *type_out, void **out)
{

	int ii;

	int nu = 3;
	int nx = 30;

	int nh = nx;

	// fun
	struct blasfeo_dvec_args *fun_args = out[0];
	struct blasfeo_dvec *fun = fun_args->x;
	int xi = fun_args->xi;
	struct blasfeo_dvec *ux = in[0];
	blasfeo_dveccp(nx, ux, nu, fun, xi);

	// jacobian
	struct blasfeo_dmat_args *jac_args = out[1];
	struct blasfeo_dmat *jac = jac_args->A;
	int ai = jac_args->ai;
	int aj = jac_args->aj;
	blasfeo_dgese(nu+nx, nh, 0.0, jac, ai, aj);
	for (ii=0; ii<nh; ii++)
		BLASFEO_DMATEL(jac, ai+nu+ii, aj+ii) = 1.0;

	return;

}



/************************************************
* main
************************************************/

int main()
{
    // _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);

	const int NMF = 3;  // number of free masses: actually one more is used: possible values are 1,2,3,4,5

    print_problem_info(NMF);

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
	int nz[NN+1] = {0};

    nx[0] = NX;
    nu[0] = NU;
#if CONSTRAINTS==0 // box
    nbx[0] = nx[0];
    nbu[0] = nu[0];
    nb[0] = nbu[0]+nbx[0];
	ng[0] = 0;
	nh[0] = 0;
#elif CONSTRAINTS==1 // general
    nbx[0] = 0;
    nbu[0] = 0;
	nb[0] = 0;
    ng[0] = nu[0]+nx[0];
	nh[0] = 0;
#else // general+nonlinear constraints
    nbx[0] = 0;
    nbu[0] = 0;
	nb[0] = 0;
    ng[0] = nu[0];
	nh[0] = nx[0];
#endif
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
	int *idxb0 = malloc(nb[0]*sizeof(int));

    for (int i = 0; i < nb[0]; i++) idxb0[i] = i;

	// idxb1
	int *idxb1 = malloc(nb[1]*sizeof(int));
    for (int i = 0; i < NU; i++) idxb1[i] = i;

    for (int i = 0; i < NMF; i++) idxb1[NU+i] = NU + 6*i + 1;

	// idxbN
	int *idxbN = malloc(nb[NN]*sizeof(int));
    for (int i = 0; i < nb[NN]; i++)
		idxbN[i] = i;

	// lb0, ub0
	double *lb0 = malloc((NX+NU)*sizeof(double));
	double *ub0 = malloc((NX+NU)*sizeof(double));

    for (int i = 0; i < NU; i++)
	{
        lb0[i] = -UMAX;
        ub0[i] = +UMAX;
    }
    read_initial_state(NX, NMF, lb0+NU);
    read_initial_state(NX, NMF, ub0+NU);

	// lb1, ub1
	double *lb1 = malloc((NMF+NU)*sizeof(double));
	double *ub1 = malloc((NMF+NU)*sizeof(double));

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
	double *lbN = malloc(NX*sizeof(double));
	double *ubN = malloc(NX*sizeof(double));

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
	plan->nlp_solver = SQP;

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
	// plan->ocp_qp_solver_plan.qp_solver = FULL_CONDENSING_HPIPM;
	// plan->ocp_qp_solver_plan.qp_solver = FULL_CONDENSING_QPOASES;
	// plan->ocp_qp_solver_plan.qp_solver = FULL_CONDENSING_OOQP;
	// plan->ocp_qp_solver_plan.qp_solver = PARTIAL_CONDENSING_OOQP;

	// NOTE(dimitris): switching between different integrators on each stage to test everything
	for (int i = 0; i < NN; i++)
	{
		if (i < NN-4)
		{
			plan->nlp_dynamics[i] = CONTINUOUS_MODEL;
			if (i < 3)
				plan->sim_solver_plan[i].sim_solver = LIFTED_IRK;
			else if (i%3 == 0)
				plan->sim_solver_plan[i].sim_solver = IRK;
			else if (i%3 == 1)
				plan->sim_solver_plan[i].sim_solver = ERK;
			else if (i%3 == 2)
				plan->sim_solver_plan[i].sim_solver = LIFTED_IRK;
		}
		else
		{
			plan->nlp_dynamics[i] = DISCRETE_MODEL;
		}
	}

	for (int i = 0; i <= NN; i++)
		plan->nlp_constraints[i] = BGH;

	// TODO(dimitris): fix minor memory leak here
	ocp_nlp_solver_config *config = ocp_nlp_config_create(*plan, NN);

    /************************************************
    * ocp_nlp_dims
    ************************************************/

	ocp_nlp_dims *dims = ocp_nlp_dims_create(config);
	ocp_nlp_dims_initialize(config, nx, nu, ny, nbx, nbu, ng, nh, nq, ns, nz, dims);

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
	external_function_casadi *expl_vde_for = malloc(NN*sizeof(external_function_casadi));

	// implicit
	external_function_casadi *impl_ode_fun = malloc(NN*sizeof(external_function_casadi));
	external_function_casadi *impl_ode_fun_jac_x_xdot = malloc(NN*sizeof(external_function_casadi));
	external_function_casadi *impl_ode_fun_jac_x_xdot_u = malloc(NN*sizeof(external_function_casadi));
	external_function_casadi *impl_ode_jac_x_xdot_u = malloc(NN*sizeof(external_function_casadi));

	// discrete model
	external_function_casadi *erk4_casadi = malloc(NN*sizeof(external_function_casadi));

	select_dynamics_casadi(NN, NMF, expl_vde_for, impl_ode_fun, impl_ode_fun_jac_x_xdot, impl_ode_fun_jac_x_xdot_u, impl_ode_jac_x_xdot_u, erk4_casadi);

	// forw_vde
	external_function_casadi_create_array(NN, expl_vde_for);

	// impl_ode
	external_function_casadi_create_array(NN, impl_ode_fun);
	//
	external_function_casadi_create_array(NN, impl_ode_fun_jac_x_xdot);
	//
	external_function_casadi_create_array(NN, impl_ode_fun_jac_x_xdot_u);
	//
	external_function_casadi_create_array(NN, impl_ode_jac_x_xdot_u);

	if (NMF<4)
	{
		// discrete model supported
		external_function_casadi_create_array(NN, erk4_casadi);
	}
	else
	{
		printf("\nERROR: in this case (Number of free masses (NMF) > 3)discrete model is not supported, commented in cmake to speed up compilation\n");
		exit(1);
	}


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
    * nonlinear constraints
    ************************************************/

#if CONSTRAINTS==2
	external_function_generic nonlin_constr_generic;

	switch(NMF)
	{
		case 1:
			nonlin_constr_generic.evaluate = &nonlin_constr_nm2;
			break;
		case 2:
			nonlin_constr_generic.evaluate = &nonlin_constr_nm3;
			break;
		case 3:
			nonlin_constr_generic.evaluate = &nonlin_constr_nm4;
			break;
		case 4:
			nonlin_constr_generic.evaluate = &nonlin_constr_nm5;
			break;
		case 5:
			nonlin_constr_generic.evaluate = &nonlin_constr_nm6;
			break;
		default:
			printf("\nnonlin constr not implemented for this numer of masses\n\n");
			exit(1);
	}
#endif

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

	// TODO(dimitris): remove after setting function via nlp interface
	ocp_nlp_dynamics_disc_model *dynamics;

	for (int i=0; i<NN; i++)
	{
		switch (plan->nlp_dynamics[i])
		{
			case CONTINUOUS_MODEL:

				if (plan->sim_solver_plan[i].sim_solver == ERK)
				{
					set_fun_status = nlp_set_model_in_stage(config, nlp_in, i, "expl_vde_for", &expl_vde_for[i]);
					if (set_fun_status != 0) exit(1);
				}
				else if (plan->sim_solver_plan[i].sim_solver == IRK)
				{
					set_fun_status = nlp_set_model_in_stage(config, nlp_in, i, "impl_ode_fun", &impl_ode_fun[i]);
					if (set_fun_status != 0) exit(1);
					set_fun_status = nlp_set_model_in_stage(config, nlp_in, i, "impl_ode_fun_jac_x_xdot", &impl_ode_fun_jac_x_xdot[i]);
					if (set_fun_status != 0) exit(1);
					set_fun_status = nlp_set_model_in_stage(config, nlp_in, i, "impl_ode_jac_x_xdot_u", &impl_ode_jac_x_xdot_u[i]);
					if (set_fun_status != 0) exit(1);
				}
				else if (plan->sim_solver_plan[i].sim_solver == LIFTED_IRK)
				{
					set_fun_status = nlp_set_model_in_stage(config, nlp_in, i, "impl_ode_fun", &impl_ode_fun[i]);
					if (set_fun_status != 0) exit(1);
					set_fun_status = nlp_set_model_in_stage(config, nlp_in, i, "impl_ode_fun_jac_x_xdot_u", &impl_ode_fun_jac_x_xdot_u[i]);
					if (set_fun_status != 0) exit(1);
				}
				break;

			case DISCRETE_MODEL:
				// TODO(dimitris): do this through the interface and remove header
				dynamics = nlp_in->dynamics[i];
				dynamics->discrete_model = (external_function_generic *) &erk4_casadi[i];
				break;
		}
	}


    /* constraints */
	ocp_nlp_constraints_bgh_model **constraints = (ocp_nlp_constraints_bgh_model **) nlp_in->constraints;
	ocp_nlp_constraints_bgh_dims **constraints_dims = (ocp_nlp_constraints_bgh_dims **) dims->constraints;

	// fist stage
#if CONSTRAINTS==0 // box constraints
	// blasfeo_pack_dvec(nb[0], lb0, &constraints[0]->d, 0);
	// blasfeo_pack_dvec(nb[0], ub0, &constraints[0]->d, nb[0]+ng[0]);
	nlp_bounds_bgh_set(constraints_dims[0], constraints[0], "lb", lb0);
	nlp_bounds_bgh_set(constraints_dims[0], constraints[0], "ub", ub0);
    constraints[0]->idxb = idxb0;
#elif CONSTRAINTS==1 // general constraints
	double *Cu0; d_zeros(&Cu0, ng[0], nu[0]);
	for (int ii=0; ii<nu[0]; ii++)
		Cu0[ii*(ng[0]+1)] = 1.0;

	double *Cx0; d_zeros(&Cx0, ng[0], nx[0]);
	for (int ii=0; ii<nx[0]; ii++)
		Cx0[nu[0]+ii*(ng[0]+1)] = 1.0;

	blasfeo_pack_tran_dmat(ng[0], nu[0], Cu0, ng[0], &constraints[0]->DCt, 0, 0);
	blasfeo_pack_tran_dmat(ng[0], nx[0], Cx0, ng[0], &constraints[0]->DCt, nu[0], 0);
	// blasfeo_pack_dvec(ng[0], lb0, &constraints[0]->d, nb[0]);
	// blasfeo_pack_dvec(ng[0], ub0, &constraints[0]->d, 2*nb[0]+ng[0]);
	nlp_bounds_bgh_set(constraints_dims[0], constraints[0], "lg", lb0);
	nlp_bounds_bgh_set(constraints_dims[0], constraints[0], "ug", ub0);

	d_free(Cu0);
	d_free(Cx0);
#else // general + nonlinear constraints
	blasfeo_dgese(nu[0]+nx[0], ng[0], 0.0, &constraints[0]->DCt, 0, 0);
	for (int ii=0; ii<ng[0]; ii++)
		BLASFEO_DMATEL(&constraints[0]->DCt, ii, ii) = 1.0;

    ocp_nlp_constraints_bgh_model **nl_constr = (ocp_nlp_constraints_bgh_model **) nlp_in->constraints;
	nl_constr[0]->h = &nonlin_constr_generic;

	blasfeo_pack_dvec(ng[0]+nh[0], lb0, &constraints[0]->d, nb[0]);
	blasfeo_pack_dvec(ng[0]+nh[0], ub0, &constraints[0]->d, 2*nb[0]+ng[0]+nh[0]);
#endif

	// other stages
    for (int i = 1; i < NN; i++)
	{
		// blasfeo_pack_dvec(nb[i], lb1, &constraints[i]->d, 0);
		// blasfeo_pack_dvec(nb[i], ub1, &constraints[i]->d, nb[i]+ng[i]);
		nlp_bounds_bgh_set(constraints_dims[i], constraints[i], "lb", lb1);
		nlp_bounds_bgh_set(constraints_dims[i], constraints[i], "ub", ub1);
        constraints[i]->idxb = idxb1;
    }
	// blasfeo_pack_dvec(nb[NN], lbN, &constraints[NN]->d, 0);
	// blasfeo_pack_dvec(nb[NN], ubN, &constraints[NN]->d, nb[NN]+ng[NN]);
	nlp_bounds_bgh_set(constraints_dims[NN], constraints[NN], "lb", lbN);
	nlp_bounds_bgh_set(constraints_dims[NN], constraints[NN], "ub", ubN);
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
		if (plan->nlp_dynamics[i] == CONTINUOUS_MODEL)
		{
			ocp_nlp_dynamics_cont_opts *dynamics_stage_opts = sqp_opts->dynamics[i];
			sim_rk_opts *sim_opts = dynamics_stage_opts->sim_solver;

			if (plan->sim_solver_plan[i].sim_solver == ERK)
			{
				sim_opts->ns = 4;
			}
			else if (plan->sim_solver_plan[i].sim_solver == IRK)
			{
				sim_opts->ns = 2;
				sim_opts->jac_reuse = true;
			}
		}
    }

    sqp_opts->maxIter = MAX_SQP_ITERS;
    sqp_opts->min_res_g = 1e-9;
    sqp_opts->min_res_b = 1e-9;
    sqp_opts->min_res_d = 1e-9;
    sqp_opts->min_res_m = 1e-9;

    /************************************************
    * ocp_nlp out
    ************************************************/

	ocp_nlp_out *nlp_out = ocp_nlp_out_create(config, dims);

	ocp_nlp_solver *solver = ocp_nlp_create(config, dims, nlp_opts);

    /************************************************
    * sqp solve
    ************************************************/

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

		// call nlp solver
        status = ocp_nlp_solve(solver, nlp_in, nlp_out);
    }

    double time = acados_toc(&timer)/NREP;

	printf("\nresiduals\n");
	ocp_nlp_res_print(dims, ((ocp_nlp_sqp_memory *)solver->mem)->nlp_res);

	printf("\nsolution\n");
	ocp_nlp_out_print(dims, nlp_out);

	ocp_nlp_sqp_memory *solver_mem = (ocp_nlp_sqp_memory *) solver->mem;

	int sqp_iter = solver_mem->sqp_iter;

    printf("\n\nstatus = %i, iterations (max %d) = %d, total time = %f ms\n", status, MAX_SQP_ITERS, sqp_iter, time*1e3);
	printf("\nlinearization time = %f ms\n", solver_mem->time_lin*1e3);
	printf("\nqp solution time   = %f ms\n", solver_mem->time_qp_sol*1e3);
	printf("\ntotal time         = %f ms\n", solver_mem->time_tot*1e3);
	printf("\n\n");

    for (int k =0; k < 3; k++) {
        printf("u[%d] = \n", k);
		blasfeo_print_tran_dvec(nu[k], nlp_out->ux+k, 0);
        printf("x[%d] = \n", k);
		blasfeo_print_tran_dvec(nx[k], nlp_out->ux+k, nu[k]);
    }
    printf("u[N-1] = \n");
	blasfeo_print_tran_dvec(nu[NN-1], nlp_out->ux+NN-1, 0);
    printf("x[N] = \n");
	blasfeo_print_tran_dvec(nx[NN], nlp_out->ux+NN, nu[NN]);

    /************************************************
    * free memory
    ************************************************/

	// TODO(dimitris): VALGRIND!
 	external_function_casadi_free(expl_vde_for);
	free(expl_vde_for);

	external_function_casadi_free(impl_ode_fun);
	external_function_casadi_free(impl_ode_fun_jac_x_xdot);
	external_function_casadi_free(impl_ode_fun_jac_x_xdot_u);
	external_function_casadi_free(impl_ode_jac_x_xdot_u);

	free(impl_ode_fun);
	free(impl_ode_fun_jac_x_xdot);
	free(impl_ode_fun_jac_x_xdot_u);
	free(impl_ode_jac_x_xdot_u);

	external_function_casadi_free(erk4_casadi);
	free(erk4_casadi);

	free(nlp_opts);
	free(nlp_in);
	free(nlp_out);
	free(solver);
	free(dims);
	free(config);

	free(xref);
	free(diag_cost_x);
	free(lb0);
	free(ub0);
	free(lb1);
	free(ub1);
	free(lbN);
	free(ubN);
	free(idxb0);
	free(idxb1);
	free(idxbN);

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


	free(ls_cost_jac_casadi);
	free(ext_cost_generic);

	free(plan);

	/************************************************
	* return
	************************************************/

	if (status == 0)
		printf("\nsuccess! (%d iter) \n\n", sqp_iter);
	else
		printf("\nfailure!\n\n");

	return 0;
}
