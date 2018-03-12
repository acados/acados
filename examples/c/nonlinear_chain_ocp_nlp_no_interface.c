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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// #include <xmmintrin.h>

#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_partial_condensing_solver.h"

#include "acados/sim/sim_common.h"
#include "acados/sim/sim_erk_integrator.h"
#include "acados/sim/sim_irk_integrator.h"
#include "acados/sim/sim_lifted_irk_integrator.h"
#include "acados/utils/mem.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"
#include "acados/utils/external_function_generic.h"

#include "acados/ocp_nlp/ocp_nlp_sqp.h"
#include "acados/ocp_nlp/ocp_nlp_cost_common.h"
#include "acados/ocp_nlp/ocp_nlp_cost_ls.h"
#include "acados/ocp_nlp/ocp_nlp_cost_nls.h"

#include "examples/c/chain_model/chain_model.h"
#include "examples/c/implicit_chain_model/chain_model_impl.h"



// temp
#include "acados/ocp_qp/ocp_qp_hpipm.h"




#define NN 15
#define TF 3.0
#define Ns 2
#define MAX_SQP_ITERS 10
#define NREP 10


// process box constraints as general constraints
#define BC_AS_GC



// dynamics: 0 erk, 1 lifted_irk, 2 irk
#define DYNAMICS 2

// cost: 0 ls, 1 nls
#define COST 1



enum sensitivities_scheme {
    EXACT_NEWTON,
    INEXACT_NEWTON,
    INIS,
    FROZEN_INEXACT_NEWTON,
    FROZEN_INIS
};



static void print_problem_info(enum sensitivities_scheme sensitivities_type,
                               const int num_free_masses, const int num_stages)
{
    char scheme_name[MAX_STR_LEN];
    switch (sensitivities_type) {
        case EXACT_NEWTON:
            snprintf(scheme_name, sizeof(scheme_name), "EXACT_NEWTON");
            break;
        case INEXACT_NEWTON:
            snprintf(scheme_name, sizeof(scheme_name), "INEXACT_NEWTON");
            break;
        case INIS:
            snprintf(scheme_name, sizeof(scheme_name), "INIS");
            break;
        case FROZEN_INEXACT_NEWTON:
            snprintf(scheme_name, sizeof(scheme_name), "FROZEN_INEXACT_NEWTON");
            break;
        case FROZEN_INIS:
            snprintf(scheme_name, sizeof(scheme_name), "FROZEN_INIS");
            break;
        default:
            printf("Chose sensitivities type not available");
            exit(1);
    }
    printf("\n----- NUMBER OF FREE MASSES = %d, stages = %d (%s) -----\n",
           num_free_masses, num_stages, scheme_name);
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



static void select_dynamics_casadi(int N, int num_free_masses, external_function_casadi *forw_vde, external_function_casadi *jac_ode, external_function_casadi *impl_ode, external_function_casadi *impl_jac_x, external_function_casadi *impl_jac_xdot, external_function_casadi *impl_jac_u)
{
	// loop index
	int ii;

	switch (num_free_masses)
	{
		case 1:
			for (ii = 0; ii < N; ii++)
			{
#if DYNAMICS==0 | DYNAMICS==1
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
#else
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
#endif
			}
			break;
		case 2:
			for (ii = 0; ii < N; ii++)
			{
#if DYNAMICS==0 | DYNAMICS==1
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
#else
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
#endif
			}
			break;
		case 3:
			for (ii = 0; ii < N; ii++)
			{
#if DYNAMICS==0 | DYNAMICS==1
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
#else
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
#endif
			}
			break;
		case 4:
			for (ii = 0; ii < N; ii++)
			{
#if DYNAMICS==0 | DYNAMICS==1
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
#else
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
#endif
			}
		case 5:
			for (ii = 0; ii < N; ii++)
			{
#if DYNAMICS==0 | DYNAMICS==1
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
#else
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
#endif
			}
			break;
		default:
			printf("Problem size not available\n");
			exit(1);
			break;
	}

	return;
}



static void select_ls_cost_jac_casadi(int N, int num_free_masses, external_function_casadi *ls_cost_jac)
{
	// loop index
	int ii;

	switch (num_free_masses)
	{
		case 1:
			for (ii = 0; ii < N; ii++)
			{
				ls_cost_jac[ii].casadi_fun = &ls_cost_nm2;
				ls_cost_jac[ii].casadi_work = &ls_cost_nm2_work;
				ls_cost_jac[ii].casadi_sparsity_in = &ls_cost_nm2_sparsity_in;
				ls_cost_jac[ii].casadi_sparsity_out = &ls_cost_nm2_sparsity_out;
				ls_cost_jac[ii].casadi_n_in = &ls_cost_nm2_n_in;
				ls_cost_jac[ii].casadi_n_out = &ls_cost_nm2_n_out;
			}
			ls_cost_jac[N].casadi_fun = &ls_costN_nm2;
			ls_cost_jac[N].casadi_work = &ls_costN_nm2_work;
			ls_cost_jac[N].casadi_sparsity_in = &ls_costN_nm2_sparsity_in;
			ls_cost_jac[N].casadi_sparsity_out = &ls_costN_nm2_sparsity_out;
			ls_cost_jac[N].casadi_n_in = &ls_costN_nm2_n_in;
			ls_cost_jac[N].casadi_n_out = &ls_costN_nm2_n_out;
			break;
		case 2:
			for (ii = 0; ii < N; ii++)
			{
				ls_cost_jac[ii].casadi_fun = &ls_cost_nm3;
				ls_cost_jac[ii].casadi_work = &ls_cost_nm3_work;
				ls_cost_jac[ii].casadi_sparsity_in = &ls_cost_nm3_sparsity_in;
				ls_cost_jac[ii].casadi_sparsity_out = &ls_cost_nm3_sparsity_out;
				ls_cost_jac[ii].casadi_n_in = &ls_cost_nm3_n_in;
				ls_cost_jac[ii].casadi_n_out = &ls_cost_nm3_n_out;
			}
			ls_cost_jac[N].casadi_fun = &ls_costN_nm3;
			ls_cost_jac[N].casadi_work = &ls_costN_nm3_work;
			ls_cost_jac[N].casadi_sparsity_in = &ls_costN_nm3_sparsity_in;
			ls_cost_jac[N].casadi_sparsity_out = &ls_costN_nm3_sparsity_out;
			ls_cost_jac[N].casadi_n_in = &ls_costN_nm3_n_in;
			ls_cost_jac[N].casadi_n_out = &ls_costN_nm3_n_out;
			break;
		case 3:
			for (ii = 0; ii < N; ii++)
			{
				ls_cost_jac[ii].casadi_fun = &ls_cost_nm4;
				ls_cost_jac[ii].casadi_work = &ls_cost_nm4_work;
				ls_cost_jac[ii].casadi_sparsity_in = &ls_cost_nm4_sparsity_in;
				ls_cost_jac[ii].casadi_sparsity_out = &ls_cost_nm4_sparsity_out;
				ls_cost_jac[ii].casadi_n_in = &ls_cost_nm4_n_in;
				ls_cost_jac[ii].casadi_n_out = &ls_cost_nm4_n_out;
			}
			ls_cost_jac[N].casadi_fun = &ls_costN_nm4;
			ls_cost_jac[N].casadi_work = &ls_costN_nm4_work;
			ls_cost_jac[N].casadi_sparsity_in = &ls_costN_nm4_sparsity_in;
			ls_cost_jac[N].casadi_sparsity_out = &ls_costN_nm4_sparsity_out;
			ls_cost_jac[N].casadi_n_in = &ls_costN_nm4_n_in;
			ls_cost_jac[N].casadi_n_out = &ls_costN_nm4_n_out;
			break;
		case 4:
			for (ii = 0; ii < N; ii++)
			{
				ls_cost_jac[ii].casadi_fun = &ls_cost_nm5;
				ls_cost_jac[ii].casadi_work = &ls_cost_nm5_work;
				ls_cost_jac[ii].casadi_sparsity_in = &ls_cost_nm5_sparsity_in;
				ls_cost_jac[ii].casadi_sparsity_out = &ls_cost_nm5_sparsity_out;
				ls_cost_jac[ii].casadi_n_in = &ls_cost_nm5_n_in;
				ls_cost_jac[ii].casadi_n_out = &ls_cost_nm5_n_out;
			}
			ls_cost_jac[N].casadi_fun = &ls_costN_nm5;
			ls_cost_jac[N].casadi_work = &ls_costN_nm5_work;
			ls_cost_jac[N].casadi_sparsity_in = &ls_costN_nm5_sparsity_in;
			ls_cost_jac[N].casadi_sparsity_out = &ls_costN_nm5_sparsity_out;
			ls_cost_jac[N].casadi_n_in = &ls_costN_nm5_n_in;
			ls_cost_jac[N].casadi_n_out = &ls_costN_nm5_n_out;
			break;
		case 5:
			for (ii = 0; ii < N; ii++)
			{
				ls_cost_jac[ii].casadi_fun = &ls_cost_nm6;
				ls_cost_jac[ii].casadi_work = &ls_cost_nm6_work;
				ls_cost_jac[ii].casadi_sparsity_in = &ls_cost_nm6_sparsity_in;
				ls_cost_jac[ii].casadi_sparsity_out = &ls_cost_nm6_sparsity_out;
				ls_cost_jac[ii].casadi_n_in = &ls_cost_nm6_n_in;
				ls_cost_jac[ii].casadi_n_out = &ls_cost_nm6_n_out;
			}
			ls_cost_jac[N].casadi_fun = &ls_costN_nm6;
			ls_cost_jac[N].casadi_work = &ls_costN_nm6_work;
			ls_cost_jac[N].casadi_sparsity_in = &ls_costN_nm6_sparsity_in;
			ls_cost_jac[N].casadi_sparsity_out = &ls_costN_nm6_sparsity_out;
			ls_cost_jac[N].casadi_n_in = &ls_costN_nm6_n_in;
			ls_cost_jac[N].casadi_n_out = &ls_costN_nm6_n_out;
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
    FILE *initial_states_file;
    switch (num_free_masses)
    {
        case 1:
            initial_states_file = fopen(X0_NM2_FILE, "r");
            break;
        case 2:
            initial_states_file = fopen(X0_NM3_FILE, "r");
            break;
        case 3:
            initial_states_file = fopen(X0_NM4_FILE, "r");
            break;
        case 4:
            initial_states_file = fopen(X0_NM5_FILE, "r");
            break;
        case 5:
            initial_states_file = fopen(X0_NM6_FILE, "r");
            break;
        // case 6:
        //     initial_states_file = fopen(X0_NM7_FILE, "r");
        //     break;
        // case 7:
        //     initial_states_file = fopen(X0_NM8_FILE, "r");
        //     break;
        default:
            initial_states_file = fopen(X0_NM2_FILE, "r");
            break;
    }
    for (int i = 0; i < nx; i++)
        if (!fscanf(initial_states_file, "%lf", &x0[i]))
            break;
    fclose(initial_states_file);
}



void read_final_state(const int nx, const int num_free_masses, double *xN)
{
    FILE *final_state_file;
    switch (num_free_masses) {
        case 1:
            final_state_file = fopen(XN_NM2_FILE, "r");
            break;
        case 2:
            final_state_file = fopen(XN_NM3_FILE, "r");
            break;
        case 3:
            final_state_file = fopen(XN_NM4_FILE, "r");
            break;
        case 4:
            final_state_file = fopen(XN_NM5_FILE, "r");
            break;
        case 5:
            final_state_file = fopen(XN_NM6_FILE, "r");
            break;
        // case 6:
        //     final_state_file = fopen(XN_NM7_FILE, "r");
        //     break;
        // case 7:
        //     final_state_file = fopen(XN_NM8_FILE, "r");
        //     break;
        default:
            final_state_file = fopen(XN_NM2_FILE, "r");
            break;
    }
    for (int i = 0; i < nx; i++)
        if (!fscanf(final_state_file, "%lf", &xN[i]))
            break;
    fclose(final_state_file);
}



int main() {
    // _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);

    enum sensitivities_scheme scheme = EXACT_NEWTON;
    const int NMF = 3;  // number of free masses
    const int d = 0;  // number of stages in integrator

    print_problem_info(scheme, NMF, d);

    // dimensions
    int NX = 6 * NMF;
    int NU = 3;

    int nx[NN + 1] = {0};
    int nu[NN + 1] = {0};
    int nbx[NN + 1] = {0};
    int nbu[NN + 1] = {0};
    int nb[NN + 1] = {0};
    int ng[NN + 1] = {0};
    int nh[NN + 1] = {0};
    int ns[NN+1] = {0};
	int ny[NN+1] = {0};

    nx[0] = NX;
    nu[0] = NU;
#ifdef BC_AS_GC
    nbx[0] = 0;
    nbu[0] = 0;
	nb[0] = 0;
    ng[0] = nu[0]+nx[0];
#else
    nbx[0] = nx[0];
    nbu[0] = nu[0];
    nb[0] = nbu[0]+nbx[0];
	ng[0] = 0;
#endif
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
    * config
    ************************************************/

	int config_size = ocp_nlp_solver_config_calculate_size(NN);
	void *config_mem = malloc(config_size);
	ocp_nlp_solver_config *config = ocp_nlp_solver_config_assign(NN, config_mem);

#if 1
	// partial condensing HPIPM
	ocp_qp_partial_condensing_solver_config_initialize_default(config->qp_solver);
	ocp_qp_hpipm_config_initialize_default(config->qp_solver->qp_solver);
#else
	// full condensing HPIPM
	ocp_qp_full_condensing_solver_config_initialize_default(config->qp_solver);
	dense_qp_hpipm_config_initialize_default(config->qp_solver->qp_solver);
#endif


	// cost: least squares
    for (int ii = 0; ii <= NN; ii++)
    {
#if COST==0
		ocp_nlp_cost_ls_config_initialize_default(config->cost[ii]);
#else
		ocp_nlp_cost_nls_config_initialize_default(config->cost[ii]);
#endif
    }

		// 4th order schemes
#if DYNAMICS==0
	// dynamics: ERK 4
    for (int ii = 0; ii < NN; ii++)
    {
		ocp_nlp_dynamics_config_initialize_default(config->dynamics[ii]);
		sim_erk_config_initialize_default(config->dynamics[ii]->sim_solver);
		config->dynamics[ii]->sim_solver->ns = 4; // number of integration stages
    }

#elif DYNAMICS==1
	// dynamics: lifted IRK GL2
    for (int ii = 0; ii < NN; ii++)
    {
		ocp_nlp_dynamics_config_initialize_default(config->dynamics[ii]);
		sim_lifted_irk_config_initialize_default(config->dynamics[ii]->sim_solver);
		config->dynamics[ii]->sim_solver->ns = 2; // number of integration stages
    }
#else
	// dynamics: IRK GL2
    for (int ii = 0; ii < NN; ii++)
    {
		ocp_nlp_dynamics_config_initialize_default(config->dynamics[ii]);
		sim_irk_config_initialize_default(config->dynamics[ii]->sim_solver);
		config->dynamics[ii]->sim_solver->ns = 2; // number of integration stages
    }
#endif

	// constraitns
    for (int ii = 0; ii <= NN; ii++)
    {
		ocp_nlp_constraints_linear_config_initialize_default(config->constraints[ii]);
    }

    /************************************************
    * ocp_nlp_dims
    ************************************************/

	int dims_size = ocp_nlp_dims_calculate_size(NN);
	void *dims_mem = malloc(dims_size);
	ocp_nlp_dims *dims = ocp_nlp_dims_assign(NN, dims_mem);
	ocp_nlp_dims_initialize(nx, nu, ny, nbx, nbu, ng, nh, ns, dims);

    /************************************************
    * dynamics
    ************************************************/

	// explicit
	external_function_casadi forw_vde_casadi[NN]; // XXX varible size array
	external_function_casadi jac_ode_casadi[NN]; // XXX varible size array
	// implicit
	external_function_casadi impl_ode_casadi[NN]; // XXX varible size array
	external_function_casadi impl_jac_x_casadi[NN]; // XXX varible size array
	external_function_casadi impl_jac_xdot_casadi[NN]; // XXX varible size array
	external_function_casadi impl_jac_u_casadi[NN]; // XXX varible size array

	select_dynamics_casadi(NN, NMF, forw_vde_casadi, jac_ode_casadi, impl_ode_casadi, impl_jac_x_casadi, impl_jac_xdot_casadi, impl_jac_u_casadi);

	int tmp_size;
	char *c_ptr;
#if DYNAMICS==0 | DYNAMICS==1

	// forw_vde
	tmp_size = 0;
	for (int ii=0; ii<NN; ii++)
	{
		tmp_size += external_function_casadi_calculate_size(forw_vde_casadi+ii);
	}
	void *forw_vde_casadi_mem = malloc(tmp_size);
	c_ptr = forw_vde_casadi_mem;
	for (int ii=0; ii<NN; ii++)
	{
		external_function_casadi_assign(forw_vde_casadi+ii, c_ptr);
		c_ptr += external_function_casadi_calculate_size(forw_vde_casadi+ii);
	}
	// jac_ode
	tmp_size = 0;
	for (int ii=0; ii<NN; ii++)
	{
		tmp_size += external_function_casadi_calculate_size(jac_ode_casadi+ii);
	}
	void *jac_ode_casadi_mem = malloc(tmp_size);
	c_ptr = jac_ode_casadi_mem;
	for (int ii=0; ii<NN; ii++)
	{
		external_function_casadi_assign(jac_ode_casadi+ii, c_ptr);
		c_ptr += external_function_casadi_calculate_size(jac_ode_casadi+ii);
	}

#else // DYNAMICS==2

	// impl_ode
	tmp_size = 0;
	for (int ii=0; ii<NN; ii++)
	{
		tmp_size += external_function_casadi_calculate_size(impl_ode_casadi+ii);
	}
	void *impl_ode_casadi_mem = malloc(tmp_size);
	c_ptr = impl_ode_casadi_mem;
	for (int ii=0; ii<NN; ii++)
	{
		external_function_casadi_assign(impl_ode_casadi+ii, c_ptr);
		c_ptr += external_function_casadi_calculate_size(impl_ode_casadi+ii);
	}
	// jac_x
	tmp_size = 0;
	for (int ii=0; ii<NN; ii++)
	{
		tmp_size += external_function_casadi_calculate_size(impl_jac_x_casadi+ii);
	}
	void *impl_jac_x_casadi_mem = malloc(tmp_size);
	c_ptr = impl_jac_x_casadi_mem;
	for (int ii=0; ii<NN; ii++)
	{
		external_function_casadi_assign(impl_jac_x_casadi+ii, c_ptr);
		c_ptr += external_function_casadi_calculate_size(impl_jac_x_casadi+ii);
	}
	// jac_xdot
	tmp_size = 0;
	for (int ii=0; ii<NN; ii++)
	{
		tmp_size += external_function_casadi_calculate_size(impl_jac_xdot_casadi+ii);
	}
	void *impl_jac_xdot_casadi_mem = malloc(tmp_size);
	c_ptr = impl_jac_xdot_casadi_mem;
	for (int ii=0; ii<NN; ii++)
	{
		external_function_casadi_assign(impl_jac_xdot_casadi+ii, c_ptr);
		c_ptr += external_function_casadi_calculate_size(impl_jac_xdot_casadi+ii);
	}
	// jac_u
	tmp_size = 0;
	for (int ii=0; ii<NN; ii++)
	{
		tmp_size += external_function_casadi_calculate_size(impl_jac_u_casadi+ii);
	}
	void *impl_jac_u_casadi_mem = malloc(tmp_size);
	c_ptr = impl_jac_u_casadi_mem;
	for (int ii=0; ii<NN; ii++)
	{
		external_function_casadi_assign(impl_jac_u_casadi+ii, c_ptr);
		c_ptr += external_function_casadi_calculate_size(impl_jac_u_casadi+ii);
	}

#endif // DYNAMICS

    /************************************************
    * nonlinear least squares
    ************************************************/

#if COST==1
	external_function_casadi ls_cost_jac_casadi[NN+1]; // XXX varible size array

	select_ls_cost_jac_casadi(NN, NMF, ls_cost_jac_casadi);

	// ls_cost_jac
	tmp_size = 0;
	for (int ii=0; ii<=NN; ii++)
	{
		tmp_size += external_function_casadi_calculate_size(ls_cost_jac_casadi+ii);
	}
	void *ls_cost_jac_casadi_mem = malloc(tmp_size);
	c_ptr = ls_cost_jac_casadi_mem;
	for (int ii=0; ii<=NN; ii++)
	{
		external_function_casadi_assign(ls_cost_jac_casadi+ii, c_ptr);
		c_ptr += external_function_casadi_calculate_size(ls_cost_jac_casadi+ii);
	}
#endif

    /************************************************
    * nlp_in (wip)
    ************************************************/

    // TODO(dimitris): clean up integrators inside
	tmp_size = ocp_nlp_in_calculate_size(config, dims);
	void *nlp_in_mem = malloc(tmp_size);
	ocp_nlp_in *nlp_in = ocp_nlp_in_assign(config, dims, nlp_in_mem);


	// sampling times
	for (int ii=0; ii<NN; ii++)
		nlp_in->Ts[ii] = TF/NN;




//	ocp_nlp_dims_print(nlp_in->dims);

    // NOTE(dimitris): use nlp_in->dims instead of &dims from now on since nb is filled with nbx+nbu!

    // Problem data
    double wall_pos = -0.01;
    double UMAX = 10;

	double x_pos_inf = +1e4;
	double x_neg_inf = -1e4;

    double xref[NX];
    read_final_state(NX, NMF, xref);
    double uref[3] = {0.0, 0.0, 0.0};
    double diag_cost_x[NX];
    for (int i = 0; i < NX; i++)
        diag_cost_x[i] = 1e-2;
    double diag_cost_u[3] = {1.0, 1.0, 1.0};



    /* least-squares cost */

	// output definition: y = [x; u]

#if COST==0
    /* linear ls */

    ocp_nlp_cost_ls_model **cost_ls = (ocp_nlp_cost_ls_model **) nlp_in->cost;

	// Cyt
	for (int i=0; i<=NN; i++)
	{
		blasfeo_dgese(nv[i], ny[i], 0.0, &cost_ls[i]->Cyt, 0, 0);
        for (int j = 0; j < nu[i]; j++)
            DMATEL_LIBSTR(&cost_ls[i]->Cyt, j, nx[i]+j) = 1.0;
        for (int j = 0; j < nx[i]; j++)
            DMATEL_LIBSTR(&cost_ls[i]->Cyt, nu[i]+j, j) = 1.0;
	}

#else
    /* nonlinear ls */

    ocp_nlp_cost_nls_model **cost_ls = (ocp_nlp_cost_nls_model **) nlp_in->cost;

	// nls_jac
	for (int i=0; i<=NN; i++)
		cost_ls[i]->nls_jac = (external_function_generic *) &ls_cost_jac_casadi[i];
#if 0
	// replace with hand-written external functions
	external_function_generic ls_cost_jac_generic[NN];
	if (NMF==3)
	{
		for (int i=0; i<NN; i++)
		{
			ls_cost_jac_generic[i].evaluate = &ls_cost_jac_nm4;
			cost_ls->nls_jac[i] = &ls_cost_jac_generic[i];
		}
	}
#endif

#endif

    /* common to linear and nonlinear ls */

	// W
	for (int i=0; i<=NN; i++)
	{
		blasfeo_dgese(ny[i], ny[i], 0.0, &cost_ls[i]->W, 0, 0);
        for (int j = 0; j < nx[i]; j++)
            DMATEL_LIBSTR(&cost_ls[i]->W, j, j) = diag_cost_x[j];
        for (int j = 0; j < nu[i]; j++)
            DMATEL_LIBSTR(&cost_ls[i]->W, nx[i]+j, nx[i]+j) = diag_cost_u[j];
	}

	// y_ref
    for (int i = 0; i < NN; i++)
	{
		blasfeo_pack_dvec(nx[i], xref, &cost_ls[i]->y_ref, 0);
		blasfeo_pack_dvec(nu[i], uref, &cost_ls[i]->y_ref, nx[i]);
    }



	/* dynamics */
#if DYNAMICS==0
	for (int i=0; i<NN; i++)
	{
		ocp_nlp_dynamics_model *dynamics = nlp_in->dynamics[i];
		erk_model *model = dynamics->sim_model;
		model->forw_vde_expl = (external_function_generic *) &forw_vde_casadi[i];
		model->jac_ode_expl = (external_function_generic *) &jac_ode_casadi[i];
	}
#elif DYNAMICS==1
	for (int i=0; i<NN; i++)
	{
		ocp_nlp_dynamics_model *dynamics = nlp_in->dynamics[i];
		lifted_irk_model *model = dynamics->sim_model;
		model->forw_vde_expl = (external_function_generic *) &forw_vde_casadi[i];
		model->jac_ode_expl = (external_function_generic *) &jac_ode_casadi[i];
	}
#else
	for (int i=0; i<NN; i++)
	{
		ocp_nlp_dynamics_model *dynamics = nlp_in->dynamics[i];
		irk_model *model = dynamics->sim_model;
		model->ode_impl = (external_function_generic *) &impl_ode_casadi[i];
		model->jac_x_ode_impl = (external_function_generic *) &impl_jac_x_casadi[i];
		model->jac_xdot_ode_impl = (external_function_generic *) &impl_jac_xdot_casadi[i];
		model->jac_u_ode_impl = (external_function_generic *) &impl_jac_u_casadi[i];
	}
#endif



    nlp_in->freezeSens = false;
    if (scheme > 2)
        nlp_in->freezeSens = true;



    /* box constraints */

	ocp_nlp_constraints_linear_model **constraints = (ocp_nlp_constraints_linear_model **) nlp_in->constraints;

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

	// stage-wise
	blasfeo_pack_dvec(nb[0], lb0, &constraints[0]->d, 0);
	blasfeo_pack_dvec(nb[0], ub0, &constraints[0]->d, nb[0]+ng[0]);
    constraints[0]->idxb = idxb0;
    for (int i = 1; i < NN; i++)
	{
		blasfeo_pack_dvec(nb[i], lb1, &constraints[i]->d, 0);
		blasfeo_pack_dvec(nb[i], ub1, &constraints[i]->d, nb[i]+ng[i]);
        constraints[i]->idxb = idxb1;
    }
	blasfeo_pack_dvec(nb[NN], lbN, &constraints[NN]->d, 0);
	blasfeo_pack_dvec(nb[NN], ubN, &constraints[NN]->d, nb[NN]+ng[NN]);
    constraints[NN]->idxb = idxbN;


	// General constraints
	if (ng[0]>0)
	{
		double *Cu0; d_zeros(&Cu0, ng[0], nu[0]);
		for (int ii=0; ii<nu[0]; ii++)
			Cu0[ii*(ng[0]+1)] = 1.0;

		double *Cx0; d_zeros(&Cx0, ng[0], nx[0]);
		for (int ii=0; ii<nx[0]; ii++)
			Cx0[nu[0]+ii*(ng[0]+1)] = 1.0;

		blasfeo_pack_tran_dmat(ng[0], nu[0], Cu0, ng[0], &constraints[0]->DCt, 0, 0);
		blasfeo_pack_tran_dmat(ng[0], nx[0], Cx0, ng[0], &constraints[0]->DCt, nu[0], 0);
		blasfeo_pack_dvec(ng[0], lb0, &constraints[0]->d, nb[0]);
		blasfeo_pack_dvec(ng[0], ub0, &constraints[0]->d, 2*nb[0]+ng[0]);

		d_free(Cu0);
		d_free(Cx0);
	}
#if 0
	blasfeo_print_dmat(nu[0]+nx[0], ng[0], constraints->DCt+0, 0, 0);
	blasfeo_print_tran_dvec(2*nb[0]+2*ng[0], constraints->d+0, 0);
//	exit(1);
#endif

    /************************************************
    * sqp opts
    ************************************************/

	tmp_size = ocp_nlp_sqp_opts_calculate_size(config, dims);
	void *nlp_opts_mem = malloc(tmp_size);
	ocp_nlp_sqp_opts *nlp_opts = ocp_nlp_sqp_opts_assign(config, dims, nlp_opts_mem);

	ocp_nlp_sqp_opts_initialize_default(config, dims, nlp_opts);

    for (int i = 0; i < NN; ++i)
	{
		ocp_nlp_dynamics_opts *dynamics_opts = nlp_opts->dynamics[i];
        sim_rk_opts *sim_opts = dynamics_opts->sim_solver;
#if DYNAMICS==2
		sim_opts->jac_reuse = true;
#endif
    }


	// XXX hack: overwrite config with hand-setted one
//	nlp_opts->qp_solver = &config_qp;
//	nlp_opts->sim_solvers = config_sim_ptrs;
//	for (int ii=0; ii<NN; ii++)
//		nlp_opts->sim_solvers[ii] = config_sim_ptrs[ii];



    nlp_opts->maxIter = MAX_SQP_ITERS;
    nlp_opts->min_res_g = 1e-9;
    nlp_opts->min_res_b = 1e-9;
    nlp_opts->min_res_d = 1e-9;
    nlp_opts->min_res_m = 1e-9;

    /************************************************
    * ocp_nlp out
    ************************************************/

	tmp_size = ocp_nlp_out_calculate_size(config, dims);
	void *nlp_out_mem = malloc(tmp_size);
	ocp_nlp_out *nlp_out = ocp_nlp_out_assign(config, dims, nlp_out_mem);

//	ocp_nlp_dims_print(nlp_out->dims);

    /************************************************
    * sqp memory
    ************************************************/

	tmp_size = ocp_nlp_sqp_memory_calculate_size(config, dims, nlp_opts);
	void *nlp_mem_mem = malloc(tmp_size);
	ocp_nlp_sqp_memory *nlp_mem = ocp_nlp_sqp_memory_assign(config, dims, nlp_opts, nlp_mem_mem);


    /************************************************
    * sqp workspace
    ************************************************/

    int workspace_size = ocp_nlp_sqp_workspace_calculate_size(config, nlp_in->dims, nlp_opts);
    void *nlp_work = acados_malloc(workspace_size, 1);

    /************************************************
    * sqp solve
    ************************************************/

    int status;

    acados_timer timer;
    acados_tic(&timer);

    for (int rep = 0; rep < NREP; rep++)
    {
		// warm start output initial guess of solution
//		if (rep==0)
//		{
			for (int i=0; i<=NN; i++)
			{
				blasfeo_pack_dvec(nu[i], uref, nlp_out->ux+i, 0);
				blasfeo_pack_dvec(nx[i], xref, nlp_out->ux+i, nu[i]);
			}
//		}

		// call nlp solver
        status = ocp_nlp_sqp(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);
    }

    double time = acados_toc(&timer)/NREP;

	printf("\nresiduals\n");
	ocp_nlp_res_print(nlp_mem->nlp_res);

	printf("\nsolution\n");
	ocp_nlp_out_print(nlp_out);

    printf("\n\nstatus = %i, iterations (max %d) = %d, total time = %f ms\n\n", status, MAX_SQP_ITERS, nlp_mem->sqp_iter, time*1e3);

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

#if DYNAMICS==0 | DYNAMICS==1

	free(forw_vde_casadi_mem);
	free(jac_ode_casadi_mem);

#else // DYNAMICS==2

	free(impl_ode_casadi_mem);
	free(impl_jac_x_casadi_mem);
	free(impl_jac_xdot_casadi_mem);
	free(impl_jac_u_casadi_mem);

#endif // DYNAMICS

	free(config_mem);
#if COST==1
	free(ls_cost_jac_casadi_mem);
#endif
	free(dims_mem);
    free(nlp_in_mem);
    free(nlp_out_mem);
    free(nlp_work);
    free(nlp_mem_mem);
    free(nlp_opts_mem);

/************************************************
* return
************************************************/

	printf("\nsuccess!\n\n");

	return 0;

}
