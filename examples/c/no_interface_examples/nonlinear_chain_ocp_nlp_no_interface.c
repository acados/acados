/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"

#include "acados/dense_qp/dense_qp_hpipm.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_partial_condensing_solver.h"
#include "acados/ocp_qp/ocp_qp_full_condensing_solver.h"
#include "acados/sim/sim_common.h"
#include "acados/sim/sim_erk_integrator.h"
#include "acados/sim/sim_irk_integrator.h"
#include "acados/sim/sim_lifted_irk_integrator.h"
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



// temp
#include "acados/ocp_qp/ocp_qp_hpipm.h"

#define NN 15
#define TF 3.75

#define MAX_SQP_ITERS 10
#define NREP 1

#define NUM_FREE_MASSES 3
#define LIFTED 1

// dynamics: 0 erk, 1 old_lifted_irk, 2 irk, 3 discrete_model, 4 lifted_irk
#define DYNAMICS 4

// cost: 0 ls, 1 nls, 2 external
#define COST 2

// constraints (at stage 0): 0 box, 1 general, 2 general+nonlinear
#define CONSTRAINTS 2

// xcond: 0 no condensing, 1 part condensing, 2 full condensing
#define XCOND 1



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



static void select_dynamics_casadi(int N, int num_free_masses, 
	external_function_casadi *forw_vde, 
	external_function_casadi *jac_ode, 
	external_function_casadi *impl_ode_fun, 
	external_function_casadi *impl_ode_fun_jac_x_xdot,
	external_function_casadi *impl_ode_fun_jac_x_xdot_u,
	external_function_casadi *impl_ode_jac_x_xdot_u,
	external_function_casadi *erk4_casadi)
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
#elif DYNAMICS==2 | DYNAMICS==4
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
#elif DYNAMICS==3
				erk4_casadi[ii].casadi_fun = &casadi_erk4_chain_nm2;
				erk4_casadi[ii].casadi_work = &casadi_erk4_chain_nm2_work;
				erk4_casadi[ii].casadi_sparsity_in = &casadi_erk4_chain_nm2_sparsity_in;
				erk4_casadi[ii].casadi_sparsity_out = &casadi_erk4_chain_nm2_sparsity_out;
				erk4_casadi[ii].casadi_n_in = &casadi_erk4_chain_nm2_n_in;
				erk4_casadi[ii].casadi_n_out = &casadi_erk4_chain_nm2_n_out;
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
#elif DYNAMICS==2 | DYNAMICS==4
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
#elif DYNAMICS == 3
				erk4_casadi[ii].casadi_fun = &casadi_erk4_chain_nm3;
				erk4_casadi[ii].casadi_work = &casadi_erk4_chain_nm3_work;
				erk4_casadi[ii].casadi_sparsity_in = &casadi_erk4_chain_nm3_sparsity_in;
				erk4_casadi[ii].casadi_sparsity_out = &casadi_erk4_chain_nm3_sparsity_out;
				erk4_casadi[ii].casadi_n_in = &casadi_erk4_chain_nm3_n_in;
				erk4_casadi[ii].casadi_n_out = &casadi_erk4_chain_nm3_n_out;
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
#elif DYNAMICS==2 | DYNAMICS==4
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
#elif DYNAMICS == 3
				erk4_casadi[ii].casadi_fun = &casadi_erk4_chain_nm4;
				erk4_casadi[ii].casadi_work = &casadi_erk4_chain_nm4_work;
				erk4_casadi[ii].casadi_sparsity_in = &casadi_erk4_chain_nm4_sparsity_in;
				erk4_casadi[ii].casadi_sparsity_out = &casadi_erk4_chain_nm4_sparsity_out;
				erk4_casadi[ii].casadi_n_in = &casadi_erk4_chain_nm4_n_in;
				erk4_casadi[ii].casadi_n_out = &casadi_erk4_chain_nm4_n_out;
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
#elif DYNAMICS==2 | DYNAMICS==4
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
#elif DYNAMICS == 3
				erk4_casadi[ii].casadi_fun = &casadi_erk4_chain_nm5;
				erk4_casadi[ii].casadi_work = &casadi_erk4_chain_nm5_work;
				erk4_casadi[ii].casadi_sparsity_in = &casadi_erk4_chain_nm5_sparsity_in;
				erk4_casadi[ii].casadi_sparsity_out = &casadi_erk4_chain_nm5_sparsity_out;
				erk4_casadi[ii].casadi_n_in = &casadi_erk4_chain_nm5_n_in;
				erk4_casadi[ii].casadi_n_out = &casadi_erk4_chain_nm5_n_out;
#endif
			}
			break;
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
#elif DYNAMICS==2 | DYNAMICS==4
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
#elif DYNAMICS == 3
				erk4_casadi[ii].casadi_fun = &casadi_erk4_chain_nm6;
				erk4_casadi[ii].casadi_work = &casadi_erk4_chain_nm6_work;
				erk4_casadi[ii].casadi_sparsity_in = &casadi_erk4_chain_nm6_sparsity_in;
				erk4_casadi[ii].casadi_sparsity_out = &casadi_erk4_chain_nm6_sparsity_out;
				erk4_casadi[ii].casadi_n_in = &casadi_erk4_chain_nm6_n_in;
				erk4_casadi[ii].casadi_n_out = &casadi_erk4_chain_nm6_n_out;
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
#if 1
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
#else
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
#endif
}



void read_final_state(const int nx, const int num_free_masses, double *xN)
{
#if 1
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
#else
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
#endif
}



// hand-generated external function for externally provided hessian and gradient
void ext_cost_nm2(void *fun, ext_fun_arg_t *type_in, void **in, ext_fun_arg_t *type_out, void **out)
{

	int ii;

	int nu = 3;
	int nx = 6;

	int nv = nu+nx;

	// ref
	double ref[nu+nx];
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

	return;

}

void ext_cost_nm3(void *fun, ext_fun_arg_t *type_in, void **in, ext_fun_arg_t *type_out, void **out)
{

	int ii;

	int nu = 3;
	int nx = 12;

	int nv = nu+nx;

	// ref
	double ref[nu+nx];
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

	return;

}

void ext_cost_nm4(void *fun, ext_fun_arg_t *type_in, void **in, ext_fun_arg_t *type_out, void **out)
{

	int ii;

	int nu = 3;
	int nx = 18;

	int nv = nu+nx;

	// ref
	double ref[nu+nx];
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

	return;

}

void ext_cost_nm5(void *fun, ext_fun_arg_t *type_in, void **in, ext_fun_arg_t *type_out, void **out)
{

	int ii;

	int nu = 3;
	int nx = 24;

	int nv = nu+nx;

	// ref
	double ref[nu+nx];
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

	return;

}

void ext_cost_nm6(void *fun, ext_fun_arg_t *type_in, void **in, ext_fun_arg_t *type_out, void **out)
{

	int ii;

	int nu = 3;
	int nx = 30;

	int nv = nu+nx;

	// ref
	double ref[nu+nx];
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



// ls hessian
void ls_cost_hess_nm2(void *evaluate, ext_fun_arg_t *type_in, void **in, ext_fun_arg_t *type_out, void **out)
{

	int ii;

	int nu = 3;
	int nx = 6;

	struct blasfeo_dmat *hess = out[0];
	blasfeo_dgese(nu+nx, nu+nx, 0.0, hess, 0, 0);

	return;

}

void ls_cost_hess_nm3(void *evaluate, ext_fun_arg_t *type_in, void **in, ext_fun_arg_t *type_out, void **out)
{

	int ii;

	int nu = 3;
	int nx = 12;

	struct blasfeo_dmat *hess = out[0];
	blasfeo_dgese(nu+nx, nu+nx, 0.0, hess, 0, 0);

	return;

}

void ls_cost_hess_nm4(void *evaluate, ext_fun_arg_t *type_in, void **in, ext_fun_arg_t *type_out, void **out)
{

	int ii;

	int nu = 3;
	int nx = 18;

	struct blasfeo_dmat *hess = out[0];
	blasfeo_dgese(nu+nx, nu+nx, 0.0, hess, 0, 0);

	return;

}

void ls_cost_hess_nm5(void *evaluate, ext_fun_arg_t *type_in, void **in, ext_fun_arg_t *type_out, void **out)
{

	int ii;

	int nu = 3;
	int nx = 24;

	struct blasfeo_dmat *hess = out[0];
	blasfeo_dgese(nu+nx, nu+nx, 0.0, hess, 0, 0);

	return;

}

void ls_cost_hess_nm6(void *evaluate, ext_fun_arg_t *type_in, void **in, ext_fun_arg_t *type_out, void **out)
{

	int ii;

	int nu = 3;
	int nx = 30;

	struct blasfeo_dmat *hess = out[0];
	blasfeo_dgese(nu+nx, nu+nx, 0.0, hess, 0, 0);

	return;

}



/************************************************
* main
************************************************/

int main() {
    // _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);

    enum sensitivities_scheme scheme = EXACT_NEWTON;
    const int NMF = NUM_FREE_MASSES;  // number of free masses
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
    * config
    ************************************************/

	int config_size = ocp_nlp_config_calculate_size(NN);
	void *config_mem = malloc(config_size);
	ocp_nlp_config *config = ocp_nlp_config_assign(NN, config_mem);

#if XCOND==2
	// full condensing HPIPM
	ocp_qp_full_condensing_solver_config_initialize_default(config->qp_solver);
	dense_qp_hpipm_config_initialize_default(config->qp_solver->qp_solver);
#else
	// no condensing or partial condensing HPIPM
	ocp_qp_partial_condensing_solver_config_initialize_default(config->qp_solver);
	ocp_qp_hpipm_config_initialize_default(config->qp_solver->qp_solver);
#endif


	// cost: least squares
#if COST==0
    for (int ii = 0; ii <= NN; ii++)
    {
		// linear ls
		ocp_nlp_cost_ls_config_initialize_default(config->cost[ii]);
    }
#elif COST==1
    for (int ii = 0; ii <= NN; ii++)
    {
		// nonlinear ls
		ocp_nlp_cost_nls_config_initialize_default(config->cost[ii]);
    }
#else
    for (int ii = 0; ii < NN; ii++)
    {
		// external cost
		ocp_nlp_cost_external_config_initialize_default(config->cost[ii]);
    }
	// linear ls
	ocp_nlp_cost_ls_config_initialize_default(config->cost[NN]);
#endif

#if DYNAMICS==0
	// dynamics: ERK
    for (int ii = 0; ii < NN; ii++)
    {
		ocp_nlp_dynamics_cont_config_initialize_default(config->dynamics[ii]);
		sim_erk_config_initialize_default(config->dynamics[ii]->sim_solver);
    }

#elif DYNAMICS==1
	// dynamics: lifted IRK
    for (int ii = 0; ii < NN; ii++)
    {
		ocp_nlp_dynamics_cont_config_initialize_default(config->dynamics[ii]);
		sim_lifted_irk_config_initialize_default(config->dynamics[ii]->sim_solver);
    }
#elif DYNAMICS==2
	// dynamics: IRK
    for (int ii = 0; ii < NN; ii++)
    {
		ocp_nlp_dynamics_cont_config_initialize_default(config->dynamics[ii]);
		sim_irk_config_initialize_default(config->dynamics[ii]->sim_solver);
    }
#elif DYNAMICS==3
	// dynamics: discrete model
    for (int ii = 0; ii < NN; ii++)
    {
		ocp_nlp_dynamics_disc_config_initialize_default(config->dynamics[ii]);
    }
#elif DYNAMICS==4
	// dynamics: lifted IRK
    for (int ii = 0; ii < NN; ii++)
    {
		ocp_nlp_dynamics_cont_config_initialize_default(config->dynamics[ii]);
		sim_lifted_irk_config_initialize_default(config->dynamics[ii]->sim_solver);
    }
#endif

	// constraitns
    for (int ii = 0; ii <= NN; ii++)
    {
		ocp_nlp_constraints_bgh_config_initialize_default(config->constraints[ii]);
    }

    /************************************************
    * ocp_nlp_dims
    ************************************************/

	int dims_size = ocp_nlp_dims_calculate_size(config);
	void *dims_mem = malloc(dims_size);
	ocp_nlp_dims *dims = ocp_nlp_dims_assign(config, dims_mem);
	ocp_nlp_dims_set_opt_vars(config, dims, "nx", nx);
	ocp_nlp_dims_set_opt_vars(config, dims, "nu", nu);
	ocp_nlp_dims_set_opt_vars(config, dims, "nz", nz);
	ocp_nlp_dims_set_opt_vars(config, dims, "ns", ns);
	for(int ii=0; ii<=NN; ii++)
		{
		ocp_nlp_dims_set_constraints(config, dims, ii, "nbx", nbx);
		ocp_nlp_dims_set_constraints(config, dims, ii, "nbu", nbu);
		ocp_nlp_dims_set_constraints(config, dims, ii, "ng", ng);
		ocp_nlp_dims_set_constraints(config, dims, ii, "nh", nh);
		ocp_nlp_dims_set_cost(config, dims, ii, "ny", ny);
		}

    /************************************************
    * dynamics
    ************************************************/

	// explicit
	external_function_casadi *expl_vde_for = malloc(NN*sizeof(external_function_casadi));
	external_function_casadi *expl_ode_jac = malloc(NN*sizeof(external_function_casadi));

	// implicit
	external_function_casadi *impl_ode_fun = malloc(NN*sizeof(external_function_casadi));
	external_function_casadi *impl_ode_fun_jac_x_xdot = malloc(NN*sizeof(external_function_casadi));
	external_function_casadi *impl_ode_fun_jac_x_xdot_u = malloc(NN*sizeof(external_function_casadi));
	external_function_casadi *impl_ode_jac_x_xdot_u = malloc(NN*sizeof(external_function_casadi));

	// casadi erk
	external_function_casadi *erk4_casadi = malloc(NN*sizeof(external_function_casadi));

	select_dynamics_casadi(NN, NMF, expl_vde_for, expl_ode_jac, impl_ode_fun, impl_ode_fun_jac_x_xdot, impl_ode_fun_jac_x_xdot_u, impl_ode_jac_x_xdot_u, erk4_casadi);

	int tmp_size;
	char *c_ptr;
#if DYNAMICS==0 | DYNAMICS==1

	// forw_vde
	tmp_size = 0;
	for (int ii=0; ii<NN; ii++)
	{
		tmp_size += external_function_casadi_calculate_size(expl_vde_for+ii);
	}
	void *forw_vde_casadi_mem = malloc(tmp_size);
	c_ptr = forw_vde_casadi_mem;
	for (int ii=0; ii<NN; ii++)
	{
		external_function_casadi_assign(expl_vde_for+ii, c_ptr);
		c_ptr += external_function_casadi_calculate_size(expl_vde_for+ii);
	}
	// jac_ode
	tmp_size = 0;
	for (int ii=0; ii<NN; ii++)
	{
		tmp_size += external_function_casadi_calculate_size(expl_ode_jac+ii);
	}
	void *jac_ode_casadi_mem = malloc(tmp_size);
	c_ptr = jac_ode_casadi_mem;
	for (int ii=0; ii<NN; ii++)
	{
		external_function_casadi_assign(expl_ode_jac+ii, c_ptr);
		c_ptr += external_function_casadi_calculate_size(expl_ode_jac+ii);
	}

#elif DYNAMICS==2 |  DYNAMICS==4

	// impl_ode
	tmp_size = 0;
	for (int ii=0; ii<NN; ii++)
	{
		tmp_size += external_function_casadi_calculate_size(impl_ode_fun+ii);
	}
	void *impl_ode_casadi_mem = malloc(tmp_size);
	c_ptr = impl_ode_casadi_mem;
	for (int ii=0; ii<NN; ii++)
	{
		external_function_casadi_assign(impl_ode_fun+ii, c_ptr);
		c_ptr += external_function_casadi_calculate_size(impl_ode_fun+ii);
	}
	//
	tmp_size = 0;
	for (int ii=0; ii<NN; ii++)
	{
		tmp_size += external_function_casadi_calculate_size(impl_ode_fun_jac_x_xdot+ii);
	}
	void *impl_ode_fun_jac_x_xdot_mem = malloc(tmp_size);
	c_ptr = impl_ode_fun_jac_x_xdot_mem;
	for (int ii=0; ii<NN; ii++)
	{
		external_function_casadi_assign(impl_ode_fun_jac_x_xdot+ii, c_ptr);
		c_ptr += external_function_casadi_calculate_size(impl_ode_fun_jac_x_xdot+ii);
	}
	//
	tmp_size = 0;
	for (int ii=0; ii<NN; ii++)
	{
		tmp_size += external_function_casadi_calculate_size(impl_ode_fun_jac_x_xdot_u+ii);
	}
	void *impl_ode_fun_jac_x_xdot_u_mem = malloc(tmp_size);
	c_ptr = impl_ode_fun_jac_x_xdot_u_mem;
	for (int ii=0; ii<NN; ii++)
	{
		external_function_casadi_assign(impl_ode_fun_jac_x_xdot_u+ii, c_ptr);
		c_ptr += external_function_casadi_calculate_size(impl_ode_fun_jac_x_xdot_u+ii);
	}
	//
	tmp_size = 0;
	for (int ii=0; ii<NN; ii++)
	{
		tmp_size += external_function_casadi_calculate_size(impl_ode_jac_x_xdot_u+ii);
	}
	void *impl_ode_jac_x_xdot_u_mem = malloc(tmp_size);
	c_ptr = impl_ode_jac_x_xdot_u_mem;
	for (int ii=0; ii<NN; ii++)
	{
		external_function_casadi_assign(impl_ode_jac_x_xdot_u+ii, c_ptr);
		c_ptr += external_function_casadi_calculate_size(impl_ode_jac_x_xdot_u+ii);
	}

#elif DYNAMICS==3

	tmp_size = 0;
	for (int ii=0; ii<NN; ii++)
	{
		tmp_size += external_function_casadi_calculate_size(erk4_casadi+ii);
	}
	void *erk4_casadi_mem = malloc(tmp_size);
	c_ptr = erk4_casadi_mem;
	for (int ii=0; ii<NN; ii++)
	{
		external_function_casadi_assign(erk4_casadi+ii, c_ptr);
		c_ptr += external_function_casadi_calculate_size(erk4_casadi+ii);
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


	external_function_generic ls_cost_hess_generic;

	switch(NMF)
	{
		case 1:
			ls_cost_hess_generic.evaluate = &ls_cost_hess_nm2;
			break;
		case 2:
			ls_cost_hess_generic.evaluate = &ls_cost_hess_nm3;
			break;
		case 3:
			ls_cost_hess_generic.evaluate = &ls_cost_hess_nm4;
			break;
		case 4:
			ls_cost_hess_generic.evaluate = &ls_cost_hess_nm5;
			break;
		case 5:
			ls_cost_hess_generic.evaluate = &ls_cost_hess_nm6;
			break;
		default:
			printf("\nls cost hess not implemented for this numer of masses\n\n");
			exit(1);
	}
#endif


#if COST==2
	external_function_generic ext_cost_generic;

	switch(NMF)
	{
		case 1:
			ext_cost_generic.evaluate = &ext_cost_nm2;
			break;
		case 2:
			ext_cost_generic.evaluate = &ext_cost_nm3;
			break;
		case 3:
			ext_cost_generic.evaluate = &ext_cost_nm4;
			break;
		case 4:
			ext_cost_generic.evaluate = &ext_cost_nm5;
			break;
		case 5:
			ext_cost_generic.evaluate = &ext_cost_nm6;
			break;
		default:
			printf("\next cost not implemented for this numer of masses\n\n");
			exit(1);
	}
#endif

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
    * nlp_in (wip)
    ************************************************/

    // TODO(dimitris): clean up integrators inside
	tmp_size = ocp_nlp_in_calculate_size(config, dims);
	void *nlp_in_mem = malloc(tmp_size);
	ocp_nlp_in *nlp_in = ocp_nlp_in_assign(config, dims, nlp_in_mem);


	// sampling times
	for (int ii=0; ii<NN; ii++)
		nlp_in->Ts[ii] = TF/NN;




// ocp_nlp_dims_print(nlp_in->dims);

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
		blasfeo_dgese(nu[i]+nx[i], ny[i], 0.0, &cost_ls[i]->Cyt, 0, 0);
        for (int j = 0; j < nu[i]; j++)
            BLASFEO_DMATEL(&cost_ls[i]->Cyt, j, nx[i]+j) = 1.0;
        for (int j = 0; j < nx[i]; j++)
            BLASFEO_DMATEL(&cost_ls[i]->Cyt, nu[i]+j, j) = 1.0;
	}

	// W
	for (int i=0; i<=NN; i++)
	{
		blasfeo_dgese(ny[i], ny[i], 0.0, &cost_ls[i]->W, 0, 0);
        for (int j = 0; j < nx[i]; j++)
            BLASFEO_DMATEL(&cost_ls[i]->W, j, j) = diag_cost_x[j];
        for (int j = 0; j < nu[i]; j++)
            BLASFEO_DMATEL(&cost_ls[i]->W, nx[i]+j, nx[i]+j) = diag_cost_u[j];
	}

	// y_ref
    for (int i=0; i<=NN; i++)
	{
		blasfeo_pack_dvec(nx[i], xref, &cost_ls[i]->y_ref, 0);
		blasfeo_pack_dvec(nu[i], uref, &cost_ls[i]->y_ref, nx[i]);
    }



#elif COST==1

    /* nonlinear ls */

    ocp_nlp_cost_nls_model **cost_ls = (ocp_nlp_cost_nls_model **) nlp_in->cost;

	// nls_res_jac
	for (int i=0; i<=NN; i++)
		cost_ls[i]->nls_res_jac = (external_function_generic *) &ls_cost_jac_casadi[i];

	// nls_hess at first stage
	cost_ls[0]->nls_hess = &ls_cost_hess_generic;
#if 0
	// replace with hand-written external functions
	external_function_generic ls_cost_jac_generic[NN];
	if (NMF==3)
	{
		for (int i=0; i<NN; i++)
		{
			ls_cost_jac_generic[i].evaluate = &ls_cost_jac_nm4;
			cost_ls->nls_res_jac[i] = &ls_cost_jac_generic[i];
		}
	}
#endif

	// W
	for (int i=0; i<=NN; i++)
	{
		blasfeo_dgese(ny[i], ny[i], 0.0, &cost_ls[i]->W, 0, 0);
        for (int j = 0; j < nx[i]; j++)
            BLASFEO_DMATEL(&cost_ls[i]->W, j, j) = diag_cost_x[j];
        for (int j = 0; j < nu[i]; j++)
            BLASFEO_DMATEL(&cost_ls[i]->W, nx[i]+j, nx[i]+j) = diag_cost_u[j];
	}

	// y_ref
    for (int i=0; i<=NN; i++)
	{
		blasfeo_pack_dvec(nx[i], xref, &cost_ls[i]->y_ref, 0);
		blasfeo_pack_dvec(nu[i], uref, &cost_ls[i]->y_ref, nx[i]);
    }

#else
	/* external cost */

    ocp_nlp_cost_external_model **cost_external = (ocp_nlp_cost_external_model **) nlp_in->cost;

	// ext_cost
	for (int i=0; i<NN; i++)
		cost_external[i]->ext_cost = &ext_cost_generic;

    /* linear ls */

    ocp_nlp_cost_ls_model **cost_ls = (ocp_nlp_cost_ls_model **) nlp_in->cost;

	// last stage

	// Cyt
	blasfeo_dgese(nu[NN]+nx[NN], ny[NN], 0.0, &cost_ls[NN]->Cyt, 0, 0);
	for (int j = 0; j < nu[NN]; j++)
		BLASFEO_DMATEL(&cost_ls[NN]->Cyt, j, nx[NN]+j) = 1.0;
	for (int j = 0; j < nx[NN]; j++)
		BLASFEO_DMATEL(&cost_ls[NN]->Cyt, nu[NN]+j, j) = 1.0;

	// W
	blasfeo_dgese(ny[NN], ny[NN], 0.0, &cost_ls[NN]->W, 0, 0);
	for (int j = 0; j < nx[NN]; j++)
		BLASFEO_DMATEL(&cost_ls[NN]->W, j, j) = diag_cost_x[j];
	for (int j = 0; j < nu[NN]; j++)
		BLASFEO_DMATEL(&cost_ls[NN]->W, nx[NN]+j, nx[NN]+j) = diag_cost_u[j];

	// y_ref
	blasfeo_pack_dvec(nx[NN], xref, &cost_ls[NN]->y_ref, 0);
	blasfeo_pack_dvec(nu[NN], uref, &cost_ls[NN]->y_ref, nx[NN]);



#endif



	/* dynamics */
#if DYNAMICS==0
	for (int i=0; i<NN; i++)
	{
		ocp_nlp_dynamics_cont_model *dynamics = nlp_in->dynamics[i];
		erk_model *model = dynamics->sim_model;
		model->expl_vde_for = (external_function_generic *) &expl_vde_for[i];
		model->expl_ode_jac = (external_function_generic *) &expl_ode_jac[i];
	}
#elif DYNAMICS==1
	for (int i=0; i<NN; i++)
	{
		ocp_nlp_dynamics_cont_model *dynamics = nlp_in->dynamics[i];
		lifted_irk_model *model = dynamics->sim_model;
		model->expl_vde_for = (external_function_generic *) &expl_vde_for[i];
		model->expl_ode_jac = (external_function_generic *) &expl_ode_jac[i];
	}
#elif DYNAMICS==2
	for (int i=0; i<NN; i++)
	{
		ocp_nlp_dynamics_cont_model *dynamics = nlp_in->dynamics[i];
		irk_model *model = dynamics->sim_model;
		model->impl_ode_fun = (external_function_generic *) &impl_ode_fun[i];
		model->impl_ode_fun_jac_x_xdot = (external_function_generic *) &impl_ode_fun_jac_x_xdot[i];
		model->impl_ode_jac_x_xdot_u = (external_function_generic *) &impl_ode_jac_x_xdot_u[i];
	}
#elif DYNAMICS==3
	for (int i=0; i<NN; i++)
	{
		ocp_nlp_dynamics_disc_model *dynamics = nlp_in->dynamics[i];
		dynamics->disc_dyn_fun_jac = (external_function_generic *) &erk4_casadi[i];
	}
#elif DYNAMICS==4
	for (int i=0; i<NN; i++)
	{
		ocp_nlp_dynamics_cont_model *dynamics = nlp_in->dynamics[i];
		lifted_irk_model *model = dynamics->sim_model;
		model->impl_ode_fun = (external_function_generic *) &impl_ode_fun[i];
		model->impl_ode_fun_jac_x_xdot_u = (external_function_generic *) &impl_ode_fun_jac_x_xdot_u[i];
	}
#endif



    /* box constraints */

	ocp_nlp_constraints_bgh_model **constraints = (ocp_nlp_constraints_bgh_model **) nlp_in->constraints;

	// idxb0
    // int idxb0[nb[0]];
    // for (int i = 0; i < nb[0]; i++) idxb0[i] = i;

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

	// fist stage
#if CONSTRAINTS==0 // box constraints
	blasfeo_pack_dvec(nb[0], lb0, &constraints[0]->d, 0);
	blasfeo_pack_dvec(nb[0], ub0, &constraints[0]->d, nb[0]+ng[0]);
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
	blasfeo_pack_dvec(ng[0], lb0, &constraints[0]->d, nb[0]);
	blasfeo_pack_dvec(ng[0], ub0, &constraints[0]->d, 2*nb[0]+ng[0]);

	d_free(Cu0);
	d_free(Cx0);
#else // general+nonlinear constraints
	blasfeo_dgese(nu[0]+nx[0], ng[0], 0.0, &constraints[0]->DCt, 0, 0);
	for (int ii=0; ii<ng[0]; ii++)
		BLASFEO_DMATEL(&constraints[0]->DCt, ii, ii) = 1.0;

    ocp_nlp_constraints_bgh_model **nl_constr = (ocp_nlp_constraints_bgh_model **) nlp_in->constraints;
	nl_constr[0]->nl_constr_h_fun_jac = &nonlin_constr_generic;

	blasfeo_pack_dvec(ng[0]+nh[0], lb0, &constraints[0]->d, nb[0]);
	blasfeo_pack_dvec(ng[0]+nh[0], ub0, &constraints[0]->d, 2*nb[0]+ng[0]+nh[0]);
#endif

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

	tmp_size = ocp_nlp_sqp_opts_calculate_size(config, dims);
	void *nlp_opts_mem = malloc(tmp_size);
	ocp_nlp_sqp_opts *nlp_opts = ocp_nlp_sqp_opts_assign(config, dims, nlp_opts_mem);

	ocp_nlp_sqp_opts_initialize_default(config, dims, nlp_opts);

#if XCOND==1
	// partial condensing
	ocp_qp_partial_condensing_solver_opts *pcond_solver_opts = nlp_opts->qp_solver_opts;
	pcond_solver_opts->pcond_opts->N2 = 5; // set partial condensing horizon
#endif

    for (int i = 0; i < NN; ++i)
	{
#if DYNAMICS==0
		ocp_nlp_dynamics_cont_opts *dynamics_opts = nlp_opts->dynamics[i];
        sim_opts *sim_opts = dynamics_opts->sim_solver;
		// dynamics: ERK 4
		sim_opts->ns = 4;
		sim_opts->num_steps = 10;
#elif DYNAMICS==1
		ocp_nlp_dynamics_cont_opts *dynamics_opts = nlp_opts->dynamics[i];
        sim_opts *sim_opts = dynamics_opts->sim_solver;
		// dynamics: lifted IRK GL2
		sim_opts->ns = 3;
#elif DYNAMICS==2
		ocp_nlp_dynamics_cont_opts *dynamics_opts = nlp_opts->dynamics[i];
        sim_opts *sim_opts = dynamics_opts->sim_solver;
		// dynamics: IRK GL2
		sim_opts->ns = 3;
		sim_opts->jac_reuse = true;
#elif DYNAMICS==3
		// dynamics: discrete model
		// no options
#elif DYNAMICS==4
		ocp_nlp_dynamics_cont_opts *dynamics_opts = nlp_opts->dynamics[i];
        sim_opts *sim_opts = dynamics_opts->sim_solver;
		// dynamics: lifterd IRK GL2
		sim_opts->ns = 4;
#endif
    }



#if COST==1
	// exact hessian of cost function at first stage
	ocp_nlp_cost_nls_opts *cost_opts = nlp_opts->cost[0];
	cost_opts->gauss_newton_hess = 1;
#endif



	// XXX hack: overwrite config with hand-setted one
// nlp_opts->qp_solver = &config_qp;
// nlp_opts->sim_solvers = config_sim_ptrs;
// for (int ii=0; ii<NN; ii++)
// 	nlp_opts->sim_solvers[ii] = config_sim_ptrs[ii];



    nlp_opts->max_iter = MAX_SQP_ITERS;
    nlp_opts->tol_stat = 1e-9;
    nlp_opts->tol_eq = 1e-9;
    nlp_opts->tol_ineq = 1e-9;
    nlp_opts->tol_comp = 1e-9;

	// update after user-defined options
	ocp_nlp_sqp_opts_update(config, dims, nlp_opts);

    /************************************************
    * ocp_nlp out
    ************************************************/

	tmp_size = ocp_nlp_out_calculate_size(config, dims);
	void *nlp_out_mem = malloc(tmp_size);
	ocp_nlp_out *nlp_out = ocp_nlp_out_assign(config, dims, nlp_out_mem);

// ocp_nlp_dims_print(nlp_out->dims);

    /************************************************
    * sqp memory
    ************************************************/

	tmp_size = ocp_nlp_sqp_memory_calculate_size(config, dims, nlp_opts);
	void *nlp_mem_mem = malloc(tmp_size);
	ocp_nlp_sqp_memory *nlp_mem = ocp_nlp_sqp_memory_assign(config, dims, nlp_opts, nlp_mem_mem);


    /************************************************
    * sqp workspace
    ************************************************/

    int workspace_size = ocp_nlp_sqp_workspace_calculate_size(config, dims, nlp_opts);
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
// 	if (rep==0)
// 	{
			for (int i=0; i<=NN; i++)
			{
				blasfeo_pack_dvec(nu[i], uref, nlp_out->ux+i, 0);
				blasfeo_pack_dvec(nx[i], xref, nlp_out->ux+i, nu[i]);
			}
// 	}

		// call nlp solver
        status = ocp_nlp_sqp(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);
    }

    double time = acados_toc(&timer)/NREP;

	// printf("\nresiduals\n");
	// ocp_nlp_res_print(dims, nlp_mem->nlp_res);

	printf("\nsolution\n");
	ocp_nlp_out_print(dims, nlp_out);

    int sqp_iter;
    ocp_nlp_sqp_get(config, dims, nlp_mem, "sqp_iter", &sqp_iter);

    printf("\n\nstatus = %i, iterations (max %d) = %d, total time = %f ms\n\n", status, MAX_SQP_ITERS, time*1e3);

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

#elif DYNAMICS==2 | DYNAMICS==4

	free(impl_ode_casadi_mem);
	free(impl_ode_fun_jac_x_xdot_mem);
	free(impl_ode_fun_jac_x_xdot_u_mem);
	free(impl_ode_jac_x_xdot_u_mem);

#elif DYNAMICS==3

	free(erk4_casadi_mem);

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

	if (status == 0)
		printf("\nsuccess!\n\n");
	else
		printf("\nfailure!\n\n");

	return 0;

}
