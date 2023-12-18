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


// standard
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// blasfeo
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"

// acados
#include "acados_c/external_function_interface.h"
#include "acados_c/ocp_nlp_interface.h"

// TODO(dimitris): use only the strictly necessary includes here
#include "acados/utils/mem.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"

#include "acados/ocp_nlp/ocp_nlp_cost_common.h"
#include "acados/ocp_nlp/ocp_nlp_constraints_bgh.h"

// model
#include "examples/c/chain_model/chain_model.h"
#include "examples/c/implicit_chain_model/chain_model_impl.h"

// x0
#include "examples/c/chain_model/x0_nm2.c"
#include "examples/c/chain_model/x0_nm3.c"
#include "examples/c/chain_model/x0_nm4.c"
#include "examples/c/chain_model/x0_nm5.c"
// #include "examples/c/chain_model/x0_nm6.c"

// xN
#include "examples/c/chain_model/xN_nm2.c"
#include "examples/c/chain_model/xN_nm3.c"
#include "examples/c/chain_model/xN_nm4.c"
#include "examples/c/chain_model/xN_nm5.c"
// #include "examples/c/chain_model/xN_nm6.c"

#define NN 15
#define TF 3.75
#define MAX_SQP_ITERS 10
#define NREP 10

// constraints (at stage 0): 0 box, 1 general
#define CONSTRAINTS 1


// TODO(dimitris): DOES THIS EVEN WORK ATM?
//enum sensitivities_scheme {
//    EXACT_NEWTON,
//    INEXACT_NEWTON,
//    INIS,
//    FROZEN_INEXACT_NEWTON,
//    FROZEN_INIS
//};



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
		// case 5:
		// 	for (int ii = 0; ii < N; ii++)
		// 	{
		// 		forw_vde[ii].casadi_fun = &vde_chain_nm6;
		// 		forw_vde[ii].casadi_work = &vde_chain_nm6_work;
		// 		forw_vde[ii].casadi_sparsity_in = &vde_chain_nm6_sparsity_in;
		// 		forw_vde[ii].casadi_sparsity_out = &vde_chain_nm6_sparsity_out;
		// 		forw_vde[ii].casadi_n_in = &vde_chain_nm6_n_in;
		// 		forw_vde[ii].casadi_n_out = &vde_chain_nm6_n_out;

		// 		impl_ode_fun[ii].casadi_fun = &casadi_impl_ode_fun_chain_nm6;
		// 		impl_ode_fun[ii].casadi_work = &casadi_impl_ode_fun_chain_nm6_work;
		// 		impl_ode_fun[ii].casadi_sparsity_in = &casadi_impl_ode_fun_chain_nm6_sparsity_in;
		// 		impl_ode_fun[ii].casadi_sparsity_out = &casadi_impl_ode_fun_chain_nm6_sparsity_out;
		// 		impl_ode_fun[ii].casadi_n_in = &casadi_impl_ode_fun_chain_nm6_n_in;
		// 		impl_ode_fun[ii].casadi_n_out = &casadi_impl_ode_fun_chain_nm6_n_out;

		// 		impl_ode_fun_jac_x_xdot[ii].casadi_fun = &casadi_impl_ode_fun_jac_x_xdot_chain_nm6;
		// 		impl_ode_fun_jac_x_xdot[ii].casadi_work = &casadi_impl_ode_fun_jac_x_xdot_chain_nm6_work;
		// 		impl_ode_fun_jac_x_xdot[ii].casadi_sparsity_in = &casadi_impl_ode_fun_jac_x_xdot_chain_nm6_sparsity_in;
		// 		impl_ode_fun_jac_x_xdot[ii].casadi_sparsity_out = &casadi_impl_ode_fun_jac_x_xdot_chain_nm6_sparsity_out;
		// 		impl_ode_fun_jac_x_xdot[ii].casadi_n_in = &casadi_impl_ode_fun_jac_x_xdot_chain_nm6_n_in;
		// 		impl_ode_fun_jac_x_xdot[ii].casadi_n_out = &casadi_impl_ode_fun_jac_x_xdot_chain_nm6_n_out;

		// 		impl_ode_fun_jac_x_xdot_u[ii].casadi_fun = &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm6;
		// 		impl_ode_fun_jac_x_xdot_u[ii].casadi_work = &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm6_work;
		// 		impl_ode_fun_jac_x_xdot_u[ii].casadi_sparsity_in = &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm6_sparsity_in;
		// 		impl_ode_fun_jac_x_xdot_u[ii].casadi_sparsity_out = &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm6_sparsity_out;
		// 		impl_ode_fun_jac_x_xdot_u[ii].casadi_n_in = &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm6_n_in;
		// 		impl_ode_fun_jac_x_xdot_u[ii].casadi_n_out = &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm6_n_out;

		// 		impl_ode_jac_x_xdot_u[ii].casadi_fun = &casadi_impl_ode_jac_x_xdot_u_chain_nm6;
		// 		impl_ode_jac_x_xdot_u[ii].casadi_work = &casadi_impl_ode_jac_x_xdot_u_chain_nm6_work;
		// 		impl_ode_jac_x_xdot_u[ii].casadi_sparsity_in = &casadi_impl_ode_jac_x_xdot_u_chain_nm6_sparsity_in;
		// 		impl_ode_jac_x_xdot_u[ii].casadi_sparsity_out = &casadi_impl_ode_jac_x_xdot_u_chain_nm6_sparsity_out;
		// 		impl_ode_jac_x_xdot_u[ii].casadi_n_in = &casadi_impl_ode_jac_x_xdot_u_chain_nm6_n_in;
		// 		impl_ode_jac_x_xdot_u[ii].casadi_n_out = &casadi_impl_ode_jac_x_xdot_u_chain_nm6_n_out;

		// 		// erk4_casadi[ii].casadi_fun = &casadi_erk4_chain_nm6;
		// 		// erk4_casadi[ii].casadi_work = &casadi_erk4_chain_nm6_work;
		// 		// erk4_casadi[ii].casadi_sparsity_in = &casadi_erk4_chain_nm6_sparsity_in;
		// 		// erk4_casadi[ii].casadi_sparsity_out = &casadi_erk4_chain_nm6_sparsity_out;
		// 		// erk4_casadi[ii].casadi_n_in = &casadi_erk4_chain_nm6_n_in;
		// 		// erk4_casadi[ii].casadi_n_out = &casadi_erk4_chain_nm6_n_out;
		// 	}
			break;
		default:
			printf("Problem size not available\n");
			exit(1);
			break;
	}
	return;
}


static void select_external_stage_cost_casadi(int indx, int N, int num_free_masses, external_function_casadi *external_cost)
{
	switch (num_free_masses)
	{
		case 1:
			if (indx < N)
			{
				external_cost->casadi_fun = &chain_nm_2_external_cost;
				external_cost->casadi_work = &chain_nm_2_external_cost_work;
				external_cost->casadi_sparsity_in = &chain_nm_2_external_cost_sparsity_in;
				external_cost->casadi_sparsity_out = &chain_nm_2_external_cost_sparsity_out;
				external_cost->casadi_n_in = &chain_nm_2_external_cost_n_in;
				external_cost->casadi_n_out = &chain_nm_2_external_cost_n_out;
			}
			else
			{
				printf("external cost not implemented for final stage");
				exit(1);
			}
			break;
		case 2:
			if (indx < N)
			{
				external_cost->casadi_fun = &chain_nm_3_external_cost;
				external_cost->casadi_work = &chain_nm_3_external_cost_work;
				external_cost->casadi_sparsity_in = &chain_nm_3_external_cost_sparsity_in;
				external_cost->casadi_sparsity_out = &chain_nm_3_external_cost_sparsity_out;
				external_cost->casadi_n_in = &chain_nm_3_external_cost_n_in;
				external_cost->casadi_n_out = &chain_nm_3_external_cost_n_out;
			}
			else
			{
				printf("external cost not implemented for final stage");
				exit(1);
			}
			break;
		case 3:
			if (indx < N)
			{
				external_cost->casadi_fun = &chain_nm_4_external_cost;
				external_cost->casadi_work = &chain_nm_4_external_cost_work;
				external_cost->casadi_sparsity_in = &chain_nm_4_external_cost_sparsity_in;
				external_cost->casadi_sparsity_out = &chain_nm_4_external_cost_sparsity_out;
				external_cost->casadi_n_in = &chain_nm_4_external_cost_n_in;
				external_cost->casadi_n_out = &chain_nm_4_external_cost_n_out;
			}
			else
			{
				printf("external cost not implemented for final stage");
				exit(1);
			}
			break;
		case 4:
			if (indx < N)
			{
				external_cost->casadi_fun = &chain_nm_5_external_cost;
				external_cost->casadi_work = &chain_nm_5_external_cost_work;
				external_cost->casadi_sparsity_in = &chain_nm_5_external_cost_sparsity_in;
				external_cost->casadi_sparsity_out = &chain_nm_5_external_cost_sparsity_out;
				external_cost->casadi_n_in = &chain_nm_5_external_cost_n_in;
				external_cost->casadi_n_out = &chain_nm_5_external_cost_n_out;
			}
			else
			{
				printf("external cost not implemented for final stage");
				exit(1);
			}
			break;
		case 5:
			if (indx < N)
			{
				external_cost->casadi_fun = &chain_nm_6_external_cost;
				external_cost->casadi_work = &chain_nm_6_external_cost_work;
				external_cost->casadi_sparsity_in = &chain_nm_6_external_cost_sparsity_in;
				external_cost->casadi_sparsity_out = &chain_nm_6_external_cost_sparsity_out;
				external_cost->casadi_n_in = &chain_nm_6_external_cost_n_in;
				external_cost->casadi_n_out = &chain_nm_6_external_cost_n_out;
			}
			else
			{
				printf("external cost not implemented for final stage");
				exit(1);
			}
			break;
	}
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
		// case 5:
		// 	if (indx < N)
		// 	{
		// 		ls_cost_jac->casadi_fun = &ls_cost_nm6;
		// 		ls_cost_jac->casadi_work = &ls_cost_nm6_work;
		// 		ls_cost_jac->casadi_sparsity_in = &ls_cost_nm6_sparsity_in;
		// 		ls_cost_jac->casadi_sparsity_out = &ls_cost_nm6_sparsity_out;
		// 		ls_cost_jac->casadi_n_in = &ls_cost_nm6_n_in;
		// 		ls_cost_jac->casadi_n_out = &ls_cost_nm6_n_out;
		// 	}
		// 	else
		// 	{
		// 		ls_cost_jac->casadi_fun = &ls_costN_nm6;
		// 		ls_cost_jac->casadi_work = &ls_costN_nm6_work;
		// 		ls_cost_jac->casadi_sparsity_in = &ls_costN_nm6_sparsity_in;
		// 		ls_cost_jac->casadi_sparsity_out = &ls_costN_nm6_sparsity_out;
		// 		ls_cost_jac->casadi_n_in = &ls_costN_nm6_n_in;
		// 		ls_cost_jac->casadi_n_out = &ls_costN_nm6_n_out;
		// 	}
		// 	break;
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
        // case 5:
        //     ptr = x0_nm6;
        //     break;
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
        // case 5:
        //     ptr = xN_nm6;
        //     break;
        default:
            printf("\nwrong number of free masses\n");
			exit(1);
            break;
    }
    for (int i = 0; i < nx; i++)
		xN[i] = ptr[i];
}





// // hand-wirtten box constraints on states as nonlinear constraints
// void nonlin_constr_nm2(void *evaluate, ext_fun_arg_t *type_in, void **in, ext_fun_arg_t *type_out, void **out)
// {

// 	int ii;

// 	int nu = 3;
// 	int nx = 6;

// 	int nh = nx;

// 	// fun
// 	struct blasfeo_dvec_args *fun_args = out[0];
// 	struct blasfeo_dvec *fun = fun_args->x;
// 	int xi = fun_args->xi;
// 	struct blasfeo_dvec *ux = in[0];
// 	blasfeo_dveccp(nx, ux, nu, fun, xi);

// 	// jacobian
// 	struct blasfeo_dmat_args *jac_args = out[1];
// 	struct blasfeo_dmat *jac = jac_args->A;
// 	int ai = jac_args->ai;
// 	int aj = jac_args->aj;
// 	blasfeo_dgese(nu+nx, nh, 0.0, jac, ai, aj);
// 	for (ii=0; ii<nh; ii++)
// 		BLASFEO_DMATEL(jac, ai+nu+ii, aj+ii) = 1.0;

// 	return;

// }

// void nonlin_constr_nm3(void *evaluate, ext_fun_arg_t *type_in, void **in, ext_fun_arg_t *type_out, void **out)
// {

// 	int ii;

// 	int nu = 3;
// 	int nx = 12;

// 	int nh = nx;

// 	// fun
// 	struct blasfeo_dvec_args *fun_args = out[0];
// 	struct blasfeo_dvec *fun = fun_args->x;
// 	int xi = fun_args->xi;
// 	struct blasfeo_dvec *ux = in[0];
// 	blasfeo_dveccp(nx, ux, nu, fun, xi);

// 	// jacobian
// 	struct blasfeo_dmat_args *jac_args = out[1];
// 	struct blasfeo_dmat *jac = jac_args->A;
// 	int ai = jac_args->ai;
// 	int aj = jac_args->aj;
// 	blasfeo_dgese(nu+nx, nh, 0.0, jac, ai, aj);
// 	for (ii=0; ii<nh; ii++)
// 		BLASFEO_DMATEL(jac, ai+nu+ii, aj+ii) = 1.0;

// 	return;

// }



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

	double y_ref[NX+3];
	for (int ii = 0; ii < NX; ii++) {
		y_ref[ii] = xref[ii];
	}
	for (int ii = 0; ii < 3; ii++) {
		y_ref[ii+NX] = uref[ii];
	}

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

	// Cyt
	double *Cyt = malloc((NX+NU)*(NX+NU)*sizeof(double));
	for (int j = 0; j < (NX+NU)*(NX+NU); j++)
		Cyt[j] = 0.0;
	for (int j = 0; j < NU; j++)
		Cyt[j+(NX+NU)*(j+NX)] = 1.0;
	for (int j = 0; j < NX; j++)
		Cyt[NU+j+(NX+NU)*j] = 1.0;

	// CytN
	double *CytN = malloc((NX)*(NX)*sizeof(double));
	for (int j = 0; j < NX; j++)
		CytN[j+NX*j] = 1.0;

	// W
	double *W = malloc((NX+NU)*(NX+NU)*sizeof(double));
	for (int j = 0; j < (NX+NU)*(NX+NU); j++)
		W[j] = 0.0;
	for (int j = 0; j < NX; j++)
		W[j + (NX+NU) * j] = diag_cost_x[j];
	for (int j = 0; j < NU; j++)
		W[NX+j + (NX+NU) * (NX+j)] = diag_cost_u[j];

	// WN
	double *WN = malloc(NX*NX*sizeof(double));
	for (int j = 0; j < NX*NX; j++)
		WN[j] = 0.0;
	for (int j = 0; j < NX; j++)
		WN[j + NX * j] = diag_cost_x[j];

	double *specific_u = malloc(NU*sizeof(double));
	double *specific_x = malloc(NX*sizeof(double));

    /************************************************
    * plan + config
    ************************************************/

	ocp_nlp_plan_t *plan = ocp_nlp_plan_create(NN);

	// TODO(dimitris): not necessarily GN, depends on cost module
	plan->nlp_solver = SQP;

	plan->regularization = NO_REGULARIZE;

	// NOTE(dimitris): switching between different objectives on each stage to test everything
	for (int i = 0; i <= NN; i++)
	{
		if (i < 3)
			plan->nlp_cost[i] = EXTERNAL;  // also implements linear LS for this example
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
	// plan->ocp_qp_solver_plan.qp_solver = PARTIAL_CONDENSING_OSQP;

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

	ocp_nlp_config *config = ocp_nlp_config_create(*plan);

    /************************************************
    * ocp_nlp_dims
    ************************************************/

	ocp_nlp_dims *dims = ocp_nlp_dims_create(config);

    ocp_nlp_dims_set_opt_vars(config, dims, "nx", nx);
    ocp_nlp_dims_set_opt_vars(config, dims, "nu", nu);
    ocp_nlp_dims_set_opt_vars(config, dims, "nz", nz);
    ocp_nlp_dims_set_opt_vars(config, dims, "ns", ns);

	for (int i = 0; i <= NN; i++)
    {
		if (plan->nlp_cost[i] != EXTERNAL)
		{
	        ocp_nlp_dims_set_cost(config, dims, i, "ny", &ny[i]);
		}

        ocp_nlp_dims_set_constraints(config, dims, i, "nbx", &nbx[i]);
        ocp_nlp_dims_set_constraints(config, dims, i, "nbu", &nbu[i]);
        ocp_nlp_dims_set_constraints(config, dims, i, "ng", &ng[i]);
        ocp_nlp_dims_set_constraints(config, dims, i, "nh", &nh[i]);
    }

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

	select_dynamics_casadi(NN, NMF, expl_vde_for, impl_ode_fun, impl_ode_fun_jac_x_xdot,
	                       impl_ode_fun_jac_x_xdot_u, impl_ode_jac_x_xdot_u, erk4_casadi);

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
	external_function_casadi *external_cost = malloc(NN*sizeof(external_function_casadi));

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

			case EXTERNAL:
				select_external_stage_cost_casadi(i, NN, NMF, &external_cost[i]);
				external_function_casadi_create(&external_cost[i]);
				break;

			default:
				printf("\ninvalid cost module\n\n");
				exit(1);	
		}
	}

//     /************************************************
//     * nonlinear constraints
//     ************************************************/

// #if CONSTRAINTS==2
// 	external_function_generic nonlin_constr_generic;

// 	switch(NMF)
// 	{
// 		case 1:
// 			nonlin_constr_generic.evaluate = &nonlin_constr_nm2;
// 			break;
// 		case 2:
// 			nonlin_constr_generic.evaluate = &nonlin_constr_nm3;
// 			break;
// 		case 3:
// 			nonlin_constr_generic.evaluate = &nonlin_constr_nm4;
// 			break;
// 		case 4:
// 			nonlin_constr_generic.evaluate = &nonlin_constr_nm5;
// 			break;
// 		case 5:
// 			nonlin_constr_generic.evaluate = &nonlin_constr_nm6;
// 			break;
// 		default:
// 			printf("\nnonlin constr not implemented for this numer of masses\n\n");
// 			exit(1);
// 	}
// #endif

    /************************************************
    * nlp_in
    ************************************************/

	ocp_nlp_in *nlp_in = ocp_nlp_in_create(config, dims);

	// sampling times
	for (int ii=0; ii<NN; ii++)
		nlp_in->Ts[ii] = TF/NN;

	// output definition: y = [x; u]

	/* cost */

	for (int i = 0; i <= NN; i++)
	{
		switch (plan->nlp_cost[i])
		{
			case LINEAR_LS:
				// Cyt
				if (i < NN)
		            ocp_nlp_cost_model_set(config, dims, nlp_in, i, "Cyt", Cyt);
				else
					ocp_nlp_cost_model_set(config, dims, nlp_in, i, "Cyt", CytN);

				// W
				if (i < NN)
		            ocp_nlp_cost_model_set(config, dims, nlp_in, i, "W", W);
				else
					ocp_nlp_cost_model_set(config, dims, nlp_in, i, "W", WN);

				// y_ref
				ocp_nlp_cost_model_set(config, dims, nlp_in, i, "yref", y_ref);
				break;

			case NONLINEAR_LS:
				// nls_res_jac
				ocp_nlp_cost_model_set(config, dims, nlp_in, i, "nls_res_jac", &ls_cost_jac_casadi[i]);

				// W
				if (i < NN)
		            ocp_nlp_cost_model_set(config, dims, nlp_in, i, "W", W);
				else
					ocp_nlp_cost_model_set(config, dims, nlp_in, i, "W", WN);

				// y_ref
				ocp_nlp_cost_model_set(config, dims, nlp_in, i, "yref", y_ref);
				break;

			case EXTERNAL:

				ocp_nlp_cost_model_set(config, dims, nlp_in, i, "ext_cost_fun_jac_hes", &external_cost[i]);

				assert(i < NN && "externally provided cost not implemented for last stage!");

				break;

			default:
				printf("\ninvalid cost module\n\n");
				exit(1);	
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
					set_fun_status = ocp_nlp_dynamics_model_set(config, dims, nlp_in, i, "expl_vde_for", &expl_vde_for[i]);
					if (set_fun_status != 0) exit(1);
				}
				else if (plan->sim_solver_plan[i].sim_solver == IRK)
				{
					set_fun_status = ocp_nlp_dynamics_model_set(config, dims, nlp_in, i, "impl_ode_fun", &impl_ode_fun[i]);
					if (set_fun_status != 0) exit(1);
					set_fun_status = ocp_nlp_dynamics_model_set(config, dims, nlp_in, i, "impl_ode_fun_jac_x_xdot", &impl_ode_fun_jac_x_xdot[i]);
					if (set_fun_status != 0) exit(1);
					set_fun_status = ocp_nlp_dynamics_model_set(config, dims, nlp_in, i, "impl_ode_jac_x_xdot_u", &impl_ode_jac_x_xdot_u[i]);
					if (set_fun_status != 0) exit(1);
				}
				else if (plan->sim_solver_plan[i].sim_solver == LIFTED_IRK)
				{
					set_fun_status = ocp_nlp_dynamics_model_set(config, dims, nlp_in, i, "impl_ode_fun", &impl_ode_fun[i]);
					if (set_fun_status != 0) exit(1);
					set_fun_status = ocp_nlp_dynamics_model_set(config, dims, nlp_in, i, "impl_ode_fun_jac_x_xdot_u", &impl_ode_fun_jac_x_xdot_u[i]);
					if (set_fun_status != 0) exit(1);
				}
				break;

			case DISCRETE_MODEL:
				set_fun_status = ocp_nlp_dynamics_model_set(config, dims, nlp_in, i, "disc_dyn_fun_jac", &erk4_casadi[i]);
				if (set_fun_status != 0) exit(1);
				break;

			default:
				printf("\ninvalid dynamics module\n\n");
				exit(1);	
		}
	}


    /* constraints */
	ocp_nlp_constraints_bgh_model **constraints = (ocp_nlp_constraints_bgh_model **) nlp_in->constraints;

	// fist stage
#if CONSTRAINTS==0 // box constraints
	ocp_nlp_constraints_model_set(config, dims, nlp_in, 0, "lb", lb0);
	ocp_nlp_constraints_model_set(config, dims, nlp_in, 0, "ub", ub0);
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
	ocp_nlp_constraints_model_set(config, dims, nlp_in, 0, "lg", lb0);
	ocp_nlp_constraints_model_set(config, dims, nlp_in, 0, "ug", ub0);

	d_free(Cu0);
	d_free(Cx0);
// #else // general + nonlinear constraints
// 	blasfeo_dgese(nu[0]+nx[0], ng[0], 0.0, &constraints[0]->DCt, 0, 0);
// 	for (int ii=0; ii<ng[0]; ii++)
// 		BLASFEO_DMATEL(&constraints[0]->DCt, ii, ii) = 1.0;

//     ocp_nlp_constraints_bgh_model **nl_constr = (ocp_nlp_constraints_bgh_model **) nlp_in->constraints;
// 	nl_constr[0]->nl_constr_h_fun_jac = &nonlin_constr_generic;

// 	ocp_nlp_constraints_model_set(config, dims, nlp_in, 0, "lg", lb0);
// 	ocp_nlp_constraints_model_set(config, dims, nlp_in, 0, "ug", ub0);
// 	ocp_nlp_constraints_model_set(config, dims, nlp_in, 0, "lh", &lb0[ng[0]]);
// 	ocp_nlp_constraints_model_set(config, dims, nlp_in, 0, "uh", &ub0[ng[0]]);
#endif

	// other stages
    for (int i = 1; i < NN; i++)
	{
		ocp_nlp_constraints_model_set(config, dims, nlp_in, i, "lb", lb1);
		ocp_nlp_constraints_model_set(config, dims, nlp_in, i, "ub", ub1);
        constraints[i]->idxb = idxb1;
    }
	ocp_nlp_constraints_model_set(config, dims, nlp_in, NN, "lb", lbN);
	ocp_nlp_constraints_model_set(config, dims, nlp_in, NN, "ub", ubN);

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

	void *nlp_opts = ocp_nlp_solver_opts_create(config, dims);

    for (int i = 0; i < NN; ++i)
	{
		if (plan->nlp_dynamics[i] == CONTINUOUS_MODEL)
		{
			if (plan->sim_solver_plan[i].sim_solver == ERK)
			{
				int ns = 4;

				ocp_nlp_solver_opts_set_at_stage(config, nlp_opts, i, "dynamics_ns", &ns);
			}
			else if (plan->sim_solver_plan[i].sim_solver == IRK)
			{
				int ns = 2;
				bool jac_reuse = true;

				ocp_nlp_solver_opts_set_at_stage(config, nlp_opts, i, "dynamics_ns", &ns);
				ocp_nlp_solver_opts_set_at_stage(config, nlp_opts, i, "dynamics_jac_reuse", &jac_reuse);
			}
		}
    }

    int max_iter = MAX_SQP_ITERS;
    double tol_stat = 1e-9;
    double tol_eq   = 1e-9;
    double tol_ineq = 1e-9;
    double tol_comp = 1e-9;

    ocp_nlp_solver_opts_set(config, nlp_opts, "max_iter", &max_iter);
    ocp_nlp_solver_opts_set(config, nlp_opts, "tol_stat", &tol_stat);
    ocp_nlp_solver_opts_set(config, nlp_opts, "tol_eq", &tol_eq);
    ocp_nlp_solver_opts_set(config, nlp_opts, "tol_ineq", &tol_ineq);
    ocp_nlp_solver_opts_set(config, nlp_opts, "tol_comp", &tol_comp);

    /************************************************
    * ocp_nlp out
    ************************************************/

	ocp_nlp_out *nlp_out = ocp_nlp_out_create(config, dims);

	ocp_nlp_solver *solver = ocp_nlp_solver_create(config, dims, nlp_opts);
	int status = ocp_nlp_precompute(solver, nlp_in, nlp_out);

    /************************************************
    * sqp solve
    ************************************************/

    acados_timer timer;
    acados_tic(&timer);

    for (int rep = 0; rep < NREP; rep++)
    {
		// warm start output initial guess of solution
		for (int i=0; i<=NN; i++)
		{
			blasfeo_pack_dvec(nu[i], uref, 1, nlp_out->ux+i, 0);
			blasfeo_pack_dvec(nx[i], xref, 1, nlp_out->ux+i, nu[i]);
		}

		// call nlp solver
        status = ocp_nlp_solve(solver, nlp_in, nlp_out);
    }

    double time = acados_toc(&timer)/NREP;

    ocp_nlp_res *residual;
    ocp_nlp_get(config, solver, "nlp_res", &residual);
    printf("\nresiduals\n");
    ocp_nlp_res_print(dims, residual);

	printf("\nsolution\n");
	ocp_nlp_out_print(dims, nlp_out);

	int sqp_iter;
    double time_lin, time_qp_sol, time_tot;

    ocp_nlp_get(config, solver, "sqp_iter", &sqp_iter);
    ocp_nlp_get(config, solver, "time_tot", &time_tot);
    ocp_nlp_get(config, solver, "time_qp_sol", &time_qp_sol);
    ocp_nlp_get(config, solver, "time_lin", &time_lin);

    printf("\n\nstatus = %i, iterations (max %d) = %d, total time = %f ms\n", status, MAX_SQP_ITERS, sqp_iter, time*1e3);
	printf("\nlinearization time = %f ms\n", time_lin*1e3);
	printf("\nqp solution time   = %f ms\n", time_qp_sol*1e3);
	printf("\ntotal time         = %f ms\n", time_tot*1e3);
	printf("\n\n");

    for (int k =0; k < 3; k++) {
        printf("u[%d] = \n", k);
        ocp_nlp_out_get(config, dims, nlp_out, k, "u", specific_u);
        d_print_mat(1, NU, specific_u, 1);

        printf("x[%d] = \n", k);
        ocp_nlp_out_get(config, dims, nlp_out, k, "x", specific_x);
        d_print_mat(1, NX, specific_x, 1);
    }
    printf("u[N-1] = \n");
    ocp_nlp_out_get(config, dims, nlp_out, NN-1, "u", specific_u);
    d_print_mat(1, NU, specific_u, 1);

    printf("x[N] = \n");
    ocp_nlp_out_get(config, dims, nlp_out, NN, "x", specific_x);
    d_print_mat(1, NX, specific_x, 1);

    /************************************************
    * free memory
    ************************************************/

	// free memory of the external functions
 	external_function_casadi_free(expl_vde_for);
	external_function_casadi_free_array(NN, impl_ode_fun);
	external_function_casadi_free_array(NN, impl_ode_fun_jac_x_xdot);
	external_function_casadi_free_array(NN, impl_ode_fun_jac_x_xdot_u);
	external_function_casadi_free_array(NN, impl_ode_jac_x_xdot_u);
	external_function_casadi_free_array(NN, erk4_casadi);

	// free pointers to the external functions
	free(expl_vde_for);
	free(impl_ode_fun);
	free(impl_ode_fun_jac_x_xdot);
	free(impl_ode_fun_jac_x_xdot_u);
	free(impl_ode_jac_x_xdot_u);
	free(erk4_casadi);

	// free ocp_nlp module
	ocp_nlp_solver_opts_destroy(nlp_opts);
	ocp_nlp_in_destroy(nlp_in);
	ocp_nlp_out_destroy(nlp_out);
	ocp_nlp_solver_destroy(solver);
	ocp_nlp_dims_destroy(dims);
	ocp_nlp_config_destroy(config);

	// free memory allocated in main
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

	free(CytN);
	free(Cyt);
	free(WN);
	free(W);
    free(specific_x);
    free(specific_u);

	for (int i = 0; i <= NN; i++)
	{
		switch (plan->nlp_cost[i])
		{
			case NONLINEAR_LS:
				external_function_casadi_free(&ls_cost_jac_casadi[i]);
				break;
			case EXTERNAL:
				external_function_casadi_free(&external_cost[i]);
			default:
				break;
		}
	}

	ocp_nlp_plan_destroy(plan);

	free(ls_cost_jac_casadi);
	free(external_cost);


	/************************************************
	* return
	************************************************/

	if (status == 0)
		printf("\nsuccess! (%d iter) \n\n", sqp_iter);
	else
		printf("\nfailure!\n\n");

	return 0;
}
