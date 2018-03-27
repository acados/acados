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

// #define WT_MODEL

#ifdef WT_MODEL

#include "examples/c/wt_model_nx6_expl/wt_model.h"
#include "examples/c/wt_model_nx6_expl/u_x0.c"

#else

#include "examples/c/chain_model/chain_model.h"

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

#endif

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


#ifdef WT_MODEL

static void select_dynamics_wt_casadi(int N, external_function_casadi *forw_vde)
{
	for (int ii = 0; ii < N; ii++)
	{
		forw_vde[ii].casadi_fun = &expl_forw_vde;
		forw_vde[ii].casadi_work = &expl_forw_vde_work;
		forw_vde[ii].casadi_sparsity_in = &expl_forw_vde_sparsity_in;
		forw_vde[ii].casadi_sparsity_out = &expl_forw_vde_sparsity_out;
		forw_vde[ii].casadi_n_in = &expl_forw_vde_n_in;
		forw_vde[ii].casadi_n_out = &expl_forw_vde_n_out;
	}
}



static void read_initial_state_wt(const int nx, double *x0)
{
	double *ptr = x_ref;

    for (int i = 0; i < nx; i++)
		x0[i] = ptr[i];
}



static void read_final_state_wt(const int nx, double *xN)
{
	// TODO
    for (int i = 0; i < nx; i++)
		xN[i] = 0.0;
}

#else

static void select_dynamics_chain_casadi(int N, int num_free_masses, external_function_casadi *forw_vde)
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
			}
			break;
		default:
			printf("Problem size not available\n");
			exit(1);
			break;
	}
	return;
}

static void read_initial_state_chain(const int nx, const int num_free_masses, double *x0)
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



void read_final_state_chain(const int nx, const int num_free_masses, double *xN)
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


#endif


/************************************************
* main
************************************************/

// TODO(dimitris): compile on windows

int main()
{
    // _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);

#ifdef WT_MODEL
	// TODO(eco4wind): ACADO formulation has 8 states with control of input rates
	int NX = 6;
#else
	const int NMF = 3;  // number of free masses
    int NX = 6 * NMF;
#endif

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

	// TODO(eco4wind): setup number of bounds on states and control
    nx[0] = NX;
    nu[0] = NU;
    nbx[0] = nx[0];
    nbu[0] = nu[0];
    nb[0] = nbu[0]+nbx[0];
	ng[0] = 0;

	// TODO(eco4wind): add bilinear constraints in nh
	nh[0] = 0;
#ifdef WT_MODEL
	ny[0] = 2;
#else
	ny[0] = nx[0]+nu[0];
#endif

    for (int i = 1; i < NN; i++)
    {
        nx[i] = NX;
        nu[i] = NU;
#ifdef WT_MODEL
        nbx[i] = NX;
#else
        nbx[i] = NMF;
#endif
        nbu[i] = NU;
		nb[i] = nbu[i]+nbx[i];
		ng[i] = 0;
		nh[i] = 0;
#ifdef WT_MODEL
		ny[i] = 2;
#else
		ny[i] = nx[i]+nu[i];
#endif
    }

    nx[NN] = NX;
    nu[NN] = 0;
    nbx[NN] = NX;
    nbu[NN] = 0;
    nb[NN] = nbu[NN]+nbx[NN];
	ng[NN] = 0;
	nh[NN] = 0;
#ifdef WT_MODEL
	ny[NN] = 2;
#else
	ny[NN] = nx[NN]+nu[NN];
#endif
    /************************************************
    * problem data
    ************************************************/

    double *x_end = malloc(sizeof(double)*NX);
    double *u_end = malloc(sizeof(double)*NU);

	// value of last stage when shifting states and controls
	for (int i = 0; i < NX; i++) x_end[i] = 0;
	for (int i = 0; i < NU; i++) u_end[i] = 0;

#ifdef WT_MODEL
	// TODO(eco4wind): setup bounds on controls
    double UMAX = 5;
#else
    double wall_pos = -0.01;
    double UMAX = 10;
#endif

	double x_pos_inf = +1e4;
	double x_neg_inf = -1e4;

    double *xref = malloc(NX*sizeof(double));

#ifdef WT_MODEL
    read_final_state_wt(NX, xref);
    double uref[3] = {0.0, 0.0, 0.0};
	// NOTE(dimitris): following values do not give nans at first iteration, so the problem must be
	// the correct xref/uref values that are missing
	// read_initial_state_wt(NX, xref);
	// double uref[3] = {8.169651470932033e+00, 4.024634365037572e+00, 1.399449920654297e+01};
#else
    read_final_state_chain(NX, NMF, xref);

    double uref[3] = {0.0, 0.0, 0.0};

    double *diag_cost_x = malloc(NX*sizeof(double));

    for (int i = 0; i < NX; i++)
        diag_cost_x[i] = 1e-2;

    double diag_cost_u[3] = {1.0, 1.0, 1.0};

#endif

	// idxb0
    int idxb0[nb[0]];
    for (int i = 0; i < nb[0]; i++) idxb0[i] = i;

	// idxb1
	int idxb1[nb[1]];
    for (int i = 0; i < NU; i++) idxb1[i] = i;
#ifdef WT_MODEL
    for (int i = 0; i < NX; i++) idxb1[NU+i] = NU + i;
#else
    for (int i = 0; i < NMF; i++) idxb1[NU+i] = NU + 6*i + 1;
#endif

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
#ifdef WT_MODEL
    read_initial_state_wt(NX, lb0+NU);
    read_initial_state_wt(NX, ub0+NU);
#else
    read_initial_state_chain(NX, NMF, lb0+NU);
    read_initial_state_chain(NX, NMF, ub0+NU);
#endif

	// lb1, ub1
#ifdef WT_MODEL
    double lb1[NX+NU], ub1[NX+NU];
#else
    double lb1[NMF+NU], ub1[NMF+NU];
#endif

    for (int j = 0; j < NU; j++)
	{
        lb1[j] = -UMAX;  // umin
        ub1[j] = +UMAX;  // umax
    }

#ifdef WT_MODEL
    for (int j = 0; j < NX; j++)
	{
        lb1[NU+j] = x_neg_inf;
        ub1[NU+j] = x_pos_inf;
    }
#else
    for (int j = 0; j < NMF; j++)
	{
        lb1[NU+j] = wall_pos;  // wall position
        ub1[NU+j] = x_pos_inf;
    }
#endif

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

	plan->nlp_solver = SQP_GN;

	for (int i = 0; i <= NN; i++)
		plan->nlp_cost[i] = LINEAR_LS;

	plan->ocp_qp_solver_plan.qp_solver = PARTIAL_CONDENSING_HPIPM;

	for (int i = 0; i < NN; i++)
	{
		plan->nlp_dynamics[i] = CONTINUOUS_MODEL;
		plan->sim_solver_plan[i].sim_solver = ERK;
	}

	ocp_nlp_solver_config *config = ocp_nlp_config_create(*plan, NN);

    /************************************************
    * ocp_nlp_dims
    ************************************************/

	ocp_nlp_dims *dims = ocp_nlp_dims_create(config);
	ocp_nlp_dims_initialize(config, nx, nu, ny, nbx, nbu, ng, nh, ns, nq, dims);

    /************************************************
    * dynamics
    ************************************************/

	external_function_casadi *forw_vde_casadi = malloc(NN*sizeof(external_function_casadi));

#ifdef WT_MODEL
	select_dynamics_wt_casadi(NN, forw_vde_casadi);
#else
	select_dynamics_chain_casadi(NN, NMF, forw_vde_casadi);
#endif
	external_function_casadi_create_array(NN, forw_vde_casadi);

    /************************************************
    * nlp_in
    ************************************************/

	ocp_nlp_in *nlp_in = ocp_nlp_in_create(config, dims);

	// sampling times
	for (int ii=0; ii<NN; ii++)
	{
#ifdef WT_MODEL
    	nlp_in->Ts[ii] = 0.2;
#else
		nlp_in->Ts[ii] = TF/NN;
#endif
	}

	// output definition: y = [x; u]

	/* cost */
	ocp_nlp_cost_ls_model *stage_cost_ls;

	for (int i = 0; i <= NN; i++)
	{
		switch (plan->nlp_cost[i])
		{
			case LINEAR_LS:

				stage_cost_ls = (ocp_nlp_cost_ls_model *) nlp_in->cost[i];

#ifdef WT_MODEL
				// TODO(eco4wind): set Cyt, W, y_ref for WT_MODEL

				// Cyt
				blasfeo_dgese(nu[i]+nx[i], ny[i], 0.0, &stage_cost_ls->Cyt, 0, 0);
				// penalize only 1st and 5th state
				BLASFEO_DMATEL(&stage_cost_ls->Cyt, 0, 0) = 1.0;
				BLASFEO_DMATEL(&stage_cost_ls->Cyt, 4, 1) = 1.0;

				// W
				blasfeo_dgese(ny[i], ny[i], 0.0, &stage_cost_ls->W, 0, 0);

				BLASFEO_DMATEL(&stage_cost_ls->W, 0, 0) = 1.5114;
				BLASFEO_DMATEL(&stage_cost_ls->W, 0, 1) = -0.0649;
				BLASFEO_DMATEL(&stage_cost_ls->W, 1, 0) = -0.0649;
				BLASFEO_DMATEL(&stage_cost_ls->W, 1, 1) = 0.0180;

#else
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
#endif
				break;
			default:
				break;
		}
	}

	/* dynamics */
	int set_fun_status;

	for (int i=0; i<NN; i++)
	{
		set_fun_status = nlp_set_model_in_stage(config, nlp_in, i, "forward_vde", &forw_vde_casadi[i]);
		if (set_fun_status != 0) exit(1);
		// set_fun_status = nlp_set_model_in_stage(config, nlp_in, i, "explicit_jacobian", &jac_ode_casadi[i]);
		// if (set_fun_status != 0) exit(1);
	}

    nlp_in->freezeSens = false;

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

	// TODO(eco4wind): setup bilinear constraint

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
	config->opts_update(config, dims, nlp_opts);

    /************************************************
    * ocp_nlp out
    ************************************************/

	ocp_nlp_out *nlp_out = ocp_nlp_out_create(config, dims);

	ocp_nlp_solver *solver = ocp_nlp_create(config, dims, nlp_opts);

    /************************************************
    * sqp solve
    ************************************************/

	int nmpc_problems = 3;

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

#ifdef WT_MODEL
		for (int i = 0; i <= NN; i++)
		{
			stage_cost_ls = (ocp_nlp_cost_ls_model *) nlp_in->cost[i];

			// set up y_ref
			// TODO(eco4wind): load from c file
			BLASFEO_DVECEL(&stage_cost_ls->y_ref, 0) = xref[0];
			BLASFEO_DVECEL(&stage_cost_ls->y_ref, 1) = xref[4];
		}
#endif

   	 	for (int idx = 0; idx < nmpc_problems; idx++)
		{
        	status = ocp_nlp_solve(solver, nlp_in, nlp_out);

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
	free(forw_vde_casadi);

	free(nlp_opts);
	free(nlp_in);
	free(nlp_out);
	free(solver);
	free(dims);
	free(config);
	free(plan);

	free(xref);
#ifndef WT_MODEL
	free(diag_cost_x);
#endif

	free(x_end);
	free(u_end);

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
