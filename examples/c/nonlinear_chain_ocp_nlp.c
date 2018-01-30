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

#include "acados/ocp_nlp/ocp_nlp_gn_sqp.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/sim/sim_casadi_wrapper.h"

#include "acados/sim/sim_common.h"
#include "acados/sim/sim_erk_integrator.h"
#include "acados/sim/sim_lifted_irk_integrator.h"
#include "acados/utils/mem.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"
#include "examples/c/chain_model/chain_model.h"

#include <acados_c/legacy_create.h>

#define NN 15
#define TF 3.0
#define Ns 2
#define MAX_SQP_ITERS 10
#define NREP 1


#define BC_AS_GC


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

static void select_model(const int num_free_masses, ocp_nlp_in *nlp)
{
	ocp_nlp_model_expl *model = (ocp_nlp_model_expl *) nlp->model;
    for (int ii = 0; ii < nlp->dims->N; ii++)
    {
        switch (num_free_masses)
        {
            case 1:
                model->vde[ii] = &vde_chain_nm2;
                model->jac[ii] = &jac_chain_nm2;
                model->vde_adj[ii] = &vde_hess_chain_nm2;
                break;
            case 2:
                model->vde[ii] = &vde_chain_nm3;
                model->jac[ii] = &jac_chain_nm3;
                model->vde_adj[ii] = &vde_hess_chain_nm3;
                break;
            case 3:
                model->vde[ii] = &vde_chain_nm4;
                model->jac[ii] = &jac_chain_nm4;
                model->vde_adj[ii] = &vde_hess_chain_nm4;
                break;
            default:
                printf("Problem size not available\n");
                exit(1);
                break;
        }
    }
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
        // case 4:
        //     initial_states_file = fopen(X0_NM5_FILE, "r");
        //     break;
        // case 5:
        //     initial_states_file = fopen(X0_NM6_FILE, "r");
        //     break;
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
        // case 4:
        //     final_state_file = fopen(XN_NM5_FILE, "r");
        //     break;
        // case 5:
        //     final_state_file = fopen(XN_NM6_FILE, "r");
        //     break;
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
    int ns[NN+1] = {0};
    int nh[NN+1] = {0};
	int nv[NN+1] = {0};
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
	nv[0] = nx[0]+nu[0];
	ny[0] = nx[0]+nu[0];

    for (int i = 1; i < NN; i++)
    {
        nx[i] = NX;
        nu[i] = NU;
        nbx[i] = NMF;
        nbu[i] = NU;
		nb[i] = nbu[i]+nbx[i];
		ng[i] = 0;
		nv[i] = nx[i]+nu[i];
		ny[i] = nx[i]+nu[i];
    }

    nx[NN] = NX;
    nu[NN] = 0;
    nbx[NN] = NX;
    nbu[NN] = 0;
    nb[NN] = nbu[NN]+nbx[NN];
	ng[NN] = 0;
	nv[NN] = nx[NN]+nu[NN];
	ny[NN] = nx[NN]+nu[NN];

    /************************************************
    * ocp_nlp_dims
    ************************************************/

	/* ocp_nlp_cost_ls_dims */

	int cost_dims_size = ocp_nlp_cost_ls_dims_calculate_size(NN);
	void *cost_dims_mem = malloc(cost_dims_size);
	ocp_nlp_cost_ls_dims *cost_dims = ocp_nlp_cost_ls_dims_assign(NN, cost_dims_mem);
	ocp_nlp_cost_ls_dims_init(nv, ny, cost_dims);

	/* ocp_nlp_dims */

	int dims_size = ocp_nlp_dims_calculate_size(NN);
	void *dims_mem = malloc(dims_size);
	ocp_nlp_dims *dims = ocp_nlp_dims_assign(NN, dims_mem);
	ocp_nlp_dims_init(nx, nu, nbx, nbu, ng, nh, ns, cost_dims, dims);

//	ocp_nlp_dims_print(dims);


    /************************************************
    * nlp_in (wip)
    ************************************************/

    // TODO(dimitris): clean up integrators inside
    ocp_nlp_in *nlp_in = create_ocp_nlp_in(dims, d);

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
    ocp_nlp_cost_ls *cost_ls = (ocp_nlp_cost_ls *) nlp_in->cost;

	// nls mask
	for (int i=0; i<=NN; i++)
		cost_ls->nls_mask[i] = 0;

	// W
	for (int i=0; i<=NN; i++)
	{
		blasfeo_dgese(ny[i], ny[i], 0.0, cost_ls->W+i, 0, 0);
        for (int j = 0; j < nu[i]; j++)
            DMATEL_LIBSTR(cost_ls->W+i, j, j) = diag_cost_u[j];
        for (int j = 0; j < nx[i]; j++)
            DMATEL_LIBSTR(cost_ls->W+i, nu[i]+j, nu[i]+j) = diag_cost_x[j];
	}

	// Cyt
	for (int i=0; i<=NN; i++)
	{
		blasfeo_dgese(nv[i], ny[i], 0.0, cost_ls->Cyt+i, 0, 0);
		int n_min = ny[i]<nv[i] ? ny[i] : nv[i];
        for (int j = 0; j < n_min; j++)
            DMATEL_LIBSTR(cost_ls->Cyt+i, j, j) = 1.0;
	}

	// y_ref
    for (int i = 0; i < NN; i++)
	{
		blasfeo_pack_dvec(nu[i], uref, cost_ls->y_ref+i, 0);
		blasfeo_pack_dvec(nx[i], xref, cost_ls->y_ref+i, nu[i]);
    }

    for (int jj = 0; jj < NN; jj++)
    {
        select_model(NMF, nlp_in);
    }



    nlp_in->freezeSens = false;
    if (scheme > 2)
        nlp_in->freezeSens = true;



    // Box constraints

	// idxb0
    int idxb0[nb[0]];
    for (int i = 0; i < nb[0]; i++)
        idxb0[i] = i;

	// idxb1
	int idxb1[nb[1]];
    for (int i = 0; i < NU; i++)
        idxb1[i] = i;
    for (int i = 0; i < NMF; i++)
        idxb1[NU+i] = NU + 6*i + 1;

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
	blasfeo_pack_dvec(nb[0], lb0, nlp_in->d+0, 0);
	blasfeo_pack_dvec(nb[0], ub0, nlp_in->d+0, nb[0]+ng[0]);
    nlp_in->idxb[0] = idxb0;
    for (int i = 1; i < NN; i++)
	{
		blasfeo_pack_dvec(nb[i], lb1, nlp_in->d+i, 0);
		blasfeo_pack_dvec(nb[i], ub1, nlp_in->d+i, nb[i]+ng[i]);
        nlp_in->idxb[i] = idxb1;
    }
	blasfeo_pack_dvec(nb[NN], lbN, nlp_in->d+NN, 0);
	blasfeo_pack_dvec(nb[NN], ubN, nlp_in->d+NN, nb[NN]+ng[NN]);
    nlp_in->idxb[NN] = idxbN;


	// General constraints
	if (ng[0]>0)
	{
		double *Cu0; d_zeros(&Cu0, ng[0], nu[0]);
		for (int ii=0; ii<nu[0]; ii++)
			Cu0[ii*(ng[0]+1)] = 1.0;

		double *Cx0; d_zeros(&Cx0, ng[0], nx[0]);
		for (int ii=0; ii<nx[0]; ii++)
			Cx0[nu[0]+ii*(ng[0]+1)] = 1.0;

		blasfeo_pack_tran_dmat(ng[0], nu[0], Cu0, ng[0], nlp_in->DCt+0, 0, 0);
		blasfeo_pack_tran_dmat(ng[0], nx[0], Cx0, ng[0], nlp_in->DCt+0, nu[0], 0);
		blasfeo_pack_dvec(ng[0], lb0, nlp_in->d+0, nb[0]);
		blasfeo_pack_dvec(ng[0], ub0, nlp_in->d+0, 2*nb[0]+ng[0]);

		d_free(Cu0);
		d_free(Cx0);
	}
#if 0
	blasfeo_print_dmat(nu[0]+nx[0], ng[0], nlp_in->DCt+0, 0, 0);
	blasfeo_print_tran_dvec(2*nb[0]+2*ng[0], nlp_in->d+0, 0);
//	exit(1);
#endif

    /************************************************
    * gn_sqp args
    ************************************************/

    // choose QP solver
    ocp_qp_solver_t qp_solver_name = PARTIAL_CONDENSING_HPIPM;
//    ocp_qp_solver_t qp_solver_name = FULL_CONDENSING_HPIPM;

    // set up args with nested structs
    sim_solver_t sim_solver_names[NN];
    int num_stages[NN];

    for (int ii = 0; ii < NN; ii++)
    {
        sim_solver_names[ii] = LIFTED_IRK;
        num_stages[ii] = 4;
    }

    nlp_in->dims->num_stages = num_stages;

    ocp_nlp_gn_sqp_args *nlp_args = ocp_nlp_gn_sqp_create_args(nlp_in->dims, qp_solver_name, sim_solver_names);
    for (int i = 0; i < NN; ++i) {
        sim_rk_opts *sim_opts = nlp_args->sim_solvers_args[i];
        sim_opts->interval = TF/NN;
    }

    nlp_args->maxIter = MAX_SQP_ITERS;
    nlp_args->min_res_g = 1e-9;
    nlp_args->min_res_b = 1e-9;
    nlp_args->min_res_d = 1e-9;
    nlp_args->min_res_m = 1e-9;

    /************************************************
    * ocp_nlp out
    ************************************************/

    ocp_nlp_out *nlp_out = create_ocp_nlp_out(nlp_in->dims);

//	ocp_nlp_dims_print(nlp_out->dims);

    /************************************************
    * gn_sqp memory
    ************************************************/

    ocp_nlp_gn_sqp_memory *nlp_mem = ocp_nlp_gn_sqp_create_memory(nlp_in->dims, nlp_args);

    /************************************************
    * gn_sqp workspace
    ************************************************/

    int workspace_size = ocp_nlp_gn_sqp_calculate_workspace_size(nlp_in->dims, nlp_args);
    void *nlp_work = acados_malloc(workspace_size, 1);

    /************************************************
    * gn_sqp solve
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
        status = ocp_nlp_gn_sqp(nlp_in, nlp_out, nlp_args, nlp_mem, nlp_work);
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

	free(cost_dims_mem);
	free(dims_mem);
    free(nlp_in);
    free(nlp_out);
    free(nlp_work);
    free(nlp_mem);
    free(nlp_args);

	return 0;

}
