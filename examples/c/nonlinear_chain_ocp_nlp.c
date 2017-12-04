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
#include "acados/utils/create.h"
#include "acados/utils/mem.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"
#include "examples/c/chain_model/chain_model.h"

#define NN 15
#define TF 3.0
#define Ns 2
#define MAX_SQP_ITERS 20
#define NREP 1

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
    for (int ii = 0; ii < nlp->dims->N; ii++)
    {
        switch (num_free_masses)
        {
            case 1:
                nlp->vde[ii] = &vde_chain_nm2;
                nlp->jac[ii] = &jac_chain_nm2;
                nlp->vde_adj[ii] = &vde_hess_chain_nm2;
                break;
            case 2:
                nlp->vde[ii] = &vde_chain_nm3;
                nlp->jac[ii] = &jac_chain_nm3;
                nlp->vde_adj[ii] = &vde_hess_chain_nm3;
                break;
            case 3:
                nlp->vde[ii] = &vde_chain_nm4;
                nlp->jac[ii] = &jac_chain_nm4;
                nlp->vde_adj[ii] = &vde_hess_chain_nm4;
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
    int nc[NN + 1] = {0};
    int ng[NN + 1] = {0};
    int ns[NN+1] = {0};
    int nh[NN+1] = {0};

    nx[0] = NX;
    nu[0] = NU;
    nbx[0] = nx[0];
    nbu[0] = nu[0];

    for (int i = 1; i < NN; i++)
    {
        nx[i] = NX;
        nu[i] = NU;
        nbx[i] = NMF;
        nbu[i] = NU;
    }

    nx[NN] = NX;
    nbx[NN] = NX;
    nbu[NN] = 0;

    // TODO(dimitris): if dims were defined stage-wise, we could do directly ocp_nlp_dims dims[N]..
    ocp_nlp_dims dims;
    dims.N  = NN;
    dims.nx = nx;
    dims.nu = nu;
    dims.ng = nc;
    dims.nbx = nbx;
    dims.nbu = nbu;
    dims.nh = nh;
    dims.ns = ns;

    /************************************************
    * nlp_in (wip)
    ************************************************/

    // TODO(dimitris): clean up integrators inside
    ocp_nlp_in *nlp = create_ocp_nlp_in(&dims, d);

    // NOTE(dimitris): use nlp->dims instead of &dims from now on since nb is filled with nbx+nbu!

    // Problem data
    double wall_pos = -0.01;
    double UMAX = 10;
    double lb0[NX+NU], ub0[NX+NU];
    read_initial_state(NX, NMF, lb0);
    read_initial_state(NX, NMF, ub0);
    for (int i = NX; i < NX+NU; i++) {
        lb0[i] = -UMAX;
        ub0[i] = +UMAX;
    }
    double lb[NMF+NU], ub[NMF+NU];
    double xref[NX];
    read_final_state(NX, NMF, xref);
    double uref[3] = {0.0, 0.0, 0.0};
    double diag_cost_x[NX];
    for (int i = 0; i < NX; i++)
        diag_cost_x[i] = 1e-2;
    double diag_cost_u[3] = {1.0, 1.0, 1.0};
    double x_pos_inf[NX], x_neg_inf[NX];
    for (int i = 0; i < NX; i++) {
        x_pos_inf[i] = +1e4;
        x_neg_inf[i] = -1e4;
    }

    // Least-squares cost
    ocp_nlp_ls_cost *ls_cost = (ocp_nlp_ls_cost *) nlp->cost;
    for (int i = 0; i < NN; i++) {
        for (int j = 0; j < NX; j++)
            ls_cost->W[i][j * (NX + NU + 1)] = diag_cost_x[j];
        for (int j = 0; j < NU; j++)
            ls_cost->W[i][(NX + j) * (NX + NU + 1)] = diag_cost_u[j];
        for (int j = 0; j < NX; j++)
            ls_cost->y_ref[i][j] = xref[j];
        for (int j = 0; j < NU; j++)
            ls_cost->y_ref[i][NX+j] = uref[j];
    }
    for (int j = 0; j < NX; j++)
        ((ocp_nlp_ls_cost *) nlp->cost)->W[NN][j * (NX + 1)] = diag_cost_x[j];

    for (int jj = 0; jj < NN; jj++)
    {
        select_model(NMF, nlp);
    }

    nlp->freezeSens = false;
    if (scheme > 2)
        nlp->freezeSens = true;

    // Box constraints
    int *nb = nlp->dims->nb;

    int idxb_0[nb[0]], idxb_1[nb[1]], idxb_N[nb[NN]];
    for (int i = 0; i < nb[0]; i++)
        idxb_0[i] = i;
    for (int i = 0; i < NMF; i++)
        idxb_1[i] = 6*i + 1;
    for (int i = 0; i < NU; i++)
        idxb_1[NMF+i] = NX+i;
    for (int i = 0; i < nb[NN]; i++)
        idxb_N[i] = i;
    nlp->lb[0] = lb0;
    nlp->ub[0] = ub0;
    nlp->idxb[0] = idxb_0;
    for (int j = 0; j < NMF; j++) {
        lb[j] = wall_pos;  // wall position
        ub[j] = 1e4;
    }
    for (int j = 0; j < NU; j++) {
        lb[NMF+j] = -UMAX;  // umin
        ub[NMF+j] = +UMAX;  // umax
    }
    for (int i = 1; i < NN; i++) {
        nlp->lb[i] = lb;
        nlp->ub[i] = ub;
        nlp->idxb[i] = idxb_1;
    }
    nlp->lb[NN] = x_neg_inf;
    nlp->ub[NN] = x_pos_inf;
    nlp->idxb[NN] = idxb_N;

    /************************************************
    * gn_sqp args
    ************************************************/

    // choose QP solver
    qp_solver_t qp_solver_name = HPIPM;

    // set up args with nested structs
    sim_solver_t sim_solver_names[NN];
    int num_stages[NN];

    for (int ii = 0; ii < NN; ii++)
    {
        sim_solver_names[ii] = LIFTED_IRK;
        num_stages[ii] = 4;
    }

    nlp->dims->num_stages = num_stages;

    ocp_nlp_gn_sqp_args *nlp_args = ocp_nlp_gn_sqp_create_args(nlp->dims, qp_solver_name, sim_solver_names);
    for (int i = 0; i < NN; ++i) {
        sim_rk_opts *sim_opts = nlp_args->sim_solvers_args[i];
        sim_opts->interval = TF/NN;
    }

    nlp_args->maxIter = MAX_SQP_ITERS;

    /************************************************
    * ocp_nlp out
    ************************************************/

    ocp_nlp_out *nlp_out = create_ocp_nlp_out(nlp->dims);

    /************************************************
    * gn_sqp memory
    ************************************************/

    ocp_nlp_gn_sqp_memory *nlp_mem = ocp_nlp_gn_sqp_create_memory(nlp->dims, nlp_args);

    // TODO(dimitris): users shouldn't write directly on memory..
    for (int i = 0; i < NN; i++) {
        for (int j = 0; j < NX; j++)
            nlp_mem->x[i][j] = xref[j];  // resX(j,i)
        for (int j = 0; j < NU; j++)
            nlp_mem->u[i][j] = uref[j];  // resU(j, i)
    }
    for (int j = 0; j < NX; j++)
        nlp_mem->x[NN][j] = xref[j];  // resX(j, NN)

    /************************************************
    * gn_sqp workspace
    ************************************************/

    int workspace_size = ocp_nlp_gn_sqp_calculate_workspace_size(nlp->dims, nlp_args);
    void *nlp_work = acados_malloc(workspace_size, 1);

    /************************************************
    * gn_sqp solve
    ************************************************/

    int status;

    acados_timer timer;
    acados_tic(&timer);

    for (int rep = 0; rep < NREP; rep++)
    {
        status = ocp_nlp_gn_sqp(nlp, nlp_out, nlp_args, nlp_mem, nlp_work);
    }

    double time = acados_toc(&timer)/NREP;

    printf("\n\nstatus = %i, iterations (fixed) = %d total time = %f ms\n\n", status, MAX_SQP_ITERS, time*1e3);

    for (int k =0; k < 3; k++) {
        printf("u[%d] = \n", k);
        d_print_mat(1, nu[k], nlp_out->u[k], 1);
        printf("x[%d] = \n", k);
        d_print_mat(1, nx[k], nlp_out->x[k], 1);
    }
    printf("u[N-1] = \n");
    d_print_mat(1, nu[NN-1], nlp_out->u[NN-1], 1);
    printf("x[N] = \n");
    d_print_mat(1, nx[NN], nlp_out->x[NN], 1);

    /************************************************
    * free memory
    ************************************************/

    free(nlp);
    free(nlp_out);
    free(nlp_work);
    free(nlp_mem);
    free(nlp_args);
}
