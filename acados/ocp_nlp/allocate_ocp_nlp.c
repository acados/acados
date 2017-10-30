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

#include "acados/ocp_nlp/allocate_ocp_nlp.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "acados/ocp_nlp/ocp_nlp_sm_common.h"
#include "acados/ocp_nlp/ocp_nlp_sm_gn.h"

#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"

static void allocate_ocp_nlp_in_basic(int_t N, int_t *nx, int_t *nu, ocp_nlp_in *const nlp) {
    int_zeros((int_t **)&nlp->nx, N + 1, 1);
    int_zeros((int_t **)&nlp->nu, N + 1, 1);
    int_zeros((int_t **)&nlp->nb, N + 1, 1);
    int_zeros((int_t **)&nlp->ng, N + 1, 1);

    nlp->N = N;
    memcpy((void *)nlp->nx, (void *)nx, sizeof(*nx) * (N + 1));
    memcpy((void *)nlp->nu, (void *)nu, sizeof(*nu) * (N + 1));
}

static void free_ocp_nlp_in_basic(ocp_nlp_in *const nlp) {
    int_free((int_t *)nlp->nx);
    int_free((int_t *)nlp->nu);
    int_free((int_t *)nlp->nb);
    int_free((int_t *)nlp->ng);
}

static void allocate_ocp_nlp_in_bounds(int_t N, int_t *nb, ocp_nlp_in *const nlp) {
    nlp->N = N;
    nlp->lb = (const real_t **)calloc(N + 1, sizeof(*nlp->lb));
    nlp->ub = (const real_t **)calloc(N + 1, sizeof(*nlp->ub));
    nlp->idxb = (const int_t **)calloc(N + 1, sizeof(*nlp->idxb));

    memcpy((void *)nlp->nb, (void *)nb, sizeof(*nb) * (N + 1));

    for (int_t i = 0; i <= N; i++) {
        int_zeros((int_t **)&nlp->idxb[i], nb[i], 1);
        d_zeros((real_t **)&nlp->lb[i], nb[i], 1);
        d_zeros((real_t **)&nlp->ub[i], nb[i], 1);
    }
}

static void free_ocp_nlp_in_bounds(ocp_nlp_in *const nlp) {
    for (int_t i = 0; i <= nlp->N; i++) {
        d_free((real_t *)nlp->lb[i]);
        d_free((real_t *)nlp->ub[i]);
        int_free((int_t *)nlp->idxb[i]);
    }
    free((real_t **)nlp->lb);
    free((real_t **)nlp->ub);
    free((int_t **)nlp->idxb);
}

static void allocate_ocp_nlp_in_nonlinear_constraints(int_t N, int_t *ng, ocp_nlp_in *const nlp) {
    nlp->path_constraints = (void **) malloc((N+1)*sizeof(ocp_nlp_function *));
}

static void free_ocp_nlp_in_nonlinear_constraints(ocp_nlp_in *const nlp) {
    free(nlp->path_constraints);
}

static void allocate_ocp_nlp_in_sim_solver(int_t N, int_t *nx, int_t *nu, ocp_nlp_in *const nlp) {
    nlp->sim = (void **)calloc(N, sizeof(sim_solver *));
    sim_solver **simulators = (sim_solver **) nlp->sim;
    for (int_t i = 0; i < N; i++) {
        simulators[i] = (sim_solver *) malloc(sizeof(sim_solver));
        int_t nx_i = nx[i];
        int_t nu_i = nu[i];
        simulators[i]->in = (sim_in *)malloc(sizeof(sim_in));
        d_zeros(&simulators[i]->in->x, nx_i, 1);
        d_zeros(&simulators[i]->in->u, nu_i, 1);
        d_zeros(&simulators[i]->in->S_forw, nx_i, nx_i + nu_i);
        for (int_t j = 0; j < nx_i; j++)
            simulators[i]->in->S_forw[j * (nx_i + 1)] = 1.0;

        d_zeros(&simulators[i]->in->S_adj, nx_i + nu_i, 1);
        d_zeros(&simulators[i]->in->grad_K, nx_i, 1);

        int_t nx_i1 = nx[i + 1];
        simulators[i]->out = (sim_out *)malloc(sizeof(sim_out));
        d_zeros(&simulators[i]->out->xn, nx_i1, 1);
        d_zeros(&simulators[i]->out->S_forw, nx_i1, nx_i + nu_i);
        d_zeros(&simulators[i]->out->grad, nx_i + nu_i, 1);
        simulators[i]->out->info = (sim_info *)malloc(sizeof(sim_info));

        simulators[i]->mem = NULL;
    }
}

static void free_ocp_nlp_in_sim_solver(ocp_nlp_in *const nlp) {
    sim_solver **simulators = (sim_solver **) nlp->sim;
    for (int_t i = 0; i < nlp->N; i++) {
        free(simulators[i]->in->x);
        free(simulators[i]->in->u);
        free(simulators[i]->in->S_forw);
        free(simulators[i]->in->S_adj);
        free(simulators[i]->in->grad_K);
        free(simulators[i]->in);

        free(simulators[i]->out->xn);
        free(simulators[i]->out->S_forw);
        free(simulators[i]->out->info);
        free(simulators[i]->out->grad);
        free(simulators[i]->out);
        free(nlp->sim[i]);
    }
    free(nlp->sim);
}

void allocate_ocp_nlp_in(int_t N, int_t *nx, int_t *nu, int_t *nb, int_t *ng,
                         ocp_nlp_in *const nlp) {
    allocate_ocp_nlp_in_basic(N, nx, nu, nlp);
    allocate_ocp_nlp_in_bounds(N, nb, nlp);
    allocate_ocp_nlp_in_nonlinear_constraints(N, ng, nlp);
    allocate_ocp_nlp_in_sim_solver(N, nx, nu, nlp);
}

void free_ocp_nlp_in(ocp_nlp_in *const nlp) {
    free_ocp_nlp_in_basic(nlp);
    free_ocp_nlp_in_bounds(nlp);
    free_ocp_nlp_in_nonlinear_constraints(nlp);
    free_ocp_nlp_in_sim_solver(nlp);
}

void allocate_ocp_nlp_out(ocp_nlp_in *const in, ocp_nlp_out *out) {
    int_t N = in->N;
    out->x = (real_t **)calloc(N + 1, sizeof(*out->x));
    out->u = (real_t **)calloc(N + 1, sizeof(*out->u));
    out->pi = (real_t **)calloc(N + 1, sizeof(*out->pi));
    out->lam = (real_t **)calloc(N + 1, sizeof(*out->lam));
    for (int_t k = 0; k < N; k++) {
        d_zeros(&out->x[k], in->nx[k], 1);
        d_zeros(&out->u[k], in->nu[k], 1);
        d_zeros(&out->pi[k], in->nx[k+1], 1);
        d_zeros(&out->lam[k], 2 * (in->nb[k] + in->ng[k]), 1);
    }
    d_zeros(&out->x[N], in->nx[N], 1);
    d_zeros(&out->u[N], in->nu[N], 1);
    d_zeros(&out->pi[N], 0, 1);
    d_zeros(&out->lam[N], 2 * (in->nb[N] + in->ng[N]), 1);
}

void free_ocp_nlp_out(int_t N, ocp_nlp_out *out) {
    for (int_t k = 0; k < N; k++) {
        d_free(out->x[k]);
        d_free(out->u[k]);
        d_free(out->pi[k]);
        d_free(out->lam[k]);
    }
    d_free(out->x[N]);
    d_free(out->u[N]);
    d_free(out->pi[N]);
    d_free(out->lam[N]);
    free(out->x);
    free(out->u);
    free(out->pi);
    free(out->lam);
}

void allocate_ls_cost(int_t N, int_t *nx, int_t *nu, int_t *ny, ocp_nlp_ls_cost *ls_cost) {
    ls_cost->N = N;
    ls_cost->W = calloc(N + 1, sizeof(*ls_cost->W));
    ls_cost->y_ref = calloc(N + 1, sizeof(*ls_cost->y_ref));
    ls_cost->fun = calloc(N + 1, sizeof(*ls_cost->fun));
    for (int_t i = 0; i <= N; i++) {
        ls_cost->W[i] = calloc(ny[i]*ny[i], sizeof(*ls_cost->W[i]));
        ls_cost->y_ref[i] = calloc(ny[i], sizeof(*ls_cost->y_ref[i]));
        ls_cost->fun[i] = malloc(sizeof(ocp_nlp_function));
        // Initialize LS cost
        ls_cost->fun[i]->nx = nx[i];
        ls_cost->fun[i]->nu = nu[i];
        ls_cost->fun[i]->np = 0;
        ls_cost->fun[i]->ny = ny[i];
        ls_cost->fun[i]->in = malloc(sizeof(casadi_wrapper_in));
        ls_cost->fun[i]->in->compute_jac = true;
        ls_cost->fun[i]->in->compute_hess = false;
        ls_cost->fun[i]->out = malloc(sizeof(casadi_wrapper_out));
        ls_cost->fun[i]->args = casadi_wrapper_create_arguments();
    }
}

void free_ls_cost(int_t N, ocp_nlp_ls_cost *ls_cost) {
    for (int_t i = 0; i <= N; i++) {
        free(ls_cost->fun[i]);
        free(ls_cost->fun[i]->in);
        free(ls_cost->fun[i]->out);
        free(ls_cost->fun[i]->args);
        casadi_wrapper_destroy(ls_cost->fun[i]->work);
        free(ls_cost->y_ref[i]);
    }
    free(ls_cost->W);
    free(ls_cost->y_ref);
    free(ls_cost->fun);
}
