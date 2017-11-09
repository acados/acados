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

#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"

// static void allocate_ocp_nlp_in_basic(int_t N, int_t *nx, int_t *nu, ocp_nlp_in *const nlp) {
//     int_zeros((int_t **)&nlp->nx, N + 1, 1);
//     int_zeros((int_t **)&nlp->nu, N + 1, 1);
//     int_zeros((int_t **)&nlp->nb, N + 1, 1);
//     int_zeros((int_t **)&nlp->nc, N + 1, 1);
//     int_zeros((int_t **)&nlp->ng, N + 1, 1);

//     nlp->N = N;
//     memcpy((void *)nlp->nx, (void *)nx, sizeof(*nx) * (N + 1));
//     memcpy((void *)nlp->nu, (void *)nu, sizeof(*nu) * (N + 1));
// }

// static void free_ocp_nlp_in_basic(ocp_nlp_in *const nlp) {
//     int_free((int_t *)nlp->nx);
//     int_free((int_t *)nlp->nu);
//     int_free((int_t *)nlp->nb);
//     int_free((int_t *)nlp->nc);
//     int_free((int_t *)nlp->ng);
// }

// static void allocate_ocp_nlp_in_bounds(int_t N, int_t *nb, ocp_nlp_in *const nlp) {
//     nlp->N = N;
//     nlp->lb = (const real_t **)calloc(N + 1, sizeof(*nlp->lb));
//     nlp->ub = (const real_t **)calloc(N + 1, sizeof(*nlp->ub));
//     nlp->idxb = (const int_t **)calloc(N + 1, sizeof(*nlp->idxb));

//     memcpy((void *)nlp->nb, (void *)nb, sizeof(*nb) * (N + 1));

//     for (int_t i = 0; i <= N; i++) {
//         int_zeros((int_t **)&nlp->idxb[i], nb[i], 1);
//         d_zeros((real_t **)&nlp->lb[i], nb[i], 1);
//         d_zeros((real_t **)&nlp->ub[i], nb[i], 1);
//     }
// }

// static void free_ocp_nlp_in_bounds(ocp_nlp_in *const nlp) {
//     for (int_t i = 0; i <= nlp->N; i++) {
//         d_free((real_t *)nlp->lb[i]);
//         d_free((real_t *)nlp->ub[i]);
//         int_free((int_t *)nlp->idxb[i]);
//     }
//     free((real_t **)nlp->lb);
//     free((real_t **)nlp->ub);
//     free((int_t **)nlp->idxb);
// }

// static void allocate_ocp_nlp_in_polyhedral(int_t N, int_t *nc,
//                                            ocp_nlp_in *const nlp) {
//     nlp->lc = (const real_t **)calloc(N + 1, sizeof(*nlp->lc));
//     nlp->uc = (const real_t **)calloc(N + 1, sizeof(*nlp->uc));
//     nlp->Cx = (const real_t **)calloc(N + 1, sizeof(*nlp->Cx));
//     nlp->Cu = (const real_t **)calloc(N + 1, sizeof(*nlp->Cu));

//     nlp->N = N;
//     memcpy((void *)nlp->nc, (void *)nc, sizeof(*nc) * (N + 1));

//     for (int_t i = 0; i <= N; i++) {
//         d_zeros((real_t **)&nlp->lc[i], nc[i], 1);
//         d_zeros((real_t **)&nlp->uc[i], nc[i], 1);
//         d_zeros((real_t **)&nlp->Cx[i], nc[i], nlp->nx[i]);
//         d_zeros((real_t **)&nlp->Cu[i], nc[i], nlp->nu[i]);
//     }
// }

// static void free_ocp_nlp_in_polyhedral(ocp_nlp_in *const nlp) {
//     for (int_t i = 0; i <= nlp->N; i++) {
//         d_free((real_t *)nlp->lc[i]);
//         d_free((real_t *)nlp->uc[i]);
//         d_free((real_t *)nlp->Cx[i]);
//         d_free((real_t *)nlp->Cu[i]);
//     }

//     free((real_t **)nlp->lc);
//     free((real_t **)nlp->uc);
//     free((real_t **)nlp->Cx);
//     free((real_t **)nlp->Cu);
// }

// static void allocate_ocp_nlp_in_nonlinear_constraints(int_t N, int_t *ng, ocp_nlp_in *const nlp) {
//     nlp->lg = (const real_t **)calloc(N + 1, sizeof(*nlp->lc));
//     nlp->ug = (const real_t **)calloc(N + 1, sizeof(*nlp->uc));

//     nlp->N = N;
//     memcpy((void *)nlp->ng, (void *)ng, sizeof(*ng) * (N + 1));
//     for (int_t i = 0; i <= N; i++) {
//         d_zeros((real_t **)&nlp->lg[i], ng[i], 1);
//         d_zeros((real_t **)&nlp->ug[i], ng[i], 1);
//     }
// }

// static void free_ocp_nlp_in_nonlinear_constraints(ocp_nlp_in *const nlp) {
//     for (int_t i = 0; i <= nlp->N; i++) {
//         d_free((real_t *)nlp->lg[i]);
//         d_free((real_t *)nlp->ug[i]);
//     }
//     free((real_t **)nlp->lg);
//     free((real_t **)nlp->ug);
// }

// static void allocate_ocp_nlp_in_cost(int_t N, int_t *nx, int_t *nu, ocp_nlp_in *const nlp) {
//     ocp_nlp_ls_cost *cost = (ocp_nlp_ls_cost *)malloc(sizeof(ocp_nlp_ls_cost));
//     cost->W = (real_t **)calloc(N + 1, sizeof(*cost->W));
//     cost->y_ref = (real_t **)calloc(N + 1, sizeof(*cost->y_ref));
//     for (int_t k = 0; k <= N; k++) {
//         d_zeros(&cost->W[k], nx[k] + nu[k], nx[k] + nu[k]);
//         d_zeros(&cost->y_ref[k], nx[k] + nu[k], 1);
//     }
//     nlp->cost = cost;
// }

// static void free_ocp_nlp_in_cost(ocp_nlp_in *const nlp) {
//     for (int_t k = 0; k <= nlp->N; k++) {
//         d_free(((ocp_nlp_ls_cost *)nlp->cost)->W[k]);
//         d_free(((ocp_nlp_ls_cost *)nlp->cost)->y_ref[k]);
//     }
//     free(((ocp_nlp_ls_cost *)nlp->cost)->W);
//     free(((ocp_nlp_ls_cost *)nlp->cost)->y_ref);
//     free(nlp->cost);
// }

// static void allocate_ocp_nlp_in_sim_solver(int_t N, int_t *nx, int_t *nu, ocp_nlp_in *const nlp) {
//     nlp->sim = (sim_solver *)calloc(N, sizeof(sim_solver));
//     for (int_t i = 0; i < N; i++) {
//         int_t nx_i = nx[i];
//         int_t nu_i = nu[i];
//         nlp->sim[i].in = (sim_in *)malloc(sizeof(sim_in));
//         d_zeros(&nlp->sim[i].in->x, nx_i, 1);
//         d_zeros(&nlp->sim[i].in->u, nu_i, 1);
//         d_zeros(&nlp->sim[i].in->S_forw, nx_i, nx_i + nu_i);
//         for (int_t j = 0; j < nx_i; j++)
//             nlp->sim[i].in->S_forw[j * (nx_i + 1)] = 1.0;

//         d_zeros(&nlp->sim[i].in->S_adj, nx_i + nu_i, 1);
//         d_zeros(&nlp->sim[i].in->grad_K, nx_i, 1);

//         int_t nx_i1 = nx[i + 1];
//         nlp->sim[i].out = (sim_out *)malloc(sizeof(sim_out));
//         d_zeros(&nlp->sim[i].out->xn, nx_i1, 1);
//         d_zeros(&nlp->sim[i].out->S_forw, nx_i1, nx_i + nu_i);
//         d_zeros(&nlp->sim[i].out->grad, nx_i + nu_i, 1);
//         nlp->sim[i].out->info = (sim_info *)malloc(sizeof(sim_info));

//         nlp->sim[i].mem = NULL;
//     }
// }

// static void free_ocp_nlp_in_sim_solver(ocp_nlp_in *const nlp) {
//     for (int_t i = 0; i < nlp->N; i++) {
//         free(nlp->sim[i].in->x);
//         free(nlp->sim[i].in->u);
//         free(nlp->sim[i].in->S_forw);
//         free(nlp->sim[i].in->S_adj);
//         free(nlp->sim[i].in->grad_K);
//         free(nlp->sim[i].in);

//         free(nlp->sim[i].out->xn);
//         free(nlp->sim[i].out->S_forw);
//         free(nlp->sim[i].out->info);
//         free(nlp->sim[i].out->grad);
//         free(nlp->sim[i].out);
//     }
//     free(nlp->sim);
// }

// void allocate_ocp_nlp_in(int_t N, int_t *nx, int_t *nu, int_t *nb, int_t *nc,
//                          int_t *ng, ocp_nlp_in *const nlp) {
//     allocate_ocp_nlp_in_basic(N, nx, nu, nlp);
//     allocate_ocp_nlp_in_bounds(N, nb, nlp);
//     allocate_ocp_nlp_in_polyhedral(N, nc, nlp);
//     allocate_ocp_nlp_in_nonlinear_constraints(N, ng, nlp);
//     allocate_ocp_nlp_in_cost(N, nx, nu, nlp);
//     allocate_ocp_nlp_in_sim_solver(N, nx, nu, nlp);
// }

// void free_ocp_nlp_in(ocp_nlp_in *const nlp) {
//     free_ocp_nlp_in_basic(nlp);
//     free_ocp_nlp_in_bounds(nlp);
//     free_ocp_nlp_in_polyhedral(nlp);
//     free_ocp_nlp_in_nonlinear_constraints(nlp);
//     free_ocp_nlp_in_cost(nlp);
//     free_ocp_nlp_in_sim_solver(nlp);
// }

// void allocate_ocp_nlp_out(ocp_nlp_in *const in, ocp_nlp_out *out) {
//     int_t N = in->N;
//     out->x = (real_t **)calloc(N + 1, sizeof(*out->x));
//     out->u = (real_t **)calloc(N + 1, sizeof(*out->u));
//     out->pi = (real_t **)calloc(N, sizeof(*out->pi));
//     out->lam = (real_t **)calloc(N + 1, sizeof(*out->lam));
//     for (int_t k = 0; k < N; k++) {
//         d_zeros(&out->x[k], in->nx[k], 1);
//         d_zeros(&out->u[k], in->nu[k], 1);
//         d_zeros(&out->pi[k], in->nx[k+1], 1);
//         d_zeros(&out->lam[k], 2 * (in->nb[k] + in->nc[k] + in->ng[k]), 1);
//     }
//     d_zeros(&out->x[N], in->nx[N], 1);
//     d_zeros(&out->u[N], in->nu[N], 1);
//     d_zeros(&out->lam[N], 2 * (in->nb[N] + in->nc[N] + in->ng[N]), 1);
// }

// void free_ocp_nlp_out(int_t N, ocp_nlp_out *out) {
//     for (int_t k = 0; k < N; k++) {
//         d_free(out->x[k]);
//         d_free(out->u[k]);
//         d_free(out->pi[k]);
//         d_free(out->lam[k]);
//     }
//     d_free(out->x[N]);
//     d_free(out->u[N]);
//     d_free(out->lam[N]);
//     free(out->x);
//     free(out->u);
//     free(out->pi);
//     free(out->lam);
// }
