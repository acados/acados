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

#include "acados/ocp_nlp/ocp_nlp_sm_gn.h"
#include "acados/utils/math.h"
#include "acados/utils/types.h"

int_t ocp_nlp_sm_gn(const ocp_nlp_in *nlp_in, ocp_nlp_memory *nlp_mem, void *args_, void *memory_, void *workspace_) {
    
    ocp_nlp_sm_gn_memory *mem = (ocp_nlp_sm_gn_memory *) memory_;
    ocp_nlp_sm_gn_workspace *work = (ocp_nlp_sm_gn_workspace *) workspace_;

    const int_t N = nlp_in->N;
    const int_t *nx = nlp_in->nx;
    const int_t *nu = nlp_in->nu;
    const int_t *ng = nlp_in->ng;

    real_t **hess_l = (real_t **)nlp_mem->hess_l;
    real_t **grad_f = (real_t **)nlp_mem->grad_f;
    real_t **jac_h = (real_t **)nlp_mem->jac_h;
    real_t **jac_g = (real_t **)nlp_mem->jac_g;
    real_t **h = (real_t **)nlp_mem->h;
    real_t **g = (real_t **)nlp_mem->g;

    sim_solver *sim = mem->sim;
    ocp_nlp_ls_cost *ls_cost = mem->ls_cost;
    ocp_nlp_function *path_constraints = mem->path_constraints;

    const int_t *ny = ls_cost->ny;
    
    for (int_t i = 0; i < N; i++) {
        // Pass state and control to integrator
        for (int_t j = 0; j < nx[i]; j++) sim[i].in->x[j] = nlp_mem->x[i][j];
        for (int_t j = 0; j < nu[i]; j++) sim[i].in->u[j] = nlp_mem->x[i][j];
        sim[i].fun(sim[i].in, sim[i].out, sim[i].args, sim[i].mem, sim[i].work);

        // Sensitivities for the linearization of the system dynamics
        // TODO(rien): transition functions for changing dimensions not yet implemented!
        for (int_t j = 0; j < nx[i]; j++) {
            h[i][j] = sim[i].out->xn[j];
            for (int_t k = 0; k < nx[i] + nu[i]; k++)
                jac_h[i][k * nx[i] + j] = sim[i].out->S_forw[k * nx[i] + j];
        }
    }

    for (int_t i = 0; i <= N; i++) {
        // Sensitivities for the quadratic approximation of the objective
        for (int_t j = 0; j < nx[i]; j++) work->w[j] = nlp_mem->x[i][j];
        for (int_t j = 0; j < nu[i]; j++) work->w[nx[i] + j] = nlp_mem->u[i][j];
        // Compute residual vector F
        // ls_cost[i].fun(work->w, work->F, ...);
        for (int_t j = 0; j < ny[i]; j++) work->F[j] -= ls_cost[i].y_ref[j];
        // Compute Jacobian of residual-vector F
        //ls_cost[i].jac_fun(work->w, work->DF, ...); TODO(nielsvd): what do the arguments really mean?
        // Take transpose of DF
        for (int_t j = 0; j < nx[i] + nu[i]; j++) {
            for (int_t k = 0; k < ny[i]; k++) {
                work->DFT[k * (nx[i] + nu[i]) + j] = work->DF[j * ny[i] + k];
            }
        }
        // Compute Gauss-Newon Hessian
        for (int_t j = 0; j < (nx[i] + nu[i]) * ny[i]; j++) work->DFTW[j] = 0;
        dgemm_nn_3l(nx[i]+nu[i], ny[i], ny[i], work->DFT, nx[i]+nu[i], ls_cost[i].W, ny[i], work->DFTW, nx[i]+nu[i]);
        dgemm_nn_3l(nx[i]+nu[i], nx[i]+nu[i], ny[i], work->DFTW, nx[i]+nu[i], work->DF, ny[i], hess_l[i], nx[i]+nu[i]);
        // Compute gradient of cost
        dgemv_n_3l(nx[i]+nu[i], ny[i], work->DFTW, nx[i]+nu[i], work->F, grad_f[i]);

        // Sensitivities for the linearization of the path constraints
        path_constraints[i].fun(work->w, work->G);
        path_constraints[i].fun(work->w, work->DG);
        for (int_t j = 0; j < ng[i]; j++) {
            g[i][j] = work->G[j];
            for (int_t k = 0; k < nx[i] + nu[i]; k++) {
                jac_g[i][k * (nx[i] + nu[i]) + j] =
                    work->DG[k * (nx[i] + nu[i]) + j];
            }
        }
    }

    return 0;
}