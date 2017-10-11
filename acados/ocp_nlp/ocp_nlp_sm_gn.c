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

// TODO(nielsvd): only perform assert in debug mode?
#include <assert.h>

#include <stdlib.h>

ocp_nlp_sm_gn_args *ocp_nlp_sm_gn_create_arguments() {
    
    ocp_nlp_sm_gn_args *args = (ocp_nlp_sm_gn_args *) malloc(sizeof(ocp_nlp_sm_gn_args));
    args->dummy = 0;
    return args;
}

int_t ocp_nlp_sm_gn_calculate_memory_size(const ocp_nlp_sm_in *sm_in, void *args_) {
    
    int_t size = sizeof(ocp_nlp_sm_gn_memory);
    return size;
}

char *ocp_nlp_sm_gn_assign_memory(const ocp_nlp_sm_in *sm_in, void *args_, void **mem_, void *raw_memory) {
    
    ocp_nlp_sm_gn_memory **sm_memory = (ocp_nlp_sm_gn_memory **) mem_;

    // char pointer
    char *c_ptr = (char *)raw_memory;

    *sm_memory = (ocp_nlp_sm_gn_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_sm_gn_memory);

    return c_ptr;
}

ocp_nlp_sm_gn_memory *ocp_nlp_sm_gn_create_memory(const ocp_nlp_sm_in *sm_in, void *args_) {

    ocp_nlp_sm_gn_memory *mem;

    int_t memory_size = ocp_nlp_sm_gn_calculate_memory_size(sm_in, args_);
    void *raw_memory_ptr = malloc(memory_size);

    char *ptr_end = ocp_nlp_sm_gn_assign_memory(sm_in, args_, (void **) &mem, raw_memory_ptr);
    assert((char *)raw_memory_ptr + memory_size >= ptr_end); (void)ptr_end;

    return mem;
}

void size_of_workspace_elements(const ocp_nlp_sm_in *sm_in, int_t *size_w,
                                int_t *size_F, int_t *size_DF, int_t *size_DFT,
                                int_t *size_DFTW, int_t *size_G,
                                int_t *size_DG) {
    int_t N = sm_in->N;
    int_t *nx = (int_t *)sm_in->nx;
    int_t *nu = (int_t *)sm_in->nu;
    int_t *ng = (int_t *)sm_in->ng;

    ocp_nlp_ls_cost *cost = (ocp_nlp_ls_cost *)sm_in->cost;

    *size_w = 0;
    *size_F = 0;
    *size_DF = 0;
    *size_DFT = 0;
    *size_DFTW = 0;
    *size_G = 0;
    *size_DG = 0;

    for (int_t i = 0; i <= N; i++) {
        // w    -- max_i (nx[i]+nu[i])
        if (*size_w < nx[i] + nu[i]) *size_w = nx[i] + nu[i];
        // F    -- max_i ny[i]
        if (*size_F < cost[i].ny) *size_F = cost[i].ny;
        // DF   -- max_i ny[i]*(nx[i]+nu[i])
        if (*size_DF < cost[i].ny * (nx[i] + nu[i]))
            *size_DF = cost[i].ny * (nx[i] + nu[i]);
        // DFT  -- max_i (nx[i]*nu[i])*ny[i]
        if (*size_DFT < (nx[i] * nu[i]) * cost[i].ny)
            *size_DFT = (nx[i] * nu[i]) * cost[i].ny;
        // DFTW -- max_i (nx[i]+nu[i])*ny[i]
        if (*size_DFTW < (nx[i] + nu[i]) * cost[i].ny)
            *size_DFTW = (nx[i] + nu[i]) * cost[i].ny;
        // G    -- max_i ng[i]
        if (*size_G < ng[i]) *size_G = ng[i];
        // DG   -- max_i ng[i]*(nx[i]+nu[i])
        if (*size_DG < ng[i] * (nx[i] + nu[i]))
            *size_DG = ng[i] * (nx[i] + nu[i]);
    }

    // Size in bytes
    *size_w *= sizeof(real_t);
    *size_F *= sizeof(real_t);
    *size_DF *= sizeof(real_t);
    *size_DFT *= sizeof(real_t);
    *size_DFTW *= sizeof(real_t);
    *size_G *= sizeof(real_t);
    *size_DG *= sizeof(real_t);
}

int_t ocp_nlp_sm_gn_calculate_workspace_size(const ocp_nlp_sm_in *sm_in, void *args_) {

    int_t size = sizeof(ocp_nlp_sm_gn_workspace);

    int_t size_w, size_F, size_DF, size_DFT, size_DFTW, size_G, size_DG;
    size_of_workspace_elements(sm_in, &size_w, &size_F, &size_DF, &size_DFT,
                               &size_DFTW, &size_G, &size_DG);

    size += size_w + size_F + size_DF + size_DFT + size_DFTW + size_G + size_DG;

    return size;
}

char *ocp_nlp_sm_gn_assign_workspace(const ocp_nlp_sm_in *sm_in, void *args_, void **work_, void *raw_memory) {
    
    ocp_nlp_sm_gn_workspace **sm_workspace = (ocp_nlp_sm_gn_workspace **)work_;
    char *c_ptr = (char *)raw_memory;

    *sm_workspace = (ocp_nlp_sm_gn_workspace *)c_ptr;
    c_ptr += sizeof(ocp_nlp_sm_gn_workspace);

    int_t size_w, size_F, size_DF, size_DFT, size_DFTW, size_G, size_DG;
    size_of_workspace_elements(sm_in, &size_w, &size_F, &size_DF, &size_DFT,
                               &size_DFTW, &size_G, &size_DG);

    (*sm_workspace)->w = (real_t *) c_ptr;
    c_ptr += size_w;

    (*sm_workspace)->F = (real_t *)c_ptr;
    c_ptr += size_F;

    (*sm_workspace)->DF = (real_t *)c_ptr;
    c_ptr += size_DF;

    (*sm_workspace)->DFT = (real_t *)c_ptr;
    c_ptr += size_DFT;

    (*sm_workspace)->DFTW = (real_t *)c_ptr;
    c_ptr += size_DFTW;

    (*sm_workspace)->G = (real_t *)c_ptr;
    c_ptr += size_G;

    (*sm_workspace)->DG = (real_t *)c_ptr;
    c_ptr += size_DG;

    return c_ptr;
}

ocp_nlp_sm_gn_workspace *ocp_nlp_sm_gn_create_workspace(const ocp_nlp_sm_in *sm_in, void *args_) {
    
    ocp_nlp_sm_gn_workspace *work;

    int_t workspace_size = ocp_nlp_sm_gn_calculate_workspace_size(sm_in, args_);
    void *raw_memory_ptr = malloc(workspace_size);
    
    char *ptr_end = ocp_nlp_sm_gn_assign_workspace(sm_in, args_, (void **)&work, raw_memory_ptr);
    assert((char *)raw_memory_ptr + workspace_size >= ptr_end); (void)ptr_end;
    
    return work;
}

int_t ocp_nlp_sm_gn(const ocp_nlp_sm_in *sm_in, ocp_nlp_sm_out *sm_out, void *args_, void *memory_, void *workspace_) {
    
    ocp_nlp_sm_gn_workspace *work = (ocp_nlp_sm_gn_workspace *) workspace_;

    const int_t N = sm_in->N;
    const int_t *nx = sm_in->nx;
    const int_t *nu = sm_in->nu;
    const int_t *ng = sm_in->ng;

    real_t **hess_l = (real_t **)sm_out->hess_l;
    real_t **grad_f = (real_t **)sm_out->grad_f;
    real_t **jac_h = (real_t **)sm_out->jac_h;
    real_t **jac_g = (real_t **)sm_out->jac_g;
    real_t **h = (real_t **)sm_out->h;
    real_t **g = (real_t **)sm_out->g;

    sim_solver *sim = sm_in->sim;
    ocp_nlp_ls_cost *ls_cost = (ocp_nlp_ls_cost *) sm_in->cost;
    ocp_nlp_function *path_constraints = sm_in->path_constraints;

    for (int_t i = 0; i < N; i++) {
        // Pass state and control to integrator
        for (int_t j = 0; j < nx[i]; j++) sim[i].in->x[j] = sm_in->x[i][j];
        for (int_t j = 0; j < nu[i]; j++) sim[i].in->u[j] = sm_in->u[i][j];
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
        const int_t ny = ls_cost[i].ny;

        // Sensitivities for the quadratic approximation of the objective
        for (int_t j = 0; j < nx[i]; j++) work->w[j] = sm_in->x[i][j];
        for (int_t j = 0; j < nu[i]; j++) work->w[nx[i] + j] = sm_in->u[i][j];
        // Compute residual vector F
#pragma message "TODO(nielsvd): Add call to fun."
        // ls_cost[i].fun(work->w, work->F, ...);
        for (int_t j = 0; j < ny; j++) work->F[j] -= ls_cost[i].y_ref[j];
        // Compute Jacobian of residual-vector F
#pragma message "TODO(nielsvd): Add call to jac_fun."
        //ls_cost[i].jac_fun(work->w, work->DF, ...); TODO(nielsvd): what do the arguments really mean?
        // Take transpose of DF
        for (int_t j = 0; j < nx[i] + nu[i]; j++) {
            for (int_t k = 0; k < ny; k++) {
                work->DFT[k * (nx[i] + nu[i]) + j] = work->DF[j * ny + k];
            }
        }
        // Compute Gauss-Newon Hessian
        for (int_t j = 0; j < (nx[i] + nu[i]) * ny; j++) work->DFTW[j] = 0;
        dgemm_nn_3l(nx[i]+nu[i], ny, ny, work->DFT, nx[i]+nu[i], ls_cost[i].W, ny, work->DFTW, nx[i]+nu[i]);
        dgemm_nn_3l(nx[i]+nu[i], nx[i]+nu[i], ny, work->DFTW, nx[i]+nu[i], work->DF, ny, hess_l[i], nx[i]+nu[i]);
        // Compute gradient of cost
        dgemv_n_3l(nx[i]+nu[i], ny, work->DFTW, nx[i]+nu[i], work->F, grad_f[i]);

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

void ocp_nlp_sm_gn_initialize(const ocp_nlp_sm_in *sm_in, void *args_, void **mem, void **work) {
    *mem = ocp_nlp_sm_gn_create_memory(sm_in, args_);
    *work = ocp_nlp_sm_gn_create_workspace(sm_in, args_);
}

void ocp_nlp_sm_gn_destroy(void *mem, void *work) {
    // TODO(nielsvd): replace dummy commands once interface completed
    (void)mem;
    (void)work;
}