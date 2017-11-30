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

#include <assert.h>
#include <stdlib.h>

#include "acados/sim/sim_common.h"
#include "acados/sim/sim_rk_common.h"
#include "acados/utils/math.h"
#include "acados/utils/types.h"

ocp_nlp_sm_gn_args *ocp_nlp_sm_gn_create_arguments() {
    ocp_nlp_sm_gn_args *args =
        (ocp_nlp_sm_gn_args *)malloc(sizeof(ocp_nlp_sm_gn_args));
    args->freezeSens = false;
    return args;
}

int_t ocp_nlp_sm_gn_calculate_memory_size(const ocp_nlp_sm_in *sm_in,
                                          void *args_) {
    int_t size = sizeof(ocp_nlp_sm_gn_memory);
    return size;
}

char *ocp_nlp_sm_gn_assign_memory(const ocp_nlp_sm_in *sm_in, void *args_,
                                  void **mem_, void *raw_memory) {
    ocp_nlp_sm_gn_memory **sm_memory = (ocp_nlp_sm_gn_memory **)mem_;

    // char pointer
    char *c_ptr = (char *)raw_memory;

    *sm_memory = (ocp_nlp_sm_gn_memory *)c_ptr;
    c_ptr += sizeof(ocp_nlp_sm_gn_memory);

    return c_ptr;
}

ocp_nlp_sm_gn_memory *ocp_nlp_sm_gn_create_memory(const ocp_nlp_sm_in *sm_in,
                                                  void *args_) {
    ocp_nlp_sm_gn_memory *mem;

    int_t memory_size = ocp_nlp_sm_gn_calculate_memory_size(sm_in, args_);
    void *raw_memory_ptr = malloc(memory_size);

    char *ptr_end = ocp_nlp_sm_gn_assign_memory(sm_in, args_, (void **)&mem,
                                                raw_memory_ptr);
    assert((char *)raw_memory_ptr + memory_size >= ptr_end);
    (void)ptr_end;

    return mem;
}

void size_of_workspace_elements(const ocp_nlp_sm_in *sm_in, const int_t stage,
                                int_t *size_F, int_t *size_DF, int_t *size_DFT,
                                int_t *size_DFTW, int_t *size_G,
                                int_t *size_DG) {
    const int_t *nx = sm_in->nx;
    const int_t *nu = sm_in->nu;
    const int_t *ng = sm_in->ng;

    ocp_nlp_function **cost_fun = ((ocp_nlp_ls_cost *)sm_in->cost)->fun;

    *size_F = (cost_fun[stage]->ny) * sizeof(real_t);
    *size_DF = (cost_fun[stage]->ny * (nx[stage] + nu[stage])) * sizeof(real_t);
    *size_DFT =
        ((nx[stage] + nu[stage]) * cost_fun[stage]->ny) * sizeof(real_t);
    *size_DFTW =
        ((nx[stage] + nu[stage]) * cost_fun[stage]->ny) * sizeof(real_t);
    *size_G = (ng[stage]) * sizeof(real_t);
    *size_DG = (ng[stage] * (nx[stage] + nu[stage])) * sizeof(real_t);
}

int_t ocp_nlp_sm_gn_calculate_workspace_size(const ocp_nlp_sm_in *sm_in,
                                             void *args_) {
    int_t N = sm_in->N;

    int_t size = sizeof(ocp_nlp_sm_gn_workspace);

    size += (N + 1) * sizeof(real_t *);  // F
    size += (N + 1) * sizeof(real_t *);  // DF
    size += (N + 1) * sizeof(real_t *);  // DFT
    size += (N + 1) * sizeof(real_t *);  // DFTW
    size += (N + 1) * sizeof(real_t *);  // G
    size += (N + 1) * sizeof(real_t *);  // DG

    for (int_t i = 0; i <= N; i++) {
        int_t size_F, size_DF, size_DFT, size_DFTW, size_G, size_DG;
        size_of_workspace_elements(sm_in, i, &size_F, &size_DF, &size_DFT,
                                   &size_DFTW, &size_G, &size_DG);

        size += size_F + size_DF + size_DFT + size_DFTW + size_G + size_DG;
    }

    return size;
}

char *ocp_nlp_sm_gn_assign_workspace(const ocp_nlp_sm_in *sm_in, void *args_,
                                     void **work_, void *raw_memory) {
    int_t N = sm_in->N;

    ocp_nlp_sm_gn_workspace **sm_workspace = (ocp_nlp_sm_gn_workspace **)work_;
    char *c_ptr = (char *)raw_memory;

    *sm_workspace = (ocp_nlp_sm_gn_workspace *)c_ptr;
    c_ptr += sizeof(ocp_nlp_sm_gn_workspace);

    (*sm_workspace)->F = (real_t **)c_ptr;
    c_ptr += (N + 1) * sizeof(real_t *);

    (*sm_workspace)->DF = (real_t **)c_ptr;
    c_ptr += (N + 1) * sizeof(real_t *);

    (*sm_workspace)->DFT = (real_t **)c_ptr;
    c_ptr += (N + 1) * sizeof(real_t *);

    (*sm_workspace)->DFTW = (real_t **)c_ptr;
    c_ptr += (N + 1) * sizeof(real_t *);

    (*sm_workspace)->G = (real_t **)c_ptr;
    c_ptr += (N + 1) * sizeof(real_t *);

    (*sm_workspace)->DG = (real_t **)c_ptr;
    c_ptr += (N + 1) * sizeof(real_t *);

    for (int_t i = 0; i <= N; i++) {
        int_t size_F, size_DF, size_DFT, size_DFTW, size_G, size_DG;
        size_of_workspace_elements(sm_in, i, &size_F, &size_DF, &size_DFT,
                                   &size_DFTW, &size_G, &size_DG);

        (*sm_workspace)->F[i] = (real_t *)c_ptr;
        c_ptr += size_F;

        (*sm_workspace)->DF[i] = (real_t *)c_ptr;
        c_ptr += size_DF;

        (*sm_workspace)->DFT[i] = (real_t *)c_ptr;
        c_ptr += size_DFT;

        (*sm_workspace)->DFTW[i] = (real_t *)c_ptr;
        c_ptr += size_DFTW;

        (*sm_workspace)->G[i] = (real_t *)c_ptr;
        c_ptr += size_G;

        (*sm_workspace)->DG[i] = (real_t *)c_ptr;
        c_ptr += size_DG;
    }

    return c_ptr;
}

ocp_nlp_sm_gn_workspace *ocp_nlp_sm_gn_create_workspace(
    const ocp_nlp_sm_in *sm_in, void *args_) {
    ocp_nlp_sm_gn_workspace *work;

    int_t workspace_size = ocp_nlp_sm_gn_calculate_workspace_size(sm_in, args_);
    void *raw_memory_ptr = malloc(workspace_size);

    char *ptr_end = ocp_nlp_sm_gn_assign_workspace(sm_in, args_, (void **)&work,
                                                   raw_memory_ptr);
    assert((char *)raw_memory_ptr + workspace_size >= ptr_end);
    (void)ptr_end;

    return work;
}

int_t ocp_nlp_sm_gn(const ocp_nlp_sm_in *sm_in, ocp_nlp_sm_out *sm_out,
                    void *args_, void *memory_, void *workspace_) {
    ocp_nlp_sm_gn_workspace *work = (ocp_nlp_sm_gn_workspace *)workspace_;
    ocp_nlp_sm_gn_memory *mem = memory_;

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

    sim_solver **sim = sm_in->sim;
    ocp_nlp_ls_cost *ls_cost = (ocp_nlp_ls_cost *)sm_in->cost;
    ocp_nlp_function **path_constraints = sm_in->path_constraints;

    for (int_t i = 0; i < N; i++) {
        // Adjoint-based gradient correction (used for)
        // TODO(nielsvd): create new sensitivity methods for inexact newton methods
        sim_rk_opts *sim_opts = (sim_rk_opts *)sim[i]->args;
        if (mem->inexact_init) {
            if (sim_opts->scheme.type != exact) {
                sim[i]->in->sens_adj = true;
                sim_opts->scheme.freeze = sm_in->freezeSens;
                for (int_t j = 0; j < nx[i + 1]; j++)
                    sim[i]->in->S_adj[j] = -sm_in->pi[i][j];
            }
        } else {
            sim[i]->in->sens_adj = false;
        }

        // Pass state and control to integrator
        for (int_t j = 0; j < nx[i]; j++) sim[i]->in->x[j] = sm_in->x[i][j];
        for (int_t j = 0; j < nu[i]; j++) sim[i]->in->u[j] = sm_in->u[i][j];
        sim[i]->fun(sim[i]->in, sim[i]->out, sim[i]->args, sim[i]->mem, sim[i]->work);

        // Sensitivities for the linearization of the system dynamics
        // TODO(rien): transition functions for changing dimensions not yet
        // implemented!
        for (int_t j = 0; j < nx[i]; j++) {
            h[i][j] = sim[i]->out->xn[j];
            for (int_t k = 0; k < nx[i] + nu[i]; k++)
                jac_h[i][k * nx[i] + j] = sim[i]->out->S_forw[k * nx[i] + j];
        }
    }

    for (int_t i = 0; i <= N; i++) {
        // Least squares cost for shooting node i
        const int_t ny = ls_cost->fun[i]->ny;
        casadi_wrapper_in *ls_in = ls_cost->fun[i]->in;
        casadi_wrapper_out *ls_out = ls_cost->fun[i]->out;
        casadi_wrapper_args *ls_args = ls_cost->fun[i]->args;
        casadi_wrapper_workspace *ls_work = ls_cost->fun[i]->work;

        // Sensitivities for the quadratic approximation of the objective
        // Compute residual vector F and its Jacobian
        casadi_wrapper(ls_in, ls_out, ls_args, ls_work);
        for (int_t j = 0; j < ny; j++) work->F[i][j] -= ls_cost->y_ref[i][j];
        // Take transpose of DF
        for (int_t j = 0; j < nx[i] + nu[i]; j++) {
            for (int_t k = 0; k < ny; k++)
                work->DFT[i][k * (nx[i] + nu[i]) + j] = work->DF[i][j * ny + k];
        }

        // Compute Gauss-Newton Hessian
        for (int_t j = 0; j < (nx[i] + nu[i]) * ny; j++) work->DFTW[i][j] = 0;
        dgemm_nn_3l(nx[i] + nu[i], ny, ny, work->DFT[i], nx[i] + nu[i],
                    (real_t *)ls_cost->W[i], ny, work->DFTW[i], nx[i] + nu[i]);
        dgemm_nn_3l(nx[i] + nu[i], nx[i] + nu[i], ny, work->DFTW[i],
                    nx[i] + nu[i], work->DF[i], ny, hess_l[i], nx[i] + nu[i]);
        // Compute gradient of cost
        for (int_t j = 0; j < (nx[i] + nu[i]); j++) grad_f[i][j] = 0;
        dgemv_n_3l(nx[i] + nu[i], ny, work->DFTW[i], nx[i] + nu[i], work->F[i],
                   grad_f[i]);

        if (sm_in->ng[i] > 0) {
            // Path constraints for shooting node i
            casadi_wrapper_in *pc_in = path_constraints[i]->in;
            casadi_wrapper_out *pc_out = path_constraints[i]->out;
            casadi_wrapper_args *pc_args = path_constraints[i]->args;
            casadi_wrapper_workspace *pc_work = path_constraints[i]->work;
            // Sensitivities for the linearization of the path constraints
            casadi_wrapper(pc_in, pc_out, pc_args, pc_work);
            for (int_t j = 0; j < ng[i]; j++) {
                g[i][j] = work->G[i][j];
                for (int_t k = 0; k < nx[i] + nu[i]; k++)
                    jac_g[i][k * ng[i] + j] = work->DG[i][k * ng[i] + j];
            }
        }
    }

    // Adjoint-based gradient correction
    // TODO(nielsvd): create new sensitivity methods for inexact newton methods
    if (mem->inexact_init) {
        for (int_t i = 0; i < N; i++) {
            sim_rk_opts *sim_opts = (sim_rk_opts *)sim[i]->args;
            if (sim_opts->scheme.type != exact) {
                for (int_t j = 0; j < nx[i] + nu[i]; j++) {
                    grad_f[i][j] += sim[i]->out->grad[j];
                }
            }
        }
    } else {
        mem->inexact_init = true;
    }

    return 0;
}

void ocp_nlp_sm_gn_initialize(const ocp_nlp_sm_in *sm_in, void *args_,
                              void **mem_, void **work_) {
    ocp_nlp_sm_gn_memory **mem = (ocp_nlp_sm_gn_memory **)mem_;
    ocp_nlp_sm_gn_workspace **work = (ocp_nlp_sm_gn_workspace **)work_;

    *mem = ocp_nlp_sm_gn_create_memory(sm_in, args_);
    *work = ocp_nlp_sm_gn_create_workspace(sm_in, args_);

    (*mem)->inexact_init = false;

    int_t N = sm_in->N;
    ocp_nlp_ls_cost *ls_cost = (ocp_nlp_ls_cost *)sm_in->cost;
    ocp_nlp_function **path_constraints = sm_in->path_constraints;

    for (int_t i = 0; i <= N; i++) {
        // assign correct pointers to ls_cost-array
        ls_cost->fun[i]->in->x = sm_in->x[i];
        ls_cost->fun[i]->in->u = sm_in->u[i];
        ls_cost->fun[i]->in->p = NULL;  // TODO(nielsvd): support for parameters
        ls_cost->fun[i]->out->y = (*work)->F[i];
        ls_cost->fun[i]->out->jac_y = (*work)->DF[i];
        ls_cost->fun[i]->in->compute_jac = true;
        ls_cost->fun[i]->in->compute_hess = false;

        if (sm_in->ng[i] > 0) {
            // assign correct pointers to path_constraints-array
            path_constraints[i]->in->x = sm_in->x[i];
            path_constraints[i]->in->u = sm_in->u[i];
            path_constraints[i]->in->p = NULL;  // TODO(nielsvd): support for parameters
            path_constraints[i]->out->y = (*work)->G[i];
            path_constraints[i]->out->jac_y = (*work)->DG[i];
            path_constraints[i]->in->compute_jac = true;
            path_constraints[i]->in->compute_hess = false;
        }
    }
}

void ocp_nlp_sm_gn_destroy(void *mem_, void *work_) {
    ocp_nlp_sm_gn_memory *mem = (ocp_nlp_sm_gn_memory *)mem_;
    ocp_nlp_sm_gn_workspace *work = (ocp_nlp_sm_gn_workspace *)work_;

    free(mem);
    free(work);
}
