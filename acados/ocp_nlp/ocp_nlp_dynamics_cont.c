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

#include "acados/ocp_nlp/ocp_nlp_dynamics_cont.h"
#include "acados/ocp_nlp/ocp_nlp_common.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

// blasfeo
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_blas.h"
// acados
#include "acados/utils/mem.h"



/************************************************
 * dims
 ************************************************/

int ocp_nlp_dynamics_cont_dims_calculate_size(void *config_)
{
    ocp_nlp_dynamics_config *dyn_config = (ocp_nlp_dynamics_config *) config_;
    sim_config *sim_sol_config = (sim_config *) dyn_config->sim_solver;
    int size = 0;

    size += sizeof(ocp_nlp_dynamics_cont_dims);

    size += sim_sol_config->dims_calculate_size(sim_sol_config);

    return size;
}


void *ocp_nlp_dynamics_cont_dims_assign(void *config_, void *raw_memory)
{
    ocp_nlp_dynamics_config *dyn_config = (ocp_nlp_dynamics_config *) config_;
    sim_config *sim_sol_config = (sim_config *) dyn_config->sim_solver;

    char *c_ptr = (char *) raw_memory;

    ocp_nlp_dynamics_cont_dims *dims = (ocp_nlp_dynamics_cont_dims *) c_ptr;
    c_ptr += sizeof(ocp_nlp_dynamics_cont_dims);

    dims->sim = sim_sol_config->dims_assign(sim_sol_config, c_ptr);

    c_ptr += sim_sol_config->dims_calculate_size(sim_sol_config);

    assert((char *) raw_memory + ocp_nlp_dynamics_cont_dims_calculate_size(config_) >= c_ptr);

    return dims;
}



void ocp_nlp_dynamics_cont_dims_initialize(void *config_, void *dims_, int nx, int nu, int nx1,
                                           int nu1, int nz)
{
    ocp_nlp_dynamics_cont_dims *dims = dims_;

    dims->nx = nx;
    dims->nz = nz;
    dims->nu = nu;
    dims->nx1 = nx1;
    dims->nu1 = nu1;

    ocp_nlp_dynamics_config *dyn_config = (ocp_nlp_dynamics_config *) config_;
    sim_config *sim_config_ = (sim_config *) dyn_config->sim_solver;

    sim_config_->dims_set(sim_config_, dims->sim, "nx", &nx);
    sim_config_->dims_set(sim_config_, dims->sim, "nu", &nu);
    sim_config_->dims_set(sim_config_, dims->sim, "nz", &nz);

    return;
}


// setters
static void ocp_nlp_dynamics_cont_set_nx(void *config_, void *dims_, int *nx)
{
    ocp_nlp_dynamics_cont_dims *dims = (ocp_nlp_dynamics_cont_dims *) dims_;
    dims->nx = *nx;

    ocp_nlp_dynamics_config *dyn_config = (ocp_nlp_dynamics_config *) config_;
    sim_config *sim_config_ = (sim_config *) dyn_config->sim_solver;

    sim_config_->dims_set(sim_config_, dims->sim, "nx", nx);
}

static void ocp_nlp_dynamics_cont_set_nx1(void *config_, void *dims_, int *nx1)
{
    ocp_nlp_dynamics_cont_dims *dims = (ocp_nlp_dynamics_cont_dims *) dims_;
    dims->nx1 = *nx1;
}

static void ocp_nlp_dynamics_cont_set_nz(void *config_, void *dims_, int *nz)
{
    ocp_nlp_dynamics_cont_dims *dims = (ocp_nlp_dynamics_cont_dims *) dims_;
    dims->nz = *nz;

    ocp_nlp_dynamics_config *dyn_config = (ocp_nlp_dynamics_config *) config_;
    sim_config *sim_config_ = (sim_config *) dyn_config->sim_solver;

    sim_config_->dims_set(sim_config_, dims->sim, "nz", nz);
}

static void ocp_nlp_dynamics_cont_set_nu(void *config_, void *dims_, int *nu)
{
    ocp_nlp_dynamics_cont_dims *dims = (ocp_nlp_dynamics_cont_dims *) dims_;
    dims->nu = *nu;

    ocp_nlp_dynamics_config *dyn_config = (ocp_nlp_dynamics_config *) config_;
    sim_config *sim_config_ = (sim_config *) dyn_config->sim_solver;

    sim_config_->dims_set(sim_config_, dims->sim, "nu", nu);
}

static void ocp_nlp_dynamics_cont_set_nu1(void *config_, void *dims_, int *nu1)
{
    ocp_nlp_dynamics_cont_dims *dims = (ocp_nlp_dynamics_cont_dims *) dims_;
    dims->nu1 = *nu1;
}

void ocp_nlp_dynamics_cont_dims_set(void *config_, void *dims_, const char *field, int* value)
{
    if (!strcmp(field, "nx"))
    {
        ocp_nlp_dynamics_cont_set_nx(config_, dims_, value);
    }
    else if (!strcmp(field, "nx1"))
    {
        ocp_nlp_dynamics_cont_set_nx1(config_, dims_, value);
    }
    else if (!strcmp(field, "nz"))
    {
        ocp_nlp_dynamics_cont_set_nz(config_, dims_, value);
    }
    else if (!strcmp(field, "nu"))
    {
        ocp_nlp_dynamics_cont_set_nu(config_, dims_, value);
    }
    else if (!strcmp(field, "nu1"))
    {
        ocp_nlp_dynamics_cont_set_nu1(config_, dims_, value);
    }
    else
    {
        // set GNSF dims just within integrator module
        ocp_nlp_dynamics_config *dyn_config = (ocp_nlp_dynamics_config *) config_;
        ocp_nlp_dynamics_cont_dims *dims = (ocp_nlp_dynamics_cont_dims *) dims_;
        sim_config *sim_config_ = (sim_config *) dyn_config->sim_solver;

        sim_config_->dims_set(sim_config_, dims->sim, field, value);
    }
}

/************************************************
 * options
 ************************************************/

int ocp_nlp_dynamics_cont_opts_calculate_size(void *config_, void *dims_)
{
    ocp_nlp_dynamics_config *config = config_;
    ocp_nlp_dynamics_cont_dims *dims = dims_;

    int size = 0;

    size += sizeof(ocp_nlp_dynamics_cont_opts);

    size += config->sim_solver->opts_calculate_size(config->sim_solver, dims->sim);

    return size;
}



void *ocp_nlp_dynamics_cont_opts_assign(void *config_, void *dims_, void *raw_memory)
{
    ocp_nlp_dynamics_config *config = config_;
    ocp_nlp_dynamics_cont_dims *dims = dims_;

    char *c_ptr = (char *) raw_memory;

    ocp_nlp_dynamics_cont_opts *opts = (ocp_nlp_dynamics_cont_opts *) c_ptr;
    c_ptr += sizeof(ocp_nlp_dynamics_cont_opts);

    opts->sim_solver = config->sim_solver->opts_assign(config->sim_solver, dims->sim, c_ptr);
    c_ptr += config->sim_solver->opts_calculate_size(config->sim_solver, dims->sim);

    assert((char *) raw_memory + ocp_nlp_dynamics_cont_opts_calculate_size(config, dims) >= c_ptr);

    return opts;
}



void ocp_nlp_dynamics_cont_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dynamics_config *config = config_;
    ocp_nlp_dynamics_cont_dims *dims = dims_;
    ocp_nlp_dynamics_cont_opts *opts = opts_;

    opts->compute_adj = 1;

    config->sim_solver->opts_initialize_default(config->sim_solver, dims->sim, opts->sim_solver);

    return;
}



void ocp_nlp_dynamics_cont_opts_update(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dynamics_config *config = config_;
    ocp_nlp_dynamics_cont_dims *dims = dims_;
    ocp_nlp_dynamics_cont_opts *opts = opts_;

    config->sim_solver->opts_update(config->sim_solver, dims->sim, opts->sim_solver);

    return;
}




int ocp_nlp_dynamics_cont_opts_set(void *config_, void *opts_, const char *field, void* value)
{
    ocp_nlp_dynamics_config *config = config_;
    ocp_nlp_dynamics_cont_opts *opts = opts_;
    sim_config *sim_config_ = config->sim_solver;

    return sim_config_->opts_set(sim_config_, opts->sim_solver, field, value);

}



/************************************************
 * memory
 ************************************************/

int ocp_nlp_dynamics_cont_memory_calculate_size(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dynamics_config *config = config_;
    ocp_nlp_dynamics_cont_dims *dims = dims_;
    ocp_nlp_dynamics_cont_opts *opts = opts_;

    // extract dims
    int nx = dims->nx;
    int nu = dims->nu;
    int nx1 = dims->nx1;

    int size = 0;

    size += sizeof(ocp_nlp_dynamics_cont_memory);

    size += 1 * blasfeo_memsize_dvec(nu + nx + nx1);  // adj
    size += 1 * blasfeo_memsize_dmat(nu+nx, nu+nx);   // hes
    size += 1 * blasfeo_memsize_dvec(nx1);            // fun

    size +=
        config->sim_solver->memory_calculate_size(config->sim_solver, dims->sim, opts->sim_solver);

    size += 1*64;  // blasfeo_mem align

    return size;
}



void *ocp_nlp_dynamics_cont_memory_assign(void *config_, void *dims_, void *opts_, void *raw_memory)
{
    ocp_nlp_dynamics_config *config = config_;
    ocp_nlp_dynamics_cont_dims *dims = dims_;
    ocp_nlp_dynamics_cont_opts *opts = opts_;

    char *c_ptr = (char *) raw_memory;

    // extract dims
    int nx = dims->nx;
    int nu = dims->nu;
    int nx1 = dims->nx1;

    // struct
    ocp_nlp_dynamics_cont_memory *memory = (ocp_nlp_dynamics_cont_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_dynamics_cont_memory);

    // sim_solver
    memory->sim_solver =
        config->sim_solver->memory_assign(config->sim_solver, dims->sim, opts->sim_solver, c_ptr);
    c_ptr +=
        config->sim_solver->memory_calculate_size(config->sim_solver, dims->sim, opts->sim_solver);

    // blasfeo_mem align
    align_char_to(64, &c_ptr);

    // hes
    assign_and_advance_blasfeo_dmat_mem(nu+nx, nu+nx, &memory->hes, &c_ptr);

    // adj
    assign_and_advance_blasfeo_dvec_mem(nu + nx + nx1, &memory->adj, &c_ptr);

    // fun
    assign_and_advance_blasfeo_dvec_mem(nx1, &memory->fun, &c_ptr);

    assert((char *) raw_memory +
               ocp_nlp_dynamics_cont_memory_calculate_size(config_, dims, opts_) >=
           c_ptr);

    return memory;
}



struct blasfeo_dvec *ocp_nlp_dynamics_cont_memory_get_fun_ptr(void *memory_)
{
    ocp_nlp_dynamics_cont_memory *memory = memory_;

    return &memory->fun;
}



struct blasfeo_dvec *ocp_nlp_dynamics_cont_memory_get_adj_ptr(void *memory_)
{
    ocp_nlp_dynamics_cont_memory *memory = memory_;

    return &memory->adj;
}



void ocp_nlp_dynamics_cont_memory_set_ux_ptr(struct blasfeo_dvec *ux, void *memory_)
{
    ocp_nlp_dynamics_cont_memory *memory = memory_;

    memory->ux = ux;

    return;
}



void ocp_nlp_dynamics_cont_memory_set_ux1_ptr(struct blasfeo_dvec *ux1, void *memory_)
{
    ocp_nlp_dynamics_cont_memory *memory = memory_;

    memory->ux1 = ux1;

    return;
}



void ocp_nlp_dynamics_cont_memory_set_pi_ptr(struct blasfeo_dvec *pi, void *memory_)
{
    ocp_nlp_dynamics_cont_memory *memory = memory_;

    memory->pi = pi;

    return;
}



void ocp_nlp_dynamics_cont_memory_set_BAbt_ptr(struct blasfeo_dmat *BAbt, void *memory_)
{
    ocp_nlp_dynamics_cont_memory *memory = memory_;

    memory->BAbt = BAbt;

    return;
}

void ocp_nlp_dynamics_cont_memory_set_RSQrq_ptr(struct blasfeo_dmat *RSQrq, void *memory_)
{
    ocp_nlp_dynamics_cont_memory *memory = memory_;

    memory->RSQrq = RSQrq;

    return;
}

void ocp_nlp_dynamics_cont_memory_set_z_ptr(struct blasfeo_dvec *z, void *memory_)
{
    ocp_nlp_dynamics_cont_memory *memory = memory_;

    memory->z = z;

    return;
}

/************************************************
 * workspace
 ************************************************/

int ocp_nlp_dynamics_cont_workspace_calculate_size(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dynamics_config *config = config_;
    ocp_nlp_dynamics_cont_dims *dims = dims_;
    ocp_nlp_dynamics_cont_opts *opts = opts_;

    int size = 0;

    size += sizeof(ocp_nlp_dynamics_cont_workspace);

    size += sim_in_calculate_size(config->sim_solver, dims->sim);
    size += sim_out_calculate_size(config->sim_solver, dims->sim);
    size += config->sim_solver->workspace_calculate_size(config->sim_solver, dims->sim,
                                                         opts->sim_solver);

    return size;
}



static void ocp_nlp_dynamics_cont_cast_workspace(void *config_, void *dims_, void *opts_,
                                                 void *work_)
{
    ocp_nlp_dynamics_config *config = config_;
    ocp_nlp_dynamics_cont_dims *dims = dims_;
    ocp_nlp_dynamics_cont_opts *opts = opts_;
    ocp_nlp_dynamics_cont_workspace *work = work_;

    char *c_ptr = (char *) work_;
    c_ptr += sizeof(ocp_nlp_dynamics_cont_workspace);

    // sim in
    work->sim_in = sim_in_assign(config->sim_solver, dims->sim, c_ptr);
    c_ptr += sim_in_calculate_size(config->sim_solver, dims->sim);
    // sim out
    work->sim_out = sim_out_assign(config->sim_solver, dims->sim, c_ptr);
    c_ptr += sim_out_calculate_size(config->sim_solver, dims->sim);
    // workspace
    work->sim_solver = c_ptr;
    c_ptr += config->sim_solver->workspace_calculate_size(config->sim_solver, dims->sim,
                                                          opts->sim_solver);

    assert((char *) work + ocp_nlp_dynamics_cont_workspace_calculate_size(config, dims, opts) >=
           c_ptr);

    return;
}



/************************************************
 * model
 ************************************************/

int ocp_nlp_dynamics_cont_model_calculate_size(void *config_, void *dims_)
{
    ocp_nlp_dynamics_config *config = config_;
    ocp_nlp_dynamics_cont_dims *dims = dims_;

    // extract dims
    // int nx = dims->nx;
    // int nu = dims->nu;

    int size = 0;

    size += sizeof(ocp_nlp_dynamics_cont_model);

    size += config->sim_solver->model_calculate_size(config->sim_solver, dims->sim);

    return size;
}



void *ocp_nlp_dynamics_cont_model_assign(void *config_, void *dims_, void *raw_memory)
{
    ocp_nlp_dynamics_config *config = config_;
    ocp_nlp_dynamics_cont_dims *dims = dims_;

    char *c_ptr = (char *) raw_memory;

    // extract dims
    // int nx = dims->nx;
    // int nu = dims->nu;

    // struct
    ocp_nlp_dynamics_cont_model *model = (ocp_nlp_dynamics_cont_model *) c_ptr;
    c_ptr += sizeof(ocp_nlp_dynamics_cont_model);

    model->sim_model = config->sim_solver->model_assign(config->sim_solver, dims->sim, c_ptr);
    c_ptr += config->sim_solver->model_calculate_size(config->sim_solver, dims->sim);

    assert((char *) raw_memory + ocp_nlp_dynamics_cont_model_calculate_size(config, dims) >= c_ptr);

    return model;
}



void ocp_nlp_dynamics_cont_model_set_T(double T, void *model_)
{
    ocp_nlp_dynamics_cont_model *model = model_;

    model->T = T;

    return;
}



/************************************************
 * functions
 ************************************************/

void ocp_nlp_dynamics_cont_initialize(void *config_, void *dims_, void *model_, void *opts_,
                                      void *mem_, void *work_)
{
    return;
}



void ocp_nlp_dynamics_cont_update_qp_matrices(void *config_, void *dims_, void *model_, void *opts_,
                                              void *mem_, void *work_)
{
    ocp_nlp_dynamics_cont_cast_workspace(config_, dims_, opts_, work_);

    ocp_nlp_dynamics_config *config = config_;
    ocp_nlp_dynamics_cont_dims *dims = dims_;
    ocp_nlp_dynamics_cont_opts *opts = opts_;
    ocp_nlp_dynamics_cont_workspace *work = work_;
    ocp_nlp_dynamics_cont_memory *mem = mem_;
    ocp_nlp_dynamics_cont_model *model = model_;

    int nx = dims->nx;
    int nu = dims->nu;
    int nz = dims->nz;
    int nx1 = dims->nx1;
    int nu1 = dims->nu1;

    // setup model
    work->sim_in->model = model->sim_model;
    work->sim_in->T = model->T;

    blasfeo_unpack_dvec(nz, mem->z, 0, work->sim_in->z);

    // pass state and control to integrator
    blasfeo_unpack_dvec(nu, mem->ux, 0, work->sim_in->u);
    blasfeo_unpack_dvec(nx, mem->ux, nu, work->sim_in->x);

    // initialize seeds
    for (int jj = 0; jj < nx1 * (nx + nu); jj++) work->sim_in->S_forw[jj] = 0.0;
    for (int jj = 0; jj < nx1; jj++) work->sim_in->S_forw[jj * (nx + 1)] = 1.0;
    for (int jj = 0; jj < nx + nu; jj++) work->sim_in->S_adj[jj] = 0.0;
    blasfeo_unpack_dvec(nx1, mem->pi, 0, work->sim_in->S_adj);

    // call integrator
    config->sim_solver->evaluate(config->sim_solver, work->sim_in, work->sim_out, opts->sim_solver,
                                 mem->sim_solver, work->sim_solver);

    // TODO(rien): transition functions for changing dimensions not yet implemented!

    // B
    blasfeo_pack_tran_dmat(nx1, nu, work->sim_out->S_forw + nx1 * nx, nx1, mem->BAbt, 0, 0);
    // A
    blasfeo_pack_tran_dmat(nx1, nx, work->sim_out->S_forw + 0, nx1, mem->BAbt, nu, 0);

    // fun
    blasfeo_pack_dvec(nx1, work->sim_out->xn, &mem->fun, 0);
    blasfeo_daxpy(nx1, -1.0, mem->ux1, nu1, &mem->fun, 0, &mem->fun, 0);

    // adj TODO if not computed by the integrator
    if (opts->compute_adj)
    {
        blasfeo_dgemv_n(nu+nx, nx1, -1.0, mem->BAbt, 0, 0, mem->pi, 0, 0.0, &mem->adj, 0, &mem->adj,
                        0);
        blasfeo_dveccp(nx1, mem->pi, 0, &mem->adj, nu + nx);
    }

    if (opts->compute_hess)
    {
        // unpack d*_d2u
        blasfeo_pack_dmat(nu, nu, &work->sim_out->S_hess[(nx+nu)*nx + nx], nx + nu,
                                     &mem->hes, 0, 0);
        // unpack d*_dux: mem-hess: nx x nu
        blasfeo_pack_dmat(nx, nu, &work->sim_out->S_hess[(nx + nu)*nx], nx + nu, &mem->hes, nu, 0);
        // unpack d*_d2x
        blasfeo_pack_dmat(nx, nx, &work->sim_out->S_hess[0], nx + nu, &mem->hes, nu, nu);

        // Add hessian contribution
        blasfeo_dgead(nx+nu, nx+nu, 1.0, &mem->hes, 0, 0, mem->RSQrq, 0, 0);
    }
    return;
}

int ocp_nlp_dynamics_cont_precompute(void *config_, void *dims_, void *model_, void *opts_,
                                        void *mem_, void *work_)
{
    ocp_nlp_dynamics_cont_cast_workspace(config_, dims_, opts_, work_);

    ocp_nlp_dynamics_config *config = config_;
    // ocp_nlp_dynamics_cont_dims *dims = dims_;
    ocp_nlp_dynamics_cont_opts *opts = opts_;
    ocp_nlp_dynamics_cont_workspace *work = work_;
    ocp_nlp_dynamics_cont_memory *mem = mem_;
    ocp_nlp_dynamics_cont_model *model = model_;
    work->sim_in->model = model->sim_model;
    work->sim_in->T = model->T;

    // call integrator
    int status = config->sim_solver->precompute(config->sim_solver, work->sim_in, work->sim_out,
                                   opts->sim_solver, mem->sim_solver, work->sim_solver);
    return status;
}

void ocp_nlp_dynamics_cont_config_initialize_default(void *config_)
{
    ocp_nlp_dynamics_config *config = config_;

    config->dims_calculate_size = &ocp_nlp_dynamics_cont_dims_calculate_size;
    config->dims_assign = &ocp_nlp_dynamics_cont_dims_assign;
    config->dims_set = &ocp_nlp_dynamics_cont_dims_set;
    config->dims_initialize = &ocp_nlp_dynamics_cont_dims_initialize;
    config->model_calculate_size = &ocp_nlp_dynamics_cont_model_calculate_size;
    config->model_assign = &ocp_nlp_dynamics_cont_model_assign;
    config->model_set_T = &ocp_nlp_dynamics_cont_model_set_T;
    config->opts_calculate_size = &ocp_nlp_dynamics_cont_opts_calculate_size;
    config->opts_assign = &ocp_nlp_dynamics_cont_opts_assign;
    config->opts_initialize_default = &ocp_nlp_dynamics_cont_opts_initialize_default;
    config->opts_update = &ocp_nlp_dynamics_cont_opts_update;
    config->opts_set = &ocp_nlp_dynamics_cont_opts_set;
    config->memory_calculate_size = &ocp_nlp_dynamics_cont_memory_calculate_size;
    config->memory_assign = &ocp_nlp_dynamics_cont_memory_assign;
    config->memory_get_fun_ptr = &ocp_nlp_dynamics_cont_memory_get_fun_ptr;
    config->memory_get_adj_ptr = &ocp_nlp_dynamics_cont_memory_get_adj_ptr;
    config->memory_set_ux_ptr = &ocp_nlp_dynamics_cont_memory_set_ux_ptr;
    config->memory_set_ux1_ptr = &ocp_nlp_dynamics_cont_memory_set_ux1_ptr;
    config->memory_set_pi_ptr = &ocp_nlp_dynamics_cont_memory_set_pi_ptr;
    config->memory_set_BAbt_ptr = &ocp_nlp_dynamics_cont_memory_set_BAbt_ptr;
    config->memory_set_RSQrq_ptr = &ocp_nlp_dynamics_cont_memory_set_RSQrq_ptr;
    config->memory_set_z_ptr = &ocp_nlp_dynamics_cont_memory_set_z_ptr;
    config->workspace_calculate_size = &ocp_nlp_dynamics_cont_workspace_calculate_size;
    config->initialize = &ocp_nlp_dynamics_cont_initialize;
    config->update_qp_matrices = &ocp_nlp_dynamics_cont_update_qp_matrices;
    config->precompute = &ocp_nlp_dynamics_cont_precompute;
    config->config_initialize_default = &ocp_nlp_dynamics_cont_config_initialize_default;

    return;
}
