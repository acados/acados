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

#include "acados/ocp_nlp/ocp_nlp_dynamics_disc.h"

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

int ocp_nlp_dynamics_disc_dims_calculate_size(void *config_)
{
    int size = 0;

    size += sizeof(ocp_nlp_dynamics_disc_dims);

    return size;
}

void *ocp_nlp_dynamics_disc_dims_assign(void *config_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_nlp_dynamics_disc_dims *dims = (ocp_nlp_dynamics_disc_dims *) c_ptr;
    c_ptr += sizeof(ocp_nlp_dynamics_disc_dims);

    assert((char *) raw_memory + ocp_nlp_dynamics_disc_dims_calculate_size(config_) >= c_ptr);

    return dims;
}

void ocp_nlp_dynamics_disc_dims_initialize(void *config_, void *dims_, int nx, int nu, int nx1,
                                           int nu1)
{
    ocp_nlp_dynamics_disc_dims *dims = dims_;

    dims->nx = nx;
    dims->nu = nu;
    dims->nx1 = nx1;
    dims->nu1 = nu1;

    return;
}

/************************************************
 * options
 ************************************************/

int ocp_nlp_dynamics_disc_opts_calculate_size(void *config_, void *dims_)
{
    // ocp_nlp_dynamics_config *config = config_;
    // ocp_nlp_dynamics_disc_dims *config = dims_;

    int size = 0;

    size += sizeof(ocp_nlp_dynamics_disc_opts);

    return size;
}

void *ocp_nlp_dynamics_disc_opts_assign(void *config_, void *dims_, void *raw_memory)
{
    // ocp_nlp_dynamics_config *config = config_;
    // ocp_nlp_dynamics_disc_dims *dims = dims_;

    char *c_ptr = (char *) raw_memory;

    ocp_nlp_dynamics_disc_opts *opts = (ocp_nlp_dynamics_disc_opts *) c_ptr;
    c_ptr += sizeof(ocp_nlp_dynamics_disc_opts);

    assert((char *) raw_memory + ocp_nlp_dynamics_disc_opts_calculate_size(config_, dims_) >=
           c_ptr);

    return opts;
}

void ocp_nlp_dynamics_disc_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    // ocp_nlp_dynamics_config *config = config_;
    // ocp_nlp_dynamics_disc_opts *opts = opts_;

    return;
}

void ocp_nlp_dynamics_disc_opts_update(void *config_, void *dims_, void *opts_)
{
    // ocp_nlp_dynamics_config *config = config_;
    // ocp_nlp_dynamics_disc_opts *opts = opts_;

    return;
}

/************************************************
 * memory
 ************************************************/

int ocp_nlp_dynamics_disc_memory_calculate_size(void *config_, void *dims_, void *opts_)
{
    // ocp_nlp_dynamics_config *config = config_;
    ocp_nlp_dynamics_disc_dims *dims = dims_;
    // ocp_nlp_dynamics_disc_opts *opts = opts_;

    // extract dims
    int nx = dims->nx;
    int nu = dims->nu;
    int nx1 = dims->nx1;

    int size = 0;

    size += sizeof(ocp_nlp_dynamics_disc_memory);

    size += 1 * blasfeo_memsize_dvec(nu + nx + nx1);  // adj
    size += 1 * blasfeo_memsize_dvec(nx1);            // fun

    size += 64;  // blasfeo_mem align

    return size;
}

void *ocp_nlp_dynamics_disc_memory_assign(void *config_, void *dims_, void *opts_, void *raw_memory)
{
    // ocp_nlp_dynamics_config *config = config_;
    ocp_nlp_dynamics_disc_dims *dims = dims_;
    // ocp_nlp_dynamics_disc_opts *opts = opts_;

    char *c_ptr = (char *) raw_memory;

    // extract dims
    int nx = dims->nx;
    int nu = dims->nu;
    int nx1 = dims->nx1;

    // struct
    ocp_nlp_dynamics_disc_memory *memory = (ocp_nlp_dynamics_disc_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_dynamics_disc_memory);

    // blasfeo_mem align
    align_char_to(64, &c_ptr);

    // adj
    assign_and_advance_blasfeo_dvec_mem(nu + nx + nx1, &memory->adj, &c_ptr);
    // fun
    assign_and_advance_blasfeo_dvec_mem(nx1, &memory->fun, &c_ptr);

    assert((char *) raw_memory +
               ocp_nlp_dynamics_disc_memory_calculate_size(config_, dims, opts_) >=
           c_ptr);

    return memory;
}

struct blasfeo_dvec *ocp_nlp_dynamics_disc_memory_get_fun_ptr(void *memory_)
{
    ocp_nlp_dynamics_disc_memory *memory = memory_;

    return &memory->fun;
}

struct blasfeo_dvec *ocp_nlp_dynamics_disc_memory_get_adj_ptr(void *memory_)
{
    ocp_nlp_dynamics_disc_memory *memory = memory_;

    return &memory->adj;
}

void ocp_nlp_dynamics_disc_memory_set_ux_ptr(struct blasfeo_dvec *ux, void *memory_)
{
    ocp_nlp_dynamics_disc_memory *memory = memory_;

    memory->ux = ux;

    return;
}

void ocp_nlp_dynamics_disc_memory_set_ux1_ptr(struct blasfeo_dvec *ux1, void *memory_)
{
    ocp_nlp_dynamics_disc_memory *memory = memory_;

    memory->ux1 = ux1;

    return;
}

void ocp_nlp_dynamics_disc_memory_set_pi_ptr(struct blasfeo_dvec *pi, void *memory_)
{
    ocp_nlp_dynamics_disc_memory *memory = memory_;

    memory->pi = pi;

    return;
}

void ocp_nlp_dynamics_disc_memory_set_BAbt_ptr(struct blasfeo_dmat *BAbt, void *memory_)
{
    ocp_nlp_dynamics_disc_memory *memory = memory_;

    memory->BAbt = BAbt;

    return;
}

/************************************************
 * workspace
 ************************************************/

int ocp_nlp_dynamics_disc_workspace_calculate_size(void *config_, void *dims_, void *opts_)
{
    // ocp_nlp_dynamics_config *config = config_;
    ocp_nlp_dynamics_disc_dims *dims = dims_;
    // ocp_nlp_dynamics_disc_opts *opts = opts_;

    int nx = dims->nx;
    int nu = dims->nu;
    int nx1 = dims->nx1;

    int size = 0;

    size += sizeof(ocp_nlp_dynamics_disc_workspace);

    size += (nx + nu) * sizeof(double);                // discrete_model_in
    size += (nx1 + nx1 * (nx + nu)) * sizeof(double);  // discrete_model_out

    size += 8;

    return size;
}

static void ocp_nlp_dynamics_disc_cast_workspace(void *config_, void *dims_, void *opts_,
                                                 void *work_)
{
    // ocp_nlp_dynamics_config *config = config_;
    ocp_nlp_dynamics_disc_dims *dims = dims_;
    // ocp_nlp_dynamics_disc_opts *opts = opts_;
    ocp_nlp_dynamics_disc_workspace *work = work_;

    int nx = dims->nx;
    int nu = dims->nu;
    int nx1 = dims->nx1;

    char *c_ptr = (char *) work_;
    c_ptr += sizeof(ocp_nlp_dynamics_disc_workspace);

    align_char_to(8, &c_ptr);

    // discrete_model_in
    assign_and_advance_double(nx + nu, &work->discrete_model_in, &c_ptr);
    // discrete_model_out
    assign_and_advance_double(nx1 + nx1 * (nx + nu), &work->discrete_model_out, &c_ptr);

    assert((char *) work + ocp_nlp_dynamics_disc_workspace_calculate_size(config_, dims, opts_) >=
           c_ptr);

    return;
}

/************************************************
 * model
 ************************************************/

int ocp_nlp_dynamics_disc_model_calculate_size(void *config_, void *dims_)
{
    // ocp_nlp_dynamics_config *config = config_;

    // extract dims
    // int nx = dims->nx;
    // int nu = dims->nu;

    int size = 0;

    size += sizeof(ocp_nlp_dynamics_disc_model);

    return size;
}

void *ocp_nlp_dynamics_disc_model_assign(void *config_, void *dims_, void *raw_memory)
{
    // ocp_nlp_dynamics_config *config = config_;

    char *c_ptr = (char *) raw_memory;

    // extract dims
    // int nx = dims->nx;
    // int nu = dims->nu;

    // struct
    ocp_nlp_dynamics_disc_model *model = (ocp_nlp_dynamics_disc_model *) c_ptr;
    c_ptr += sizeof(ocp_nlp_dynamics_disc_model);

    assert((char *) raw_memory + ocp_nlp_dynamics_disc_model_calculate_size(config_, dims_) >=
           c_ptr);

    return model;
}

void ocp_nlp_dynamics_disc_model_set_T(double T, void *model_) { return; }
/************************************************
 * functions
 ************************************************/

void ocp_nlp_dynamics_disc_initialize(void *config_, void *dims_, void *model_, void *opts_,
                                      void *mem_, void *work_)
{
    return;
}

void ocp_nlp_dynamics_disc_update_qp_matrices(void *config_, void *dims_, void *model_, void *opts_,
                                              void *mem_, void *work_)
{
    ocp_nlp_dynamics_disc_cast_workspace(config_, dims_, opts_, work_);

    // ocp_nlp_dynamics_config *config = config_;
    ocp_nlp_dynamics_disc_dims *dims = dims_;
    // ocp_nlp_dynamics_disc_opts *opts = opts_;
    ocp_nlp_dynamics_disc_workspace *work = work_;
    ocp_nlp_dynamics_disc_memory *mem = mem_;
    ocp_nlp_dynamics_disc_model *model = model_;

    int nx = dims->nx;
    int nu = dims->nu;
    int nx1 = dims->nx1;
    int nu1 = dims->nu1;

    ext_fun_arg_t ext_fun_type_in[2];
    void *ext_fun_in[2];  // XXX large enough ?
    ext_fun_arg_t ext_fun_type_out[2];
    void *ext_fun_out[2];  // XXX large enough ?

    // pass state and control to integrator
    blasfeo_unpack_dvec(nu, mem->ux, 0, work->discrete_model_in + nx);
    blasfeo_unpack_dvec(nx, mem->ux, nu, work->discrete_model_in);

    ext_fun_type_in[0] = COLMAJ;
    ext_fun_in[0] = work->discrete_model_in + 0;  // x: nx
    ext_fun_type_in[1] = COLMAJ;
    ext_fun_in[1] = work->discrete_model_in + nx;  // u: nu

    ext_fun_type_out[0] = COLMAJ;
    ext_fun_out[0] = work->discrete_model_out + 0;  // fun: nx1
    ext_fun_type_out[1] = COLMAJ;
    ext_fun_out[1] = work->discrete_model_out + nx1;  // jac: nx1*(nx+nu)

    // call external function
    model->discrete_model->evaluate(model->discrete_model, ext_fun_type_in, ext_fun_in,
                                    ext_fun_type_out, ext_fun_out);

    // B
    blasfeo_pack_tran_dmat(nx1, nu, work->discrete_model_out + nx1 + nx1 * nx, nx1, mem->BAbt, 0,
                           0);
    // A
    blasfeo_pack_tran_dmat(nx1, nx, work->discrete_model_out + nx1, nx1, mem->BAbt, nu, 0);

    // fun
    blasfeo_pack_dvec(nx1, work->discrete_model_out, &mem->fun, 0);
    blasfeo_daxpy(nx1, -1.0, mem->ux1, nu1, &mem->fun, 0, &mem->fun, 0);

    // adj TODO if not computed by the external function
    blasfeo_dgemv_n(nu + nx, nx1, -1.0, mem->BAbt, 0, 0, mem->pi, 0, 0.0, &mem->adj, 0, &mem->adj,
                    0);
    blasfeo_dveccp(nx1, mem->pi, 0, &mem->adj, nu + nx);

    return;
}

void ocp_nlp_dynamics_disc_config_initialize_default(void *config_)
{
    ocp_nlp_dynamics_config *config = config_;

    config->dims_calculate_size = &ocp_nlp_dynamics_disc_dims_calculate_size;
    config->dims_assign = &ocp_nlp_dynamics_disc_dims_assign;
    config->dims_initialize = &ocp_nlp_dynamics_disc_dims_initialize;
    config->model_calculate_size = &ocp_nlp_dynamics_disc_model_calculate_size;
    config->model_assign = &ocp_nlp_dynamics_disc_model_assign;
    config->model_set_T = &ocp_nlp_dynamics_disc_model_set_T;
    config->opts_calculate_size = &ocp_nlp_dynamics_disc_opts_calculate_size;
    config->opts_assign = &ocp_nlp_dynamics_disc_opts_assign;
    config->opts_initialize_default = &ocp_nlp_dynamics_disc_opts_initialize_default;
    config->opts_update = &ocp_nlp_dynamics_disc_opts_update;
    config->memory_calculate_size = &ocp_nlp_dynamics_disc_memory_calculate_size;
    config->memory_assign = &ocp_nlp_dynamics_disc_memory_assign;
    config->memory_get_fun_ptr = &ocp_nlp_dynamics_disc_memory_get_fun_ptr;
    config->memory_get_adj_ptr = &ocp_nlp_dynamics_disc_memory_get_adj_ptr;
    config->memory_set_ux_ptr = &ocp_nlp_dynamics_disc_memory_set_ux_ptr;
    config->memory_set_ux1_ptr = &ocp_nlp_dynamics_disc_memory_set_ux1_ptr;
    config->memory_set_pi_ptr = &ocp_nlp_dynamics_disc_memory_set_pi_ptr;
    config->memory_set_BAbt_ptr = &ocp_nlp_dynamics_disc_memory_set_BAbt_ptr;
    config->workspace_calculate_size = &ocp_nlp_dynamics_disc_workspace_calculate_size;
    config->initialize = &ocp_nlp_dynamics_disc_initialize;
    config->update_qp_matrices = &ocp_nlp_dynamics_disc_update_qp_matrices;
    config->config_initialize_default = &ocp_nlp_dynamics_disc_config_initialize_default;

    return;
}
