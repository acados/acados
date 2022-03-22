/*
 * Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
 * Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
 * Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
 * Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
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

acados_size_t ocp_nlp_dynamics_cont_dims_calculate_size(void *config_)
{
    ocp_nlp_dynamics_config *dyn_config = (ocp_nlp_dynamics_config *) config_;
    sim_config *sim_sol_config = (sim_config *) dyn_config->sim_solver;
    acados_size_t size = 0;

    size += sizeof(ocp_nlp_dynamics_cont_dims);

    size += sim_sol_config->dims_calculate_size();

    size += 2*8; // align

    return size;
}


void *ocp_nlp_dynamics_cont_dims_assign(void *config_, void *raw_memory)
{
    ocp_nlp_dynamics_config *dyn_config = (ocp_nlp_dynamics_config *) config_;
    sim_config *sim_sol_config = (sim_config *) dyn_config->sim_solver;

    char *c_ptr = (char *) raw_memory;

    align_char_to(8, &c_ptr);
    ocp_nlp_dynamics_cont_dims *dims = (ocp_nlp_dynamics_cont_dims *) c_ptr;
    c_ptr += sizeof(ocp_nlp_dynamics_cont_dims);

    align_char_to(8, &c_ptr);
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

    // NOTE: this could be done just in precompute.
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



void ocp_nlp_dynamics_cont_dims_get(void *config_, void *dims_, const char *field, int* value)
{
    ocp_nlp_dynamics_cont_dims *dims = dims_;
    if (!strcmp(field, "nx"))
    {
        *value = dims->nx;
    }
    else if (!strcmp(field, "nu"))
    {
        *value = dims->nu;
    }
    else if (!strcmp(field, "nx1"))
    {
        *value = dims->nx1;
    }
    else
    {
        // get GNSF dims from integrator module
        ocp_nlp_dynamics_config *dyn_config = (ocp_nlp_dynamics_config *) config_;
        ocp_nlp_dynamics_cont_dims *dims = (ocp_nlp_dynamics_cont_dims *) dims_;
        sim_config *sim_config_ = (sim_config *) dyn_config->sim_solver;

        sim_config_->dims_get(sim_config_, dims->sim, field, value);
    }
// printf("\nexiting: ocp_nlp_dynamics_cont_dims_get, %d\n", *value);

}



/************************************************
 * options
 ************************************************/

acados_size_t ocp_nlp_dynamics_cont_opts_calculate_size(void *config_, void *dims_)
{
    ocp_nlp_dynamics_config *config = config_;
    ocp_nlp_dynamics_cont_dims *dims = dims_;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_dynamics_cont_opts);

    size += config->sim_solver->opts_calculate_size(config->sim_solver, dims->sim);
    size += 8;

    return size;
}



void *ocp_nlp_dynamics_cont_opts_assign(void *config_, void *dims_, void *raw_memory)
{
    ocp_nlp_dynamics_config *config = config_;
    ocp_nlp_dynamics_cont_dims *dims = dims_;

    char *c_ptr = (char *) raw_memory;

    ocp_nlp_dynamics_cont_opts *opts = (ocp_nlp_dynamics_cont_opts *) c_ptr;
    c_ptr += sizeof(ocp_nlp_dynamics_cont_opts);

    align_char_to(8, &c_ptr);
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

    // own opts
    opts->compute_adj = 1;
    opts->compute_hess = 0;

    // sim opts
    config->sim_solver->opts_initialize_default(config->sim_solver, dims->sim, opts->sim_solver);

	// overwrite defaults

    bool sens_forw = true;
    bool sens_adj = false;
    bool sens_hess = false;

    config->sim_solver->opts_set(config->sim_solver, opts->sim_solver, "sens_forw", &sens_forw);
    config->sim_solver->opts_set(config->sim_solver, opts->sim_solver, "sens_adj", &sens_adj);
    config->sim_solver->opts_set(config->sim_solver, opts->sim_solver, "sens_hess", &sens_hess);

    return;
}



void ocp_nlp_dynamics_cont_opts_update(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dynamics_config *config = config_;
    ocp_nlp_dynamics_cont_dims *dims = dims_;
    ocp_nlp_dynamics_cont_opts *opts = opts_;

//    if (!opts->compute_hess)
//    {
//        bool bool_false = false;
//        config->sim_solver->opts_set( config->sim_solver, opts->sim_solver, "sens_hess", &bool_false );
//    }
    config->sim_solver->opts_update(config->sim_solver, dims->sim, opts->sim_solver);

    return;
}




void ocp_nlp_dynamics_cont_opts_set(void *config_, void *opts_, const char *field, void* value)
{
    ocp_nlp_dynamics_config *config = config_;
    ocp_nlp_dynamics_cont_opts *opts = opts_;
    sim_config *sim_config_ = config->sim_solver;

    if (!strcmp(field, "compute_adj"))
    {
        int *int_ptr = value;
        opts->compute_adj = *int_ptr;
		// TODO set in the sim solver too ???
    }
    else if (!strcmp(field, "compute_hess"))
    {
        int *int_ptr = value;
        opts->compute_hess = *int_ptr;
        bool tmp_bool = true;
        if (*int_ptr==0)
        {
            tmp_bool = false;
        }
        config->sim_solver->opts_set(config->sim_solver, opts->sim_solver, "sens_adj", &tmp_bool);
        config->sim_solver->opts_set(config->sim_solver, opts->sim_solver, "sens_hess", &tmp_bool);
    }
    else
    {
        sim_config_->opts_set(sim_config_, opts->sim_solver, field, value);
    }

    return;

}



/************************************************
 * memory
 ************************************************/

acados_size_t ocp_nlp_dynamics_cont_memory_calculate_size(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dynamics_config *config = config_;
    ocp_nlp_dynamics_cont_dims *dims = dims_;
    ocp_nlp_dynamics_cont_opts *opts = opts_;

    // extract dims
    int nx = dims->nx;
    // int nz = dims->nz;
    int nu = dims->nu;
    int nx1 = dims->nx1;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_dynamics_cont_memory);

    size += 1 * blasfeo_memsize_dvec(nu + nx + nx1);  // adj
    size += 1 * blasfeo_memsize_dvec(nx1);            // fun

    size +=
        config->sim_solver->memory_calculate_size(config->sim_solver, dims->sim, opts->sim_solver);

    size += 1*64;  // blasfeo_mem align

    make_int_multiple_of(8, &size);

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
    // int nz = dims->nz;
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



void ocp_nlp_dynamics_cont_memory_set_tmp_ux_ptr(struct blasfeo_dvec *tmp_ux, void *memory_)
{
    ocp_nlp_dynamics_cont_memory *memory = memory_;

    memory->tmp_ux = tmp_ux;

    return;
}



void ocp_nlp_dynamics_cont_memory_set_ux1_ptr(struct blasfeo_dvec *ux1, void *memory_)
{
    ocp_nlp_dynamics_cont_memory *memory = memory_;

    memory->ux1 = ux1;

    return;
}



void ocp_nlp_dynamics_cont_memory_set_tmp_ux1_ptr(struct blasfeo_dvec *tmp_ux1, void *memory_)
{
    ocp_nlp_dynamics_cont_memory *memory = memory_;

    memory->tmp_ux1 = tmp_ux1;

    return;
}



void ocp_nlp_dynamics_cont_memory_set_pi_ptr(struct blasfeo_dvec *pi, void *memory_)
{
    ocp_nlp_dynamics_cont_memory *memory = memory_;

    memory->pi = pi;

    return;
}



void ocp_nlp_dynamics_cont_memory_set_tmp_pi_ptr(struct blasfeo_dvec *tmp_pi, void *memory_)
{
    ocp_nlp_dynamics_cont_memory *memory = memory_;

    memory->tmp_pi = tmp_pi;

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



void ocp_nlp_dynamics_cont_memory_set_dzduxt_ptr(struct blasfeo_dmat *mat, void *memory_)
{
    ocp_nlp_dynamics_cont_memory *memory = memory_;

    memory->dzduxt = mat;

    return;
}



void ocp_nlp_dynamics_cont_memory_set_sim_guess_ptr(struct blasfeo_dvec *vec, bool *bool_ptr, void *memory_)
{
    ocp_nlp_dynamics_cont_memory *memory = memory_;

    memory->sim_guess = vec;
    memory->set_sim_guess = bool_ptr;

    return;
}



void ocp_nlp_dynamics_cont_memory_set_z_alg_ptr(struct blasfeo_dvec *vec, void *memory_)
{
    ocp_nlp_dynamics_cont_memory *memory = memory_;

    memory->z_alg = vec;

    return;
}



void ocp_nlp_dynamics_cont_memory_get(void *config_, void *dims_, void *mem_, const char *field, void* value)
{
    ocp_nlp_dynamics_config *config = config_;
    ocp_nlp_dynamics_cont_dims *dims = dims_;
    ocp_nlp_dynamics_cont_memory *mem = mem_;

	sim_config *sim = config->sim_solver;

    if (!strcmp(field, "time_sim") || !strcmp(field, "time_sim_ad") || !strcmp(field, "time_sim_la"))
    {
		sim->memory_get(sim, dims->sim, mem->sim_solver, field, value);
    }
    else
    {
		printf("\nerror: ocp_nlp_dynamics_cont_memory_get: field %s not available\n", field);
		exit(1);
    }

}



/************************************************
 * workspace
 ************************************************/

acados_size_t ocp_nlp_dynamics_cont_workspace_calculate_size(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dynamics_config *config = config_;
    ocp_nlp_dynamics_cont_dims *dims = dims_;
    ocp_nlp_dynamics_cont_opts *opts = opts_;

    // extract dims
    int nx = dims->nx;
    int nu = dims->nu;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_dynamics_cont_workspace);

    size += sim_in_calculate_size(config->sim_solver, dims->sim);
    size += sim_out_calculate_size(config->sim_solver, dims->sim);
    size += config->sim_solver->workspace_calculate_size(config->sim_solver, dims->sim, opts->sim_solver);

    size += 1 * blasfeo_memsize_dmat(nu+nx, nu+nx);   // hess

    size += 1*64;  // blasfeo_mem align
    make_int_multiple_of(8, &size);

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

    // extract dims
    int nx = dims->nx;
    int nu = dims->nu;

    // sim in
    work->sim_in = sim_in_assign(config->sim_solver, dims->sim, c_ptr);
    c_ptr += sim_in_calculate_size(config->sim_solver, dims->sim);
    // sim out
    work->sim_out = sim_out_assign(config->sim_solver, dims->sim, c_ptr);
    c_ptr += sim_out_calculate_size(config->sim_solver, dims->sim);
    // workspace
    work->sim_solver = c_ptr;
    c_ptr += config->sim_solver->workspace_calculate_size(config->sim_solver, dims->sim, opts->sim_solver);

    // blasfeo_mem align
    align_char_to(64, &c_ptr);

    // hess
    assign_and_advance_blasfeo_dmat_mem(nu+nx, nu+nx, &work->hess, &c_ptr);

    assert((char *) work + ocp_nlp_dynamics_cont_workspace_calculate_size(config, dims, opts) >= c_ptr);

    return;
}



/************************************************
 * model
 ************************************************/

acados_size_t ocp_nlp_dynamics_cont_model_calculate_size(void *config_, void *dims_)
{
    ocp_nlp_dynamics_config *config = config_;
    ocp_nlp_dynamics_cont_dims *dims = dims_;

    // extract dims
    // int nx = dims->nx;
    // int nu = dims->nu;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_dynamics_cont_model);

    size += config->sim_solver->model_calculate_size(config->sim_solver, dims->sim);
    size += 1*8;
    make_int_multiple_of(8, &size);

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
    align_char_to(8, &c_ptr);

    model->sim_model = config->sim_solver->model_assign(config->sim_solver, dims->sim, c_ptr);
    c_ptr += config->sim_solver->model_calculate_size(config->sim_solver, dims->sim);

    assert((char *) raw_memory + ocp_nlp_dynamics_cont_model_calculate_size(config, dims) >= c_ptr);

    return model;
}



void ocp_nlp_dynamics_cont_model_set(void *config_, void *dims_, void *model_, const char *field, void *value)
{
    ocp_nlp_dynamics_config *config = config_;
    ocp_nlp_dynamics_cont_model *model = model_;

    sim_config *sim_config = config->sim_solver;

    if (!strcmp(field, "T"))
    {
        double *T = (double *) value;
        model->T = *T;
    }
    else
    {
        int status = sim_config->model_set(model->sim_model, field, value);
        if (status!=0)
        {
            printf("\nerror: field %s not available in module ocp_nlp_dynamics_cont_model_set\n", field);
            exit(1);
        }
    }

    return;
}



/************************************************
 * functions
 ************************************************/

void ocp_nlp_dynamics_cont_initialize(void *config_, void *dims_, void *model_, void *opts_, void *mem_, void *work_)
{
    return;
}



void ocp_nlp_dynamics_cont_update_qp_matrices(void *config_, void *dims_, void *model_, void *opts_, void *mem_, void *work_)
{
    ocp_nlp_dynamics_cont_cast_workspace(config_, dims_, opts_, work_);

    ocp_nlp_dynamics_config *config = config_;
    ocp_nlp_dynamics_cont_dims *dims = dims_;
    ocp_nlp_dynamics_cont_opts *opts = opts_;
    ocp_nlp_dynamics_cont_workspace *work = work_;
    ocp_nlp_dynamics_cont_memory *mem = mem_;
    ocp_nlp_dynamics_cont_model *model = model_;

    int jj;

    int nx = dims->nx;
    int nu = dims->nu;
    int nz = dims->nz;
    int nx1 = dims->nx1;
    int nu1 = dims->nu1;

    // setup model
    work->sim_in->model = model->sim_model;
    work->sim_in->T = model->T;

    // pass state and control to integrator
    blasfeo_unpack_dvec(nu, mem->ux, 0, work->sim_in->u, 1);
    blasfeo_unpack_dvec(nx, mem->ux, nu, work->sim_in->x, 1);

    // printf("sim_guess, bool %d\n", mem->set_sim_guess[0]);
    // blasfeo_print_exp_dvec(nx + nz, mem->sim_guess, 0);

    if (mem->set_sim_guess!=NULL && mem->set_sim_guess[0])
    {
        config->sim_solver->memory_set(config->sim_solver, work->sim_in->dims, mem->sim_solver,
                                        "guesses_blasfeo", mem->sim_guess);
        // only use/pass the initial guess once
        mem->set_sim_guess[0] = false;
    }

    // initialize seeds
    // TODO fix dims if nx!=nx1 !!!!!!!!!!!!!!!!!
    // set S_forw = [eye(nx), zeros(nx x nu)]
    for(jj = 0; jj < nx1 * (nx + nu); jj++)
        work->sim_in->S_forw[jj] = 0.0;
    for(jj = 0; jj < nx1; jj++)
        work->sim_in->S_forw[jj * (nx + 1)] = 1.0;
    work->sim_in->identity_seed = true;

    // adjoint seed
    for(jj = 0; jj < nx + nu; jj++)
        work->sim_in->S_adj[jj] = 0.0;
    blasfeo_unpack_dvec(nx1, mem->pi, 0, work->sim_in->S_adj, 1);

    // call integrator
    config->sim_solver->evaluate(config->sim_solver, work->sim_in, work->sim_out, opts->sim_solver,
            mem->sim_solver, work->sim_solver);

    // TODO transition functions for changing dimensions not yet implemented!

    // B
    blasfeo_pack_tran_dmat(nx1, nu, work->sim_out->S_forw + nx1 * nx, nx1, mem->BAbt, 0, 0);
    // A
    blasfeo_pack_tran_dmat(nx1, nx, work->sim_out->S_forw + 0, nx1, mem->BAbt, nu, 0);
    // dzduxt
    blasfeo_pack_tran_dmat(nz, nu, work->sim_out->S_algebraic + nx*nz, nz, mem->dzduxt, 0, 0);
    blasfeo_pack_tran_dmat(nz, nx, work->sim_out->S_algebraic + 0, nz, mem->dzduxt, nu, 0);
    // blasfeo_print_dmat(nx + nu, nz, mem->dzduxt, 0, 0);

    // function
    blasfeo_pack_dvec(nx1, work->sim_out->xn, 1, &mem->fun, 0);
    blasfeo_daxpy(nx1, -1.0, mem->ux1, nu1, &mem->fun, 0, &mem->fun, 0);
    blasfeo_pack_dvec(nz, work->sim_out->zn, 1, mem->z_alg, 0);

    // adjoint
    if (opts->compute_adj)
    {
        // check if adjoints computed in integrator
        bool adjoint_integrator;
        sim_opts_get(config->sim_solver, opts->sim_solver, "sens_adj", &adjoint_integrator);

        if (adjoint_integrator)
        {
            blasfeo_pack_dvec(nu, work->sim_out->S_adj+nx, 1, &mem->adj, 0);
            blasfeo_pack_dvec(nx, work->sim_out->S_adj+0, 1, &mem->adj, nu);
            blasfeo_dvecsc(nu+nx, -1.0, &mem->adj, 0);
        }
        // compute as forward * adj_seed
        else
        {
            blasfeo_dgemv_n(nu+nx, nx1, -1.0, mem->BAbt, 0, 0, mem->pi, 0, 0.0, &mem->adj, 0, &mem->adj, 0);
        }
        blasfeo_dveccp(nx1, mem->pi, 0, &mem->adj, nu+nx);
    }

    // hessian
    if (opts->compute_hess)
    {

//        d_print_mat(nu+nx, nu+nx, work->sim_out->S_hess, nu+nx);

        // unpack d*_d2u
        blasfeo_pack_dmat(nu, nu, &work->sim_out->S_hess[(nx+nu)*nx + nx], nx+nu, &work->hess, 0, 0);
        // unpack d*_dux: mem-hess: nx x nu
        blasfeo_pack_dmat(nx, nu, &work->sim_out->S_hess[(nx + nu)*nx], nx+nu, &work->hess, nu, 0);
        // unpack d*_d2x
        blasfeo_pack_dmat(nx, nx, &work->sim_out->S_hess[0], nx+nu, &work->hess, nu, nu);

        // Add hessian contribution
        blasfeo_dgead(nx+nu, nx+nu, 1.0, &work->hess, 0, 0, mem->RSQrq, 0, 0);
    }

    return;

}



void ocp_nlp_dynamics_cont_compute_fun(void *config_, void *dims_, void *model_, void *opts_, void *mem_, void *work_)
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
    // int nz = dims->nz;
    int nx1 = dims->nx1;
    int nu1 = dims->nu1;

    // setup model
    work->sim_in->model = model->sim_model;
    work->sim_in->T = model->T;

    // pass state and control to integrator
    blasfeo_unpack_dvec(nu, mem->tmp_ux, 0, work->sim_in->u, 1);
    blasfeo_unpack_dvec(nx, mem->tmp_ux, nu, work->sim_in->x, 1);

    // printf("sim_guess, bool %d\n", mem->set_sim_guess[0]);
    // blasfeo_print_exp_dvec(nx + nz, mem->sim_guess, 0);

    if (mem->set_sim_guess!=NULL && mem->set_sim_guess[0])
    {
        config->sim_solver->memory_set(config->sim_solver, work->sim_in->dims, mem->sim_solver,
                                       "guesses_blasfeo", mem->sim_guess);
        // only use/pass the initial guess once
        mem->set_sim_guess[0] = false;
    }

	// backup sens options
	bool sens_forw_bkp, sens_adj_bkp, sens_hess_bkp;
    config->sim_solver->opts_get(config->sim_solver, opts->sim_solver, "sens_forw", &sens_forw_bkp);
    config->sim_solver->opts_get(config->sim_solver, opts->sim_solver, "sens_adj", &sens_adj_bkp);
    config->sim_solver->opts_get(config->sim_solver, opts->sim_solver, "sens_hess", &sens_hess_bkp);

	// set all sens to false
	bool sens_all = false;
    config->sim_solver->opts_set(config->sim_solver, opts->sim_solver, "sens_forw", &sens_all);
    config->sim_solver->opts_set(config->sim_solver, opts->sim_solver, "sens_adj", &sens_all);
    config->sim_solver->opts_set(config->sim_solver, opts->sim_solver, "sens_hess", &sens_all);

    // call integrator
    config->sim_solver->evaluate(config->sim_solver, work->sim_in, work->sim_out, opts->sim_solver,
            mem->sim_solver, work->sim_solver);

	// restore sens options
    config->sim_solver->opts_set(config->sim_solver, opts->sim_solver, "sens_forw", &sens_forw_bkp);
    config->sim_solver->opts_set(config->sim_solver, opts->sim_solver, "sens_adj", &sens_adj_bkp);
    config->sim_solver->opts_set(config->sim_solver, opts->sim_solver, "sens_hess", &sens_hess_bkp);

    // function
    blasfeo_pack_dvec(nx1, work->sim_out->xn, 1, &mem->fun, 0);
    // fun -= x[next_stage]
    blasfeo_daxpy(nx1, -1.0, mem->tmp_ux1, nu1, &mem->fun, 0, &mem->fun, 0);
//    blasfeo_pack_dvec(nz, work->sim_out->zn, 1, mem->z_alg, 0);
    // printf("\ndyn_cont: compute f:\n");
    // blasfeo_print_exp_tran_dvec(nx1, &mem->fun, 0);

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
    
    config->sim_solver->memory_set_to_zero(config->sim_solver, work->sim_in->dims,
                                    opts->sim_solver, mem->sim_solver, "guesses");

    return status;
}



void ocp_nlp_dynamics_cont_config_initialize_default(void *config_)
{
    ocp_nlp_dynamics_config *config = config_;

    config->dims_calculate_size = &ocp_nlp_dynamics_cont_dims_calculate_size;
    config->dims_assign = &ocp_nlp_dynamics_cont_dims_assign;
    config->dims_set = &ocp_nlp_dynamics_cont_dims_set;
    config->dims_get = &ocp_nlp_dynamics_cont_dims_get;
    config->dims_initialize = &ocp_nlp_dynamics_cont_dims_initialize;
    config->model_calculate_size = &ocp_nlp_dynamics_cont_model_calculate_size;
    config->model_assign = &ocp_nlp_dynamics_cont_model_assign;
    config->model_set = &ocp_nlp_dynamics_cont_model_set;
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
    config->memory_set_tmp_ux_ptr = &ocp_nlp_dynamics_cont_memory_set_tmp_ux_ptr;
    config->memory_set_ux1_ptr = &ocp_nlp_dynamics_cont_memory_set_ux1_ptr;
    config->memory_set_tmp_ux1_ptr = &ocp_nlp_dynamics_cont_memory_set_tmp_ux1_ptr;
    config->memory_set_pi_ptr = &ocp_nlp_dynamics_cont_memory_set_pi_ptr;
    config->memory_set_tmp_pi_ptr = &ocp_nlp_dynamics_cont_memory_set_tmp_pi_ptr;
    config->memory_set_BAbt_ptr = &ocp_nlp_dynamics_cont_memory_set_BAbt_ptr;
    config->memory_set_RSQrq_ptr = &ocp_nlp_dynamics_cont_memory_set_RSQrq_ptr;
    config->memory_set_dzduxt_ptr = &ocp_nlp_dynamics_cont_memory_set_dzduxt_ptr;
    config->memory_set_sim_guess_ptr = &ocp_nlp_dynamics_cont_memory_set_sim_guess_ptr;
    config->memory_set_z_alg_ptr = &ocp_nlp_dynamics_cont_memory_set_z_alg_ptr;
    config->memory_get = &ocp_nlp_dynamics_cont_memory_get;
    config->workspace_calculate_size = &ocp_nlp_dynamics_cont_workspace_calculate_size;
    config->initialize = &ocp_nlp_dynamics_cont_initialize;
    config->update_qp_matrices = &ocp_nlp_dynamics_cont_update_qp_matrices;
    config->compute_fun = &ocp_nlp_dynamics_cont_compute_fun;
    config->precompute = &ocp_nlp_dynamics_cont_precompute;
    config->config_initialize_default = &ocp_nlp_dynamics_cont_config_initialize_default;

    return;
}
