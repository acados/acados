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


#include "acados/sim/sim_lifted_irk_integrator.h"

// standard
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// acados
#include "acados/utils/mem.h"
#include "acados/utils/print.h"

#include "acados/sim/sim_common.h"

#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_blas.h"
#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_v_aux_ext_dep.h"


/************************************************
* dims
************************************************/

acados_size_t sim_lifted_irk_dims_calculate_size()
{
    acados_size_t size = sizeof(sim_lifted_irk_dims);

    return size;
}

void *sim_lifted_irk_dims_assign(void *config_, void *raw_memory)
{
    char *c_ptr = raw_memory;

    sim_lifted_irk_dims *dims = (sim_lifted_irk_dims *) c_ptr;
    c_ptr += sizeof(sim_lifted_irk_dims);

    dims->nx = 0;
    dims->nu = 0;
    dims->nz = 0;

    assert((char *) raw_memory + sim_lifted_irk_dims_calculate_size() >= c_ptr);

    return dims;
}



void sim_lifted_irk_dims_set(void *config_, void *dims_, const char *field, const int *value)
{
    sim_lifted_irk_dims *dims = dims_;

    if (!strcmp(field, "nx"))
    {
        dims->nx = *value;
    }
    else if (!strcmp(field, "nu"))
    {
        dims->nu = *value;
    }
    else if (!strcmp(field, "nz"))
    {
        dims->nz = *value;
    }
    else
    {
        printf("\nerror: sim_lifted_irk_dims_set: field not available: %s\n", field);
        exit(1);
    }
}



void sim_lifted_irk_dims_get(void *config_, void *dims_, const char *field, int *value)
{
    sim_lifted_irk_dims *dims = dims_;

    if (!strcmp(field, "nx"))
    {
        *value = dims->nx;
    }
    else if (!strcmp(field, "nu"))
    {
        *value = dims->nu;
    }
    else if (!strcmp(field, "nz"))
    {
        *value = dims->nz;
    }
    else
    {
        printf("\nerror: sim_lifted_irk_dims_get: field not available: %s\n", field);
        exit(1);
    }
}



/************************************************
* model
************************************************/

acados_size_t sim_lifted_irk_model_calculate_size(void *config, void *dims)
{
    acados_size_t size = 0;

    size += sizeof(lifted_irk_model);

    return size;
}



void *sim_lifted_irk_model_assign(void *config, void *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    lifted_irk_model *data = (lifted_irk_model *) c_ptr;
    c_ptr += sizeof(lifted_irk_model);

    assert((char *) raw_memory + sim_lifted_irk_model_calculate_size(config, dims) >= c_ptr);

    return data;
}



int sim_lifted_irk_model_set(void *model_, const char *field, void *value)
{
    lifted_irk_model *model = model_;

    if (!strcmp(field, "impl_ode_fun") || !strcmp(field, "impl_dae_fun"))
    {
        model->impl_ode_fun = value;
    }
    else if (!strcmp(field, "impl_ode_fun_jac_x_xdot_u") || !strcmp(field, "impl_dae_fun_jac_x_xdot_u"))
    {
        model->impl_ode_fun_jac_x_xdot_u = value;
    }
    else
    {
        printf("\nerror: sim_lifted_irk_model_set: wrong field: %s\n", field);
        exit(1);
    }

    return ACADOS_SUCCESS;
}



/************************************************
* opts
************************************************/

acados_size_t sim_lifted_irk_opts_calculate_size(void *config_, void *dims)
{
    int ns_max = NS_MAX;

    acados_size_t size = 0;

    size += sizeof(sim_opts);

    size += ns_max * ns_max * sizeof(double);  // A_mat
    size += ns_max * sizeof(double);           // b_vec
    size += ns_max * sizeof(double);           // c_vec

    size += butcher_tableau_work_calculate_size(ns_max);

    make_int_multiple_of(8, &size);
    size += 1 * 8;

    return size;
}



void *sim_lifted_irk_opts_assign(void *config_, void *dims, void *raw_memory)
{
    int ns_max = NS_MAX;

    char *c_ptr = (char *) raw_memory;

    sim_opts *opts = (sim_opts *) c_ptr;
    c_ptr += sizeof(sim_opts);

    align_char_to(8, &c_ptr);

    // work
    opts->work = c_ptr;
    c_ptr += butcher_tableau_work_calculate_size(ns_max);

    assign_and_advance_double(ns_max * ns_max, &opts->A_mat, &c_ptr);
    assign_and_advance_double(ns_max, &opts->b_vec, &c_ptr);
    assign_and_advance_double(ns_max, &opts->c_vec, &c_ptr);

    assert((char *) raw_memory + sim_lifted_irk_opts_calculate_size(config_, dims) >= c_ptr);

    return (void *) opts;
}



void sim_lifted_irk_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    sim_opts *opts = opts_;
    sim_lifted_irk_dims *dims = (sim_lifted_irk_dims *) dims_;

    // default options
    opts->newton_iter = 3;
    // opts->scheme = NULL;
    opts->num_steps = 2;
    opts->num_forw_sens = dims->nx + dims->nu;
    opts->sens_forw = true;
    opts->sens_adj = false;
    opts->sens_hess = false;
    opts->jac_reuse = true;
    opts->exact_z_output = false;
    opts->ns = 3;
    opts->collocation_type = GAUSS_LEGENDRE;

    assert(opts->ns <= NS_MAX && "ns > NS_MAX!");

    // butcher tableau
    calculate_butcher_tableau(opts->ns, opts->collocation_type, opts->c_vec, opts->b_vec, opts->A_mat, opts->work);
    // for consistency check
    opts->tableau_size = opts->ns;

    // TODO(oj): check if constr h or cost depend on z, turn on in this case only.
    if (dims->nz > 0)
    {
        opts->output_z = true;
        opts->sens_algebraic = true;
    }
    else
    {
        opts->output_z = false;
        opts->sens_algebraic = false;
    }

    return;
}



void sim_lifted_irk_opts_update(void *config_, void *dims, void *opts_)
{
    sim_opts *opts = opts_;

    assert(opts->ns <= NS_MAX && "ns > NS_MAX!");

    calculate_butcher_tableau(opts->ns, opts->collocation_type, opts->c_vec, opts->b_vec, opts->A_mat, opts->work);

    opts->tableau_size = opts->ns;

    return;
}



void sim_lifted_irk_opts_set(void *config_, void *opts_, const char *field, void *value)
{
    sim_opts *opts = (sim_opts *) opts_;
    sim_opts_set_(opts, field, value);
}



void sim_lifted_irk_opts_get(void *config_, void *opts_, const char *field, void *value)
{
    sim_opts *opts = (sim_opts *) opts_;
    sim_opts_get_(config_, opts, field, value);
}



/************************************************
* memory
************************************************/

acados_size_t sim_lifted_irk_memory_calculate_size(void *config, void *dims_, void *opts_)
{
    sim_opts *opts = opts_;
    sim_lifted_irk_dims *dims = (sim_lifted_irk_dims *) dims_;

    int ns = opts->ns;

    int nx = dims->nx;
    int nu = dims->nu;

    int num_steps = opts->num_steps;

    acados_size_t size = sizeof(sim_lifted_irk_memory);

    size += 1 * sizeof(struct blasfeo_dmat);            // S_forw
    size += 2 * sizeof(struct blasfeo_dmat);            // JGK, JGf
    size += (num_steps) * sizeof(struct blasfeo_dmat);  // JKf
    size += (num_steps) * sizeof(struct blasfeo_dvec);  // K
    size += 2 * sizeof(struct blasfeo_dvec);            // x, u

    size += blasfeo_memsize_dmat(nx, nx + nu);                    // S_forw
    size += blasfeo_memsize_dmat(nx * ns, nx * ns);               // JGK
    size += 1 * blasfeo_memsize_dmat(nx * ns, nx + nu);           // JGf
    size += (num_steps) *blasfeo_memsize_dmat(nx * ns, nx + nu);  // JKf
    size += (num_steps) *blasfeo_memsize_dvec(nx * ns);           // K
    size += 1 * blasfeo_memsize_dvec(nx);                         // x
    size += 1 * blasfeo_memsize_dvec(nu);                         // u

    size += 1 * 8; // initial align
    make_int_multiple_of(64, &size);
    size += 1 * 64;

    return size;
}



void *sim_lifted_irk_memory_assign(void *config, void *dims_, void *opts_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    sim_opts *opts = opts_;
    sim_lifted_irk_dims *dims = (sim_lifted_irk_dims *) dims_;

    int ns = opts->ns;

    int nx = dims->nx;
    int nu = dims->nu;

    int num_steps = opts->num_steps;

    // initial align
    align_char_to(8, &c_ptr);

    sim_lifted_irk_memory *memory = (sim_lifted_irk_memory *) c_ptr;
    c_ptr += sizeof(sim_lifted_irk_memory);

    memory->S_forw = (struct blasfeo_dmat *) c_ptr;
    c_ptr += sizeof(struct blasfeo_dmat);

    memory->JGK = (struct blasfeo_dmat *) c_ptr;
    c_ptr += sizeof(struct blasfeo_dmat);

    memory->JGf = (struct blasfeo_dmat *) c_ptr;
    c_ptr += sizeof(struct blasfeo_dmat);

    assign_and_advance_blasfeo_dmat_structs(num_steps, &memory->JKf, &c_ptr);

    assign_and_advance_blasfeo_dvec_structs(num_steps, &memory->K, &c_ptr);

    memory->x = (struct blasfeo_dvec *) c_ptr;
    c_ptr += sizeof(struct blasfeo_dvec);

    memory->u = (struct blasfeo_dvec *) c_ptr;
    c_ptr += sizeof(struct blasfeo_dvec);

    align_char_to(64, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nx, nx + nu, memory->S_forw, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx * ns, nx * ns, memory->JGK, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx * ns, nx + nu, memory->JGf, &c_ptr);
    for (int i = 0; i < num_steps; i++)
    {
        assign_and_advance_blasfeo_dmat_mem(nx * ns, nx + nu, &memory->JKf[i], &c_ptr);
        blasfeo_dgese(nx * ns, nx + nu, 0.0, &memory->JKf[i], 0, 0);
    }

    for (int i = 0; i < num_steps; i++)
    {
        assign_and_advance_blasfeo_dvec_mem(nx * ns, &memory->K[i], &c_ptr);
        blasfeo_dvecse(nx * ns, 0.0, &memory->K[i], 0);
    }

    assign_and_advance_blasfeo_dvec_mem(nx, memory->x, &c_ptr);
    blasfeo_dvecse(nx, 0.0, memory->x, 0);
    assign_and_advance_blasfeo_dvec_mem(nu, memory->u, &c_ptr);
    blasfeo_dvecse(nu, 0.0, memory->u, 0);

    // memory->init_K = 0;

    // TODO(andrea): need to move this to options.
    memory->update_sens = 1;

    assert((char *) raw_memory + sim_lifted_irk_memory_calculate_size(config, dims, opts_) >=
           c_ptr);

    return (void *) memory;
}



int sim_lifted_irk_memory_set(void *config_, void *dims_, void *mem_, const char *field, void *value)
{
    // sim_config *config = config_;
    // sim_lifted_irk_memory *mem = (sim_lifted_irk_memory *) mem_;

    printf("sim_lifted_irk_memory_set field %s is not supported! \n", field);
    exit(1);

    // sim_config *config = config_;
    // sim_lifted_irk_memory *mem = (sim_lifted_irk_memory *) mem_;
    // if (!strcmp(field, "xdot_guess"))
    // {
    //     int nx;
    //     config->dims_get(config_, dims_, "nx", &nx);
    //     double *xdot = value;
    //     blasfeo_pack_dvec(nx, xdot, 0, &mem->K[0], 0);
    //     mem->init_K = 1;
    // }
    // else if (!strcmp(field, "guesses_blasfeo"))
    // {
    //     int nx, nz;
    //     config->dims_get(config_, dims_, "nx", &nx);
    //     config->dims_get(config_, dims_, "nz", &nz);

    //     struct blasfeo_dvec *sim_guess = (struct blasfeo_dvec *) value;
    //     blasfeo_dveccp(nx+nz, sim_guess, 0, &mem->K[0], 0);
    //     mem->init_K = 1;
    // }
    // else
    // {
    //     printf("sim_lifted_irk_memory_set field %s is not supported! \n", field);
    //     exit(1);
    // }
    // return ACADOS_SUCCESS;
}



int sim_lifted_irk_memory_set_to_zero(void *config_, void * dims_, void *opts_, void *mem_, const char *field)
{
    sim_config *config = config_;
    sim_lifted_irk_memory *mem = (sim_lifted_irk_memory *) mem_;
    sim_opts *opts = (sim_opts *) opts_;

    int status = ACADOS_SUCCESS;

    if (!strcmp(field, "guesses"))
    {
        int nx;
        config->dims_get(config_, dims_, "nx", &nx);
        for (int i = 0; i < opts->num_steps; i++)
        {
            blasfeo_dvecse(nx * opts->ns, 0.0, &mem->K[i], 0);
        }
    }
    else
    {
        printf("sim_lifted_irk_memory_set_to_zero field %s is not supported! \n", field);
        exit(1);
    }

    return status;
}



void sim_lifted_irk_memory_get(void *config_, void *dims_, void *mem_, const char *field, void *value)
{
    sim_lifted_irk_memory *mem = mem_;

    if (!strcmp(field, "time_sim"))
    {
		double *ptr = value;
		*ptr = mem->time_sim;
	}
    else if (!strcmp(field, "time_sim_ad"))
    {
		double *ptr = value;
		*ptr = mem->time_ad;
	}
    else if (!strcmp(field, "time_sim_la"))
    {
		double *ptr = value;
		*ptr = mem->time_la;
	}
	else
	{
		printf("sim_lifted_irk_memory_get field %s is not supported! \n", field);
		exit(1);
	}
}



/************************************************
* workspace
************************************************/

acados_size_t sim_lifted_irk_workspace_calculate_size(void *config_, void *dims_, void *opts_)
{
    sim_opts *opts = opts_;
    sim_lifted_irk_dims *dims = (sim_lifted_irk_dims *) dims_;

    int ns = opts->ns;

    int nx = dims->nx;
    int nu = dims->nu;

    acados_size_t size = sizeof(sim_lifted_irk_workspace);

    size += 3 * sizeof(struct blasfeo_dmat);  // J_temp_x, J_temp_xdot, J_temp_u

    size += 5 * sizeof(struct blasfeo_dvec);  // rG, xt, xn, xn_out, dxn
    size += 1 * sizeof(struct blasfeo_dvec);  // w ([x; u])

    size += 2 * blasfeo_memsize_dmat(nx, nx);  // J_temp_x, J_temp_xdot
    size += blasfeo_memsize_dmat(nx, nu);      // J_temp_u

    size += 1 * blasfeo_memsize_dvec(nx * ns);  // rG
    size += 4 * blasfeo_memsize_dvec(nx);       // xt, xn, xn_out, dxn
    size += blasfeo_memsize_dvec(nx + nu);      // w

    size += nx * ns * sizeof(int);  // ipiv

    make_int_multiple_of(64, &size);
    size += 1 * 64;

    return size;
}



static void *sim_lifted_irk_cast_workspace(void *config_, void *dims_, void *opts_,
                                               void *raw_memory)
{
    sim_opts *opts = opts_;
    sim_lifted_irk_dims *dims = (sim_lifted_irk_dims *) dims_;

    int ns = opts->ns;

    int nx = dims->nx;
    int nu = dims->nu;

    char *c_ptr = (char *) raw_memory;

    sim_lifted_irk_workspace *workspace = (sim_lifted_irk_workspace *) c_ptr;
    c_ptr += sizeof(sim_lifted_irk_workspace);

    workspace->J_temp_x = (struct blasfeo_dmat *) c_ptr;
    c_ptr += sizeof(struct blasfeo_dmat);

    workspace->J_temp_xdot = (struct blasfeo_dmat *) c_ptr;
    c_ptr += sizeof(struct blasfeo_dmat);
    workspace->J_temp_u = (struct blasfeo_dmat *) c_ptr;
    c_ptr += sizeof(struct blasfeo_dmat);

    workspace->rG = (struct blasfeo_dvec *) c_ptr;
    c_ptr += sizeof(struct blasfeo_dvec);

    workspace->xt = (struct blasfeo_dvec *) c_ptr;
    c_ptr += sizeof(struct blasfeo_dvec);

    workspace->xn = (struct blasfeo_dvec *) c_ptr;
    c_ptr += sizeof(struct blasfeo_dvec);

    workspace->xn_out = (struct blasfeo_dvec *) c_ptr;
    c_ptr += sizeof(struct blasfeo_dvec);

    workspace->dxn = (struct blasfeo_dvec *) c_ptr;
    c_ptr += sizeof(struct blasfeo_dvec);

    workspace->w = (struct blasfeo_dvec *) c_ptr;
    c_ptr += sizeof(struct blasfeo_dvec);

    align_char_to(64, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nx, nx, workspace->J_temp_x, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx, nx, workspace->J_temp_xdot, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx, nu, workspace->J_temp_u, &c_ptr);

    assign_and_advance_blasfeo_dvec_mem(nx * ns, workspace->rG, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nx, workspace->xt, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nx, workspace->xn, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nx, workspace->xn_out, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nx, workspace->dxn, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nx + nu, workspace->w, &c_ptr);

    assign_and_advance_int(nx * ns, &workspace->ipiv, &c_ptr);

    assert((char *) raw_memory +
               sim_lifted_irk_workspace_calculate_size(config_, dims, opts_) >=
           c_ptr);

    return (void *) workspace;
}



int sim_lifted_irk_precompute(void *config_, sim_in *in, sim_out *out, void *opts_, void *mem_,
                       void *work_)
{
    return ACADOS_SUCCESS;
}



/************************************************
* functions
************************************************/

int sim_lifted_irk(void *config_, sim_in *in, sim_out *out, void *opts_, void *mem_,
                       void *work_)
{
    // typecasting
    sim_config *config = config_;
    sim_opts *opts = opts_;
    sim_lifted_irk_memory *mem = mem_;

    void *dims_ = in->dims;
    sim_lifted_irk_dims *dims = (sim_lifted_irk_dims *) dims_;

    sim_lifted_irk_workspace *workspace =
        (sim_lifted_irk_workspace *) sim_lifted_irk_cast_workspace(config, dims, opts,
                                                                           work_);

    int nx = dims->nx;
    int nu = dims->nu;
    int nz = dims->nz;

    int ns = opts->ns;

    if ( opts->ns != opts->tableau_size )
    {
        printf("Error in sim_lifted_irk: the Butcher tableau size does not match ns");
        exit(1);
    }
    // assert - only use supported features
    if (nz != 0)
    {
        printf("nz should be zero - DAEs are not supported by the lifted IRK integrator");
        exit(1);
    }
    if (opts->output_z)
    {
        printf("opts->output_z should be false - DAEs are not supported for the lifted IRK integrator");
        exit(1);
    }
    if (opts->sens_algebraic)
    {
        printf("opts->sens_algebraic should be false - DAEs are not supported for the lifted IRK integrator");
        exit(1);
    }

    int ii, jj, ss;
    double a;


    double *x = in->x;
    double *u = in->u;
    double *S_forw_in = in->S_forw;

    // int newton_iter = opts->newton_iter; // not used; always 1 in lifted

    double *A_mat = opts->A_mat;
    double *b_vec = opts->b_vec;
    int num_steps = opts->num_steps;

    double step = in->T / num_steps;
    // TODO(FreyJo): this should be an option!
    int update_sens = mem->update_sens;

    int *ipiv = workspace->ipiv;
    struct blasfeo_dmat *JGK = mem->JGK;
    struct blasfeo_dmat *S_forw = mem->S_forw;

    struct blasfeo_dmat *J_temp_x = workspace->J_temp_x;
    struct blasfeo_dmat *J_temp_xdot = workspace->J_temp_xdot;
    struct blasfeo_dmat *J_temp_u = workspace->J_temp_u;

    struct blasfeo_dvec *rG = workspace->rG;
    struct blasfeo_dvec *K = mem->K;
    struct blasfeo_dmat *JGf = mem->JGf;
    struct blasfeo_dmat *JKf = mem->JKf;
    struct blasfeo_dvec *xt = workspace->xt;
    struct blasfeo_dvec *xn = workspace->xn;
    struct blasfeo_dvec *xn_out = workspace->xn_out;
    struct blasfeo_dvec *dxn = workspace->dxn;

    struct blasfeo_dvec *w = workspace->w;

    double *x_out = out->xn;
    double *S_forw_out = out->S_forw;

    struct blasfeo_dvec_args ext_fun_in_K;

    ext_fun_arg_t ext_fun_type_in[3];
    void *ext_fun_in[3];

    struct blasfeo_dvec_args ext_fun_out_rG;
    ext_fun_arg_t ext_fun_type_out[5];
    void *ext_fun_out[5];

    lifted_irk_model *model = in->model;

    acados_timer timer, timer_ad, timer_la;
    double timing_ad = 0.0;
    out->info->LAtime = 0.0;

    if (opts->sens_hess)
    {
        printf("LIFTED_IRK with HESSIAN PROPAGATION - NOT IMPLEMENTED YET - EXITING.");
        exit(1);
    }
    if (opts->sens_adj)
    {
        printf("LIFTED_IRK with ADJOINT SENSITIVITIES - NOT IMPLEMENTED YET - EXITING.");
        exit(1);
    }

    // if (mem->init_K)
    // {
    //     for (ss = 0; ss < num_steps; ss++)
    //     {
    //         for (ii = 0; ii < ns; ii++)
    //         {
    //             if (!(ii == 0 && ss == 0))
    //                 blasfeo_dveccp(nx, &mem->K[0], 0, &mem->K[ss], ii*nx);
    //         }
    //     }
    //     mem->init_K = 0;
    // }

    blasfeo_dgese(nx, nx, 0.0, J_temp_x, 0, 0);
    blasfeo_dgese(nx, nx, 0.0, J_temp_xdot, 0, 0);
    blasfeo_dgese(nx, nu, 0.0, J_temp_x, 0, 0);

    blasfeo_dvecse(nx * ns, 0.0, rG, 0);

    // TODO(dimitris): shouldn't this be NF instead of nx+nu??
    if (update_sens) blasfeo_pack_dmat(nx, nx + nu, S_forw_in, nx, S_forw, 0, 0);

    blasfeo_dvecse(nx * ns, 0.0, rG, 0);
    blasfeo_pack_dvec(nx, x, 1, xn, 0);
    blasfeo_pack_dvec(nx, x, 1, xn_out, 0);
    blasfeo_dvecse(nx, 0.0, dxn, 0);


    // start the loop
    acados_tic(&timer);
    for (ss = 0; ss < num_steps; ss++)
    {
        // initialize
        blasfeo_dgese(nx * ns, nx * ns, 0.0, JGK, 0, 0);
        blasfeo_dgese(nx * ns, nx + nu, 0.0, JGf, 0, 0);

        // expansion step (K variables)
        // compute x and u step
        blasfeo_pack_dvec(nx, in->x, 1, w, 0);
        blasfeo_pack_dvec(nu, in->u, 1, w, nx);

        blasfeo_daxpy(nx, -1.0, mem->x, 0, w, 0, w, 0);
        blasfeo_daxpy(nu, -1.0, mem->u, 0, w, nx, w, nx);
        blasfeo_dgemv_n(nx * ns, nx + nu, 1.0, &JKf[ss], 0, 0, w, 0, 1.0, &K[ss], 0, &K[ss], 0);

        blasfeo_pack_dvec(nx, in->x, 1, mem->x, 0);
        blasfeo_pack_dvec(nu, in->u, 1, mem->u, 0);

        // reset value of JKf
        blasfeo_dgese(nx * ns, nx + nu, 0.0, &JKf[ss], 0, 0);

        for (ii = 0; ii < ns; ii++)  // ii-th row of tableau
        {
            // take x(n); copy a strvec into a strvec
            blasfeo_dveccp(nx, xn, 0, xt, 0);

            for (jj = 0; jj < ns; jj++)
            {  // jj-th col of tableau
                a = A_mat[ii + ns * jj];
                if (a != 0)
                {
                    // xt = xt + T_int * a[i,j]*K_j
                    a *= step;
                    blasfeo_daxpy(nx, a, &K[ss], jj * nx, xt, 0, xt, 0);
                }
            }

            if (!update_sens)
            {
                // compute the residual of implicit ode at time t_ii, store value in rGt
                acados_tic(&timer_ad);

                ext_fun_type_in[0] = BLASFEO_DVEC;
                ext_fun_in[0] = xt;  // x: nx
                ext_fun_type_in[1] = BLASFEO_DVEC_ARGS;
                ext_fun_in_K.xi = ii * nx;
                ext_fun_in_K.x = &K[ss];
                ext_fun_in[1] = &ext_fun_in_K;  // K[ii*nx]: nx
                ext_fun_type_in[2] = COLMAJ;
                ext_fun_in[2] = u;  // u: nu

                ext_fun_type_out[0] = BLASFEO_DVEC_ARGS;
                ext_fun_out_rG.x = rG;
                ext_fun_out_rG.xi = ii * nx;
                ext_fun_out[0] = &ext_fun_out_rG;  // fun: nx

                model->impl_ode_fun->evaluate(model->impl_ode_fun, ext_fun_type_in, ext_fun_in,
                                              ext_fun_type_out, ext_fun_out);

                timing_ad += acados_toc(&timer_ad);
            }
            else
            {
                // compute the jacobian of implicit ode
                acados_tic(&timer_ad);

                ext_fun_type_in[0] = BLASFEO_DVEC;
                ext_fun_in[0] = xt;  // x: nx
                ext_fun_type_in[1] = BLASFEO_DVEC_ARGS;
                ext_fun_in_K.xi = ii * nx;  // K[ii*nx]: nx
                ext_fun_in_K.x = &K[ss];
                ext_fun_in[1] = &ext_fun_in_K;
                ext_fun_type_in[2] = COLMAJ;
                ext_fun_in[2] = u;  // u: nu

                ext_fun_type_out[0] = BLASFEO_DVEC_ARGS;
                ext_fun_out_rG.x = rG;
                ext_fun_out_rG.xi = ii * nx;
                ext_fun_out[0] = &ext_fun_out_rG;  // fun: nx
                ext_fun_type_out[1] = BLASFEO_DMAT;
                ext_fun_out[1] = J_temp_x;  // jac_x: nx*nx
                ext_fun_type_out[2] = BLASFEO_DMAT;
                ext_fun_out[2] = J_temp_xdot;  // jac_xdot: nx*nx
                ext_fun_type_out[3] = BLASFEO_DMAT;
                ext_fun_out[3] = J_temp_u;  // jac_u: nx*nu

                model->impl_ode_fun_jac_x_xdot_u->evaluate(model->impl_ode_fun_jac_x_xdot_u,
                                                           ext_fun_type_in, ext_fun_in,
                                                           ext_fun_type_out, ext_fun_out);

                timing_ad += acados_toc(&timer_ad);

                blasfeo_dgecp(nx, nx, J_temp_x, 0, 0, JGf, ii * nx, 0);
                blasfeo_dgecp(nx, nu, J_temp_u, 0, 0, JGf, ii * nx, nx);

                for (jj = 0; jj < ns; jj++)
                {
                    // compute the block (ii,jj)th block of JGK
                    a = A_mat[ii + ns * jj];
                    if (a != 0)
                    {
                        a *= step;
                        blasfeo_dgead(nx, nx, a, J_temp_x, 0, 0, JGK, ii * nx, jj * nx);
                    }
                    if (jj == ii)
                    {
                        blasfeo_dgead(nx, nx, 1, J_temp_xdot, 0, 0, JGK, ii * nx, jj * nx);
                    }
                }  // end jj
            }
        }  // end ii

        // obtain x(n+1) before updating K(n)
        for (ii = 0; ii < ns; ii++)
            blasfeo_daxpy(nx, step * b_vec[ii], &K[ss], ii * nx, xn, 0, xn, 0);

        // DGETRF computes an LU factorization of a general M-by-N matrix A
        // using partial pivoting with row interchanges.

        if (update_sens)
        {
            blasfeo_dgetrf_rp(nx * ns, nx * ns, JGK, 0, 0, JGK, 0, 0, ipiv);
        }

        // update r.h.s (6.23, Quirynen2017)
        blasfeo_dgemv_n(nx * ns, nx, 1.0, JGf, 0, 0, dxn, 0, 1.0, rG, 0, rG, 0);


        // permute also the r.h.s
        blasfeo_dvecpe(nx * ns, ipiv, rG, 0);

        // solve JGK * y = rG, JGK on the (l)eft, (l)ower-trian, (n)o-trans
        //                    (u)nit trian
        blasfeo_dtrsv_lnu(nx * ns, JGK, 0, 0, rG, 0, rG, 0);

        // solve JGK * x = rG, JGK on the (l)eft, (u)pper-trian, (n)o-trans
        //                    (n)o unit trian , and store x in rG
        blasfeo_dtrsv_unn(nx * ns, JGK, 0, 0, rG, 0, rG, 0);


        // scale and add a generic strmat into a generic strmat // K = K - rG, where rG is DeltaK
        blasfeo_daxpy(nx * ns, -1.0, rG, 0, &K[ss], 0, &K[ss], 0);

        // obtain dx(n)
        for (ii = 0; ii < ns; ii++)
            blasfeo_daxpy(nx, -step * b_vec[ii], rG, ii * nx, dxn, 0, dxn, 0);

        // update JKf
        // JKf[ss] = JGf * S_forw;
        if (in->identity_seed && ss == 0) // omit matrix multiplication for identity seed
            blasfeo_dgecp(nx * ns, nx + nu, JGf, 0, 0, &JKf[ss], 0, 0);
        else
        {
            blasfeo_dgemm_nn(nx * ns, nx + nu, nx, 1.0, JGf, 0, 0, S_forw, 0, 0, 0.0, &JKf[ss], 0, 0,
                            &JKf[ss], 0, 0);
            blasfeo_dgead(nx * ns, nu, 1.0, JGf, 0, nx, &JKf[ss], 0, nx);
        }

        // solve linear system
        acados_tic(&timer_la);
        blasfeo_drowpe(nx * ns, ipiv, &JKf[ss]);
        blasfeo_dtrsm_llnu(nx * ns, nx + nu, 1.0, JGK, 0, 0, &JKf[ss], 0, 0, &JKf[ss], 0, 0);
        blasfeo_dtrsm_lunn(nx * ns, nx + nu, 1.0, JGK, 0, 0, &JKf[ss], 0, 0, &JKf[ss], 0, 0);
        out->info->LAtime += acados_toc(&timer_la);

        // update forward sensitivity
        for (jj = 0; jj < ns; jj++)
            blasfeo_dgead(nx, nx + nu, -step * b_vec[jj], &JKf[ss], jj * nx, 0, S_forw, 0, 0);

        // obtain x(n+1)
        for (ii = 0; ii < ns; ii++)
            blasfeo_daxpy(nx, step * b_vec[ii], &K[ss], ii * nx, xn_out, 0, xn_out, 0);

    }  // end int step ss


    // extract output
    blasfeo_unpack_dvec(nx, xn_out, 0, x_out, 1);

    blasfeo_unpack_dmat(nx, nx + nu, S_forw, 0, 0, S_forw_out, nx);

    out->info->CPUtime = acados_toc(&timer);
    out->info->ADtime = timing_ad;

	mem->time_sim = out->info->CPUtime;
	mem->time_ad = out->info->ADtime;
	mem->time_la = out->info->LAtime;

    return 0;
}



void sim_lifted_irk_config_initialize_default(void *config_)
{
    sim_config *config = config_;

    config->evaluate = &sim_lifted_irk;
    config->precompute = &sim_lifted_irk_precompute;
    config->opts_calculate_size = &sim_lifted_irk_opts_calculate_size;
    config->opts_assign = &sim_lifted_irk_opts_assign;
    config->opts_initialize_default = &sim_lifted_irk_opts_initialize_default;
    config->opts_update = &sim_lifted_irk_opts_update;
    config->opts_set = &sim_lifted_irk_opts_set;
    config->opts_get = &sim_lifted_irk_opts_get;
    config->memory_calculate_size = &sim_lifted_irk_memory_calculate_size;
    config->memory_assign = &sim_lifted_irk_memory_assign;
    config->memory_set = &sim_lifted_irk_memory_set;
    config->memory_set_to_zero = &sim_lifted_irk_memory_set_to_zero;
    config->memory_get = &sim_lifted_irk_memory_get;
    config->workspace_calculate_size = &sim_lifted_irk_workspace_calculate_size;
    config->model_calculate_size = &sim_lifted_irk_model_calculate_size;
    config->model_assign = &sim_lifted_irk_model_assign;
    config->model_set = &sim_lifted_irk_model_set;
    config->dims_calculate_size = &sim_lifted_irk_dims_calculate_size;
    config->dims_assign = &sim_lifted_irk_dims_assign;
    config->dims_set = &sim_lifted_irk_dims_set;
    config->dims_get = &sim_lifted_irk_dims_get;
    return;
}
