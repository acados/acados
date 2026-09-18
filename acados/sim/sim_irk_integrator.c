/*
 * Copyright (c) The acados authors.
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


#include "acados/sim/sim_irk_integrator.h"

// standard
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// acados
#include "acados/utils/mem.h"
#include "acados/utils/print.h"
#include "acados/utils/math.h"
#include "acados/ocp_nlp/ocp_nlp_cost_common.h"

#include "acados/sim/sim_common.h"

#include "blasfeo_d_aux.h"
#include "blasfeo_d_blas.h"
#include "blasfeo_common.h"


/************************************************
 * dims
 ************************************************/

acados_size_t sim_irk_dims_calculate_size()
{
    acados_size_t size = sizeof(sim_irk_dims);

    return size;
}

void *sim_irk_dims_assign(void *config_, void *raw_memory)
{
    char *c_ptr = raw_memory;

    sim_irk_dims *dims = (sim_irk_dims *) c_ptr;
    c_ptr += sizeof(sim_irk_dims);

    dims->nx = 0;
    dims->nu = 0;
    dims->nz = 0;
    dims->np = 0;

    assert((char *) raw_memory + sim_irk_dims_calculate_size() >= c_ptr);

    return dims;
}



void sim_irk_dims_set(void *config_, void *dims_, const char *field, const int *value)
{
    sim_irk_dims *dims = dims_;

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
    else if (!strcmp(field, "np"))
    {
        dims->np = *value;
    }
    else if (!strcmp(field, "np_global"))
    {
        // np_global dimension not needed
    }
    else
    {
        printf("\nerror: sim_irk_dims_set: field not available: %s\n", field);
        exit(1);
    }
}



void sim_irk_dims_get(void *config_, void *dims_, const char *field, int *value)
{
    sim_irk_dims *dims = dims_;

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
    else if (!strcmp(field, "np"))
    {
        *value = dims->np;
    }
    else
    {
        printf("\nerror: sim_irk_dims_get: field not available: %s\n", field);
        exit(1);
    }
}



/************************************************
 * model
 ************************************************/

acados_size_t sim_irk_model_calculate_size(void *config, void *dims)
{
    acados_size_t size = 0;

    size += sizeof(irk_model);

    return size;
}



void *sim_irk_model_assign(void *config, void *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    irk_model *model = (irk_model *) c_ptr;
    c_ptr += sizeof(irk_model);

    model->impl_ode_fun = NULL;
    model->impl_ode_fun_jac_x_xdot_z = NULL;
    model->impl_ode_jac_x_xdot_u_z = NULL;
    model->impl_dae_jac_p = NULL;
    model->impl_ode_hess = NULL;

    assert((char *) raw_memory + sim_irk_model_calculate_size(config, dims) >= c_ptr);

    return model;
}



int sim_irk_model_set(void *model_, const char *field, void *value)
{
    irk_model *model = model_;

    if (!strcmp(field, "impl_ode_fun") || !strcmp(field, "impl_dae_fun"))
    {
        model->impl_ode_fun = value;
    }
    else if (!strcmp(field, "impl_ode_fun_jac_x_xdot_z") || !strcmp(field, "impl_dae_fun_jac_x_xdot_z") || !strcmp(field, "impl_ode_fun_jac_x_xdot") || !strcmp(field, "impl_dae_fun_jac_x_xdot"))
    {
        // TODO(oj): deprecate / remove the aliases without _z.
        model->impl_ode_fun_jac_x_xdot_z = value;
    }
    else if (!strcmp(field, "impl_ode_jac_x_xdot_u_z") || !strcmp(field, "impl_dae_jac_x_xdot_u_z") || !strcmp(field, "impl_ode_jac_x_xdot_u") || !strcmp(field, "impl_dae_jac_x_xdot_u"))
    {
        // TODO(oj): deprecate / remove the aliases without _z.
        model->impl_ode_jac_x_xdot_u_z = value;
    }
    else if (!strcmp(field, "impl_dae_jac_p"))
    {
        model->impl_dae_jac_p = value;
    }
    else if (!strcmp(field, "impl_ode_hes") || !strcmp(field, "impl_ode_hess") || !strcmp(field, "impl_dae_hess"))
    {
        model->impl_ode_hess = value;
    }
    else
    {
        printf("\nerror: sim_irk_model_set: wrong field: %s\n", field);
        exit(1);
    }

    return ACADOS_SUCCESS;
}



/************************************************
 * opts
 ************************************************/

acados_size_t sim_irk_opts_calculate_size(void *config_, void *dims)
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

void *sim_irk_opts_assign(void *config_, void *dims, void *raw_memory)
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

    assert((char *) raw_memory + sim_irk_opts_calculate_size(config_, dims) >= c_ptr);

    return (void *) opts;
}



void sim_irk_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    sim_irk_dims *dims = (sim_irk_dims *) dims_;
    sim_opts *opts = opts_;

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
    opts->newton_tol = 0.0;

    assert(opts->ns <= NS_MAX && "ns > NS_MAX!");

    // butcher tableau
    calculate_butcher_tableau(opts->ns, opts->collocation_type, opts->c_vec, opts->b_vec, opts->A_mat, opts->work);
    // for consistency check
    opts->tableau_size = opts->ns;
    opts->cost_computation = false;

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



void sim_irk_opts_update(void *config_, void *dims, void *opts_)
{
    sim_opts *opts = opts_;

    assert(opts->ns <= NS_MAX && "ns > NS_MAX!");

    calculate_butcher_tableau(opts->ns, opts->collocation_type, opts->c_vec, opts->b_vec, opts->A_mat, opts->work);

    opts->tableau_size = opts->ns;

    // for debugging: print butcher tableau
    // printf("Butcher tableau\n");
    // printf("\nc_vec:\n");
    // for (int i = 0; i < opts->ns; i++)
    // {
    //     printf("%f\t", opts->c_vec[i]);
    // }
    // printf("\nb_vec:\n");
    // for (int i = 0; i < opts->ns; i++)
    // {
    //     printf("%f\t", opts->b_vec[i]);
    // }
    // printf("\nA_mat:\n");
    // d_print_mat(opts->ns, opts->ns, opts->A_mat, opts->ns);

    return;
}



void sim_irk_opts_set(void *config_, void *opts_, const char *field, void *value)
{
    sim_opts *opts = (sim_opts *) opts_;
    sim_opts_set_(opts, field, value);
}



void sim_irk_opts_get(void *config_, void *opts_, const char *field, void *value)
{
    sim_opts *opts = (sim_opts *) opts_;
    sim_opts_get_(config_, opts, field, value);
}



/************************************************
 * memory
 ************************************************/

acados_size_t sim_irk_memory_calculate_size(void *config, void *dims_, void *opts_)
{
    // typecast
    sim_irk_dims *dims = (sim_irk_dims *) dims_;
    sim_opts *opts = opts_;

    // necessary integers
    int nx = dims->nx;
    int nz = dims->nz;
    int nu = dims->nu;

    acados_size_t size = sizeof(sim_irk_memory);

    size += nx * sizeof(double); // xdot
    size += nz * sizeof(double); // z
    size += 8;  // corresponds to memory alignment

    if (opts->cost_computation)
    {
        size += 1 * sizeof(struct blasfeo_dmat);  // cost_hess
        size += 1 * blasfeo_memsize_dmat(nx+nu, nx+nu);  // cost_hess
        size += 64;
    }

    if (opts->sens_forw_p)
    {
        size += 1 * blasfeo_memsize_dmat(nx, dims->np);
        size += 1 * sizeof(struct blasfeo_dmat);
        size += 64;
    }

    make_int_multiple_of(8, &size);

    return size;
}


void *sim_irk_memory_assign(void *config, void *dims_, void *opts_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    // typecast
    sim_irk_dims *dims = (sim_irk_dims *) dims_;
    sim_opts *opts = opts_;

    // necessary integers
    int nx = dims->nx;
    int nu = dims->nu;
    int nz = dims->nz;

    // struct
    sim_irk_memory *mem = (sim_irk_memory *) c_ptr;
    c_ptr += sizeof(sim_irk_memory);

    align_char_to(8, &c_ptr);
    if (opts->cost_computation)
    {
        assign_and_advance_blasfeo_dmat_structs(1, &mem->cost_hess, &c_ptr);
    }

    mem->S_p = NULL;

    if (opts->sens_forw_p)
    {
        assign_and_advance_blasfeo_dmat_structs(1, &mem->S_p, &c_ptr);
        align_char_to(64, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(nx, dims->np, mem->S_p, &c_ptr);
        blasfeo_dgese(nx, dims->np, 0.0, mem->S_p, 0, 0);
    }

    // assign doubles
    assign_and_advance_double(nz, &mem->z, &c_ptr);
    assign_and_advance_double(nx, &mem->xdot, &c_ptr);

    if (opts->cost_computation)
    {
        align_char_to(64, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(nx+nu, nx+nu, mem->cost_hess, &c_ptr);
    }

    // initialization of xdot, z is 0 if not changed
    for (int ii = 0; ii < nx; ii++)
        mem->xdot[ii] = 0.0;
    for (int ii = 0; ii < nz; ii++)
        mem->z[ii] = 0.0;

    return mem;
}



int sim_irk_memory_set(void *config_, void *dims_, void *mem_, const char *field, void *value)
{
    sim_config *config = config_;
    sim_irk_memory *mem = (sim_irk_memory *) mem_;

    int status = ACADOS_SUCCESS;

    if (!strcmp(field, "xdot"))
    {
        int nx;
        config->dims_get(config_, dims_, "nx", &nx);
        double *xdot = value;
        for (int ii=0; ii < nx; ii++)
            mem->xdot[ii] = xdot[ii];
    }
    else if (!strcmp(field, "z"))
    {
        int nz;
        config->dims_get(config_, dims_, "nz", &nz);
        double *z = value;
        for (int ii=0; ii < nz; ii++)
            mem->z[ii] = z[ii];
    }
    else if (!strcmp(field, "cost_capsule_ptr"))
    {
        mem->cost_capsule = value;
    }
    else if (!strcmp(field, "guesses_blasfeo"))
    {
        int nx, nz;
        config->dims_get(config_, dims_, "nx", &nx);
        config->dims_get(config_, dims_, "nz", &nz);

        struct blasfeo_dvec *sim_guess = (struct blasfeo_dvec *) value;
        blasfeo_unpack_dvec(nx, sim_guess, 0, mem->xdot, 1);
        blasfeo_unpack_dvec(nz, sim_guess, nx, mem->z, 1);
    }
    else
    {
        printf("sim_irk_memory_set: field %s is not supported! \n", field);
        exit(1);
    }

    return status;
}



int sim_irk_memory_set_to_zero(void *config_, void * dims_, void *opts_, void *mem_)
{
    sim_config *config = config_;
    sim_irk_memory *mem = (sim_irk_memory *) mem_;

    int status = ACADOS_SUCCESS;

    int nx, nz;
    config->dims_get(config_, dims_, "nz", &nz);
    config->dims_get(config_, dims_, "nx", &nx);

    for (int ii=0; ii < nz; ii++)
        mem->z[ii] = 0.0;
    for (int ii=0; ii < nx; ii++)
        mem->xdot[ii] = 0.0;

    return status;
}



void sim_irk_memory_get(void *config_, void *dims_, void *mem_, const char *field, void *value)
{
    sim_irk_memory *mem = mem_;

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
    else if (!strcmp(field, "cost_hess"))
    {
        struct blasfeo_dmat **ptr = value;
        *ptr = mem->cost_hess;
    }
    else if (!strcmp(field, "S_p"))
    {
        sim_irk_dims *dims = (sim_irk_dims *) dims_;

        if (dims->np == 0)
            return;

        if (mem->S_p == NULL)
        {
            printf("sim_irk_memory_get field %s requested but not allocated! Enable sens_forw_p.\n", field);
            exit(1);
        }

        blasfeo_unpack_dmat(dims->nx, dims->np, mem->S_p, 0, 0, value, dims->nx);
    }
    else
    {
        printf("sim_irk_memory_get field %s is not supported! \n", field);
        exit(1);
    }
}



/************************************************
 * workspace
 ************************************************/

acados_size_t sim_irk_workspace_calculate_size(void *config_, void *dims_, void *opts_)
{
    sim_irk_dims *dims = (sim_irk_dims *) dims_;
    sim_opts *opts = opts_;

    int ns = opts->ns;

    int nx = dims->nx;
    int nu = dims->nu;
    int nz = dims->nz;

    int nK = (nx + nz) * ns;

    int steps = opts->num_steps;

    acados_size_t size = sizeof(sim_irk_workspace);

    if (opts->sens_algebraic || opts->output_z)
    {
        size += (nx + nz) * sizeof(int);    // ipiv_one_stage
        size += ns * sizeof(double);        // Z_work
    }

    /* blasfeo structs */
    size += 4 * sizeof(struct blasfeo_dvec);          // rG, K, xt, xn

    if (opts->sens_adj || opts->sens_hess){
        size += 2 * steps * sizeof(struct blasfeo_dvec);  // **xn_traj, **K_traj
    }

    size += 2 * sizeof(struct blasfeo_dvec);  // lambda, lambdaK
    if (!opts->sens_hess)
    {
        size += 4 * sizeof(struct blasfeo_dmat);  // dG_dxu, dG_dK, dK_dxu, S_forw
    }
    else
    {
        size += (4 * steps + 1) * sizeof(struct blasfeo_dmat);  // dG_dxu, dG_dK, dK_dxu, S_forw
    }

    if (opts->sens_forw_p)
    {
        size += 2 * sizeof(struct blasfeo_dmat);
    }

    if (opts->cost_computation)
    {
        size += 1 * sizeof(struct blasfeo_dmat);  // S_forw_stage
    }

    /* blasfeo mem */
    if (opts->cost_computation)
    {
        size += 1 * blasfeo_memsize_dmat(nx, nx + nu);  // S_forw_stage
    }

    size += blasfeo_memsize_dvec(nK);   // K
    size += blasfeo_memsize_dvec(nK);   // rG
    size += 3 * blasfeo_memsize_dvec(nx);           // xt, xn, xtdot
    size += 1 * blasfeo_memsize_dvec(nx + nu);      // lambda
    size += 1 * blasfeo_memsize_dvec(nK);           // lambdaK

    if (!opts->sens_hess){
        size += 1 * blasfeo_memsize_dmat(nK, nx + nu);  // dG_dxu
        size += 1 * blasfeo_memsize_dmat(nK, nK);       // dG_dK
        size += 1 * blasfeo_memsize_dmat(nK, nx + nu);  // dK_dxu
        size += 1 * blasfeo_memsize_dmat(nx, nx + nu);  // S_forw
        size += nK * sizeof(int);  // ipiv
    }
    else
    {
        size += steps * blasfeo_memsize_dmat(nK, nx + nu);      // dG_dxu
        size += steps * blasfeo_memsize_dmat(nK, nK);           // dG_dK
        size += steps * blasfeo_memsize_dmat(nK, nx + nu);      // dK_dxu
        size += (steps + 1) * blasfeo_memsize_dmat(nx, nx + nu);      // S_forw
        size += steps * nK * sizeof(int);  // ipiv

        size += 1 * blasfeo_memsize_dmat(2*nx+nz+nu, 2*nx+nz+nu);  // f_hess
        size += 2 * blasfeo_memsize_dmat(2*nx+nz+nu, nx+nu);  // dxkzu_dw0, tmp_dxkzu_dw0
        size += 1 * blasfeo_memsize_dmat(nx + nu, nx + nu);  // Hess
    }

    if (opts->sens_forw_p) {
        size += 1 * blasfeo_memsize_dmat(nK, dims->np);
        size += 1 * blasfeo_memsize_dmat(nx + dims->nz, dims->np);
    }

    if ( opts->sens_adj || opts->sens_hess ){
        size += steps * blasfeo_memsize_dvec(nx);       // for xn_traj
        size += steps * blasfeo_memsize_dvec(nK);       // for K_traj
    }

    size += 2 * blasfeo_memsize_dmat(nx + nz, nx);  // df_dx, df_dxdot
    size += blasfeo_memsize_dmat(nx + nz, nu);      // df_du
    size += blasfeo_memsize_dmat(nx + nz, nz);      // df_dz

    if (opts->sens_algebraic && opts->exact_z_output)
    {
        size += blasfeo_memsize_dmat(nx + nz, nx + nz);  // df_dxdotz
        size += blasfeo_memsize_dmat(nx + nz, nx + nu);  // dk0_dxu
    }

    size += 1 * 8; // initial alignment
    make_int_multiple_of(64, &size);
    size += 1 * 64;

    return size;
}

static void *sim_irk_workspace_cast(void *config_, void *dims_, void *opts_, void *raw_memory)
{
    sim_opts *opts = opts_;
    sim_irk_dims *dims = (sim_irk_dims *) dims_;

    int ns = opts->ns;

    int nx = dims->nx;
    int nu = dims->nu;
    int nz = dims->nz;
    int nK = (nx + nz) * ns;

    int steps = opts->num_steps;

    char *c_ptr = (char *) raw_memory;

    // initial align
    align_char_to(8, &c_ptr);

    sim_irk_workspace *ws = (sim_irk_workspace *) c_ptr;
    c_ptr += sizeof(sim_irk_workspace);

    if ( opts->sens_adj || opts->sens_hess ){
        assign_and_advance_blasfeo_dvec_structs(steps, &ws->xn_traj, &c_ptr);
        assign_and_advance_blasfeo_dvec_structs(steps, &ws->K_traj, &c_ptr);
    }

    assign_and_advance_blasfeo_dvec_structs(1, &ws->rG, &c_ptr);
    assign_and_advance_blasfeo_dvec_structs(1, &ws->K, &c_ptr);
    assign_and_advance_blasfeo_dvec_structs(1, &ws->lambda, &c_ptr);
    assign_and_advance_blasfeo_dvec_structs(1, &ws->lambdaK, &c_ptr);

    // dG_dxu, dG_dK, dK_dxu, S_forw
    if (!opts->sens_hess){
        assign_and_advance_blasfeo_dmat_structs(1, &ws->dG_dxu, &c_ptr);
        assign_and_advance_blasfeo_dmat_structs(1, &ws->dG_dK, &c_ptr);
        assign_and_advance_blasfeo_dmat_structs(1, &ws->dK_dxu, &c_ptr);
        assign_and_advance_blasfeo_dmat_structs(1, &ws->S_forw, &c_ptr);
    } else {
        assign_and_advance_blasfeo_dmat_structs(steps, &ws->dG_dxu, &c_ptr);
        assign_and_advance_blasfeo_dmat_structs(steps, &ws->dG_dK, &c_ptr);
        assign_and_advance_blasfeo_dmat_structs(steps, &ws->dK_dxu, &c_ptr);
        assign_and_advance_blasfeo_dmat_structs(steps + 1, &ws->S_forw, &c_ptr);
    }
    if (opts->sens_forw_p)
    {
        assign_and_advance_blasfeo_dmat_structs(1, &ws->dK_dp, &c_ptr); // dK_dp
        assign_and_advance_blasfeo_dmat_structs(1, &ws->df_dp, &c_ptr); // df_dp
    }
    assign_and_advance_blasfeo_dvec_structs(1, &ws->xt, &c_ptr);
    assign_and_advance_blasfeo_dvec_structs(1, &ws->xn, &c_ptr);

    if (opts->cost_computation)
    {
        assign_and_advance_blasfeo_dmat_structs(1, &ws->S_forw_stage, &c_ptr);
    }

    /* algin c_ptr to 64 blasfeo_dmat_mem has to be assigned directly after that  */
    align_char_to(64, &c_ptr);

    if (opts->cost_computation)
    {
        assign_and_advance_blasfeo_dmat_mem(nx, nx+nu, ws->S_forw_stage, &c_ptr);
    }

    if (!opts->sens_hess){
        assign_and_advance_blasfeo_dmat_mem(nK, nx + nu, ws->dG_dxu, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(nK, nK,      ws->dG_dK, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(nK, nx + nu, ws->dK_dxu, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(nx, nx + nu, ws->S_forw, &c_ptr);
    }
    else
    {
        for (int ii = 0; ii < steps; ii++) {
            assign_and_advance_blasfeo_dmat_mem(nK, nx + nu, &ws->dG_dxu[ii], &c_ptr);
            assign_and_advance_blasfeo_dmat_mem(nK, nK, &ws->dG_dK[ii], &c_ptr);
            assign_and_advance_blasfeo_dmat_mem(nK, nx + nu, &ws->dK_dxu[ii], &c_ptr);
            assign_and_advance_blasfeo_dmat_mem(nx, nx + nu, &ws->S_forw[ii], &c_ptr);
        }
        assign_and_advance_blasfeo_dmat_mem(nx, nx + nu, &ws->S_forw[steps], &c_ptr);

        assign_and_advance_blasfeo_dmat_mem(2*nx+nz+nu, 2*nx+nz+nu, &ws->f_hess, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(2*nx+nz+nu, nx+nu, &ws->dxkzu_dw0, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(2*nx+nz+nu, nx+nu, &ws->tmp_dxkzu_dw0, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(nx + nu, nx + nu, &ws->Hess, &c_ptr);
    }

    if (opts->sens_forw_p)
    {
        assign_and_advance_blasfeo_dmat_mem(nK, dims->np, ws->dK_dp, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(nx + dims->nz, dims->np, ws->df_dp, &c_ptr);
    }

    assign_and_advance_blasfeo_dmat_mem(nx + nz, nx, &ws->df_dx, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx + nz, nx, &ws->df_dxdot, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx + nz, nu, &ws->df_du, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx + nz, nz, &ws->df_dz, &c_ptr);

    if (opts->sens_algebraic && opts->exact_z_output)
    {
        assign_and_advance_blasfeo_dmat_mem(nx + nz, nx + nz, &ws->df_dxdotz, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(nx + nz, nx + nu, &ws->dk0_dxu, &c_ptr);
    }

    assign_and_advance_blasfeo_dvec_mem(nK, ws->rG, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nK, ws->K, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nx, ws->xt, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nx, ws->xn, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nx, &ws->xtdot, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nx + nu, ws->lambda, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nK, ws->lambdaK, &c_ptr);


    if ( opts->sens_adj || opts->sens_hess ){
        for (int i = 0; i < steps; i++)
        {
            assign_and_advance_blasfeo_dvec_mem(nx, &ws->xn_traj[i], &c_ptr);
            assign_and_advance_blasfeo_dvec_mem(nK, &ws->K_traj[i], &c_ptr);
        }
    }

    if (opts->sens_algebraic || opts->output_z){
        assign_and_advance_double(ns, &ws->Z_work, &c_ptr);
        assign_and_advance_int((nx + nz), &ws->ipiv_one_stage, &c_ptr);
    }

    if (!opts->sens_hess){
        assign_and_advance_int(nK, &ws->ipiv, &c_ptr);
    } else {
        assign_and_advance_int(steps * nK, &ws->ipiv, &c_ptr);
    }

    // printf("\npointer moved - size calculated = %d bytes\n", c_ptr- (char*)raw_memory -
    // sim_irk_calculate_workspace_size(dims, opts_));

    assert((char *) raw_memory + sim_irk_workspace_calculate_size(config_, dims, opts_) >= c_ptr);

    return (void *) ws;
}


size_t sim_irk_get_external_fun_workspace_requirement(void *config_, void *dims_, void *opts_, void *model_)
{
    irk_model *model = model_;

    size_t size = 0;
    size_t tmp_size;

    tmp_size = external_function_get_workspace_requirement_if_defined(model->impl_ode_fun);
    size = size > tmp_size ? size : tmp_size;
    tmp_size = external_function_get_workspace_requirement_if_defined(model->impl_ode_fun_jac_x_xdot_z);
    size = size > tmp_size ? size : tmp_size;
    tmp_size = external_function_get_workspace_requirement_if_defined(model->impl_ode_hess);
    size = size > tmp_size ? size : tmp_size;
    tmp_size = external_function_get_workspace_requirement_if_defined(model->impl_ode_jac_x_xdot_u_z);
    size = size > tmp_size ? size : tmp_size;
    tmp_size = external_function_get_workspace_requirement_if_defined(model->impl_dae_jac_p);
    size = size > tmp_size ? size : tmp_size;

    return size;
}


void sim_irk_set_external_fun_workspaces(void *config_, void *dims_, void *opts_, void *model_, void *workspace_)
{
    irk_model *model = model_;

    external_function_set_fun_workspace_if_defined(model->impl_ode_fun, workspace_);
    external_function_set_fun_workspace_if_defined(model->impl_ode_fun_jac_x_xdot_z, workspace_);
    external_function_set_fun_workspace_if_defined(model->impl_ode_hess, workspace_);
    external_function_set_fun_workspace_if_defined(model->impl_dae_jac_p, workspace_);
    external_function_set_fun_workspace_if_defined(model->impl_ode_jac_x_xdot_u_z, workspace_);
}


int sim_irk_precompute(void *config_, sim_in *in, sim_out *out, void *opts_, void *mem_,
                       void *work_)
{
    return ACADOS_SUCCESS;
}



void sim_irk_compute_z_and_algebraic_sens(sim_irk_dims *dims, sim_opts *opts, sim_in *in, sim_out *out, sim_irk_memory *mem, sim_irk_workspace *ws, irk_model *model)
{
    int nx = dims->nx;
    int nz = dims->nz;
    int nu = dims->nu;
    int ns = opts->ns;

    double *u = in->u;
    double t0 = in->t0;

    double *Z_work = ws->Z_work;
    double *S_algebraic = out->S_algebraic;
    struct blasfeo_dvec *K = ws->K;

    struct blasfeo_dmat *df_dx = &ws->df_dx;
    struct blasfeo_dmat *df_dxdot = &ws->df_dxdot;
    struct blasfeo_dmat *df_du = &ws->df_du;
    struct blasfeo_dmat *df_dz = &ws->df_dz;

    struct blasfeo_dmat *df_dxdotz = &ws->df_dxdotz;
    struct blasfeo_dmat *dk0_dxu = &ws->dk0_dxu;

    int *ipiv_one_stage = ws->ipiv_one_stage;
    struct blasfeo_dvec *xtdot = &ws->xtdot;

    struct blasfeo_dvec *rG = ws->rG;

    struct blasfeo_dmat *dK_dxu = ws->dK_dxu;
    struct blasfeo_dmat *dK_dxu_ss;

    if (opts->sens_hess){
        dK_dxu_ss = &dK_dxu[0];
    }
    else
    {
        dK_dxu_ss = dK_dxu;
    }

    if (opts->sens_algebraic && !opts->exact_z_output)
    {
        double interpolated_value;
        for (int jj = 0; jj < nx+nu; jj++)
        {
            for (int ii = 0; ii < nz; ii++)
            {
                for (int kk = 0; kk < ns; kk++)
                {
                    Z_work[kk] = blasfeo_dgeex1(dK_dxu_ss, nx*ns+kk*nz+ii, jj);
                }
                neville_algorithm(0.0, ns - 1, opts->c_vec, Z_work, &interpolated_value);
                // eval polynomial through vals in Z_work at 0.
                S_algebraic[ii+jj*nz] = -interpolated_value;
            }
        }
    }

    if (opts->output_z || opts->sens_algebraic)
    {
        for (int ii = 0; ii < nz; ii++)
        {
            for (int jj = 0; jj < ns; jj++)
            {
                Z_work[jj] = blasfeo_dvecex1(K, nx * ns + nz * jj + ii);
                // copy values of z_ii in first step, into Z_work
            }
            neville_algorithm(0.0, ns - 1, opts->c_vec, Z_work, &out->zn[ii]);
            // eval polynomial through (c_jj, z_jj) at 0.
        }
    }

    if (opts->exact_z_output)
    {
        // INPUT: impl_ode
        // set input for impl_ode
        ws->impl_ode_type_in[0] = COLMAJ;
        ws->impl_ode_type_in[1] = BLASFEO_DVEC;
        ws->impl_ode_type_in[2] = COLMAJ;
        ws->impl_ode_type_in[3] = COLMAJ;
        ws->impl_ode_type_in[4] = COLMAJ;

        ws->impl_ode_in[0] = in->x;
        ws->impl_ode_in[1] = xtdot;
        ws->impl_ode_in[2] = u;
        ws->impl_ode_in[3] = &out->zn[0];
        ws->impl_ode_in[4] = &t0;

        // impl_ode_fun_jac_x_xdot_z
        ws->impl_ode_res_out.x = rG;

        ws->impl_ode_fun_jac_x_xdot_z_type_out[0] = BLASFEO_DVEC_ARGS;
        ws->impl_ode_fun_jac_x_xdot_z_out[0] = &ws->impl_ode_res_out;
        ws->impl_ode_fun_jac_x_xdot_z_type_out[1] = BLASFEO_DMAT;
        ws->impl_ode_fun_jac_x_xdot_z_out[1] = df_dx;
        ws->impl_ode_fun_jac_x_xdot_z_type_out[2] = BLASFEO_DMAT;
        ws->impl_ode_fun_jac_x_xdot_z_out[2] = df_dxdot;
        ws->impl_ode_fun_jac_x_xdot_z_type_out[3] = BLASFEO_DMAT;
        ws->impl_ode_fun_jac_x_xdot_z_out[3] = df_dz;

        // impl_ode_jac_x_xdot_u_z
        ws->impl_ode_jac_x_xdot_u_z_type_out[0] = BLASFEO_DMAT;
        ws->impl_ode_jac_x_xdot_u_z_out[0] = df_dx;
        ws->impl_ode_jac_x_xdot_u_z_type_out[1] = BLASFEO_DMAT;
        ws->impl_ode_jac_x_xdot_u_z_out[1] = df_dxdot;
        ws->impl_ode_jac_x_xdot_u_z_type_out[2] = BLASFEO_DMAT;
        ws->impl_ode_jac_x_xdot_u_z_out[2] = df_du;
        ws->impl_ode_jac_x_xdot_u_z_type_out[3] = BLASFEO_DMAT;
        ws->impl_ode_jac_x_xdot_u_z_out[3] = df_dz;

        if (opts->output_z || opts->sens_algebraic)
        {
            // initial guess for xdot0
            for (int ii = 0; ii < nx; ii++)
            {
                double interpolated_value;
                for (int jj = 0; jj < ns; jj++)
                {
                    // copy values of k_ii in first step, into Z_work
                    Z_work[jj] = blasfeo_dvecex1(K, nx * jj + ii);
                }
                neville_algorithm(0.0, ns - 1, opts->c_vec, Z_work, &interpolated_value);
                // eval polynomial through (c_jj, k_jj) at 0.
                blasfeo_pack_dvec(1, &interpolated_value, 1, xtdot, ii);
            }
            // perform extra newton iterations to get xdot0, z0 more precisely.
            for (int ii = 0; ii < opts->newton_iter; ii++)
            {
                if (ii == 0 || !opts->jac_reuse)
                {
                    // eval jacobians at interpolated values
                    acados_tic(&ws->timer_ad);
                    model->impl_ode_fun_jac_x_xdot_z->evaluate(
                        model->impl_ode_fun_jac_x_xdot_z, ws->impl_ode_type_in, ws->impl_ode_in,
                        ws->impl_ode_fun_jac_x_xdot_z_type_out, ws->impl_ode_fun_jac_x_xdot_z_out);
                    ws->timing_ad += acados_toc(&ws->timer_ad);

                    // set up df_dxdotz
                    blasfeo_dgecp(nx + nz, nx, df_dxdot, 0, 0, df_dxdotz, 0, 0);
                    blasfeo_dgecp(nx + nz, nz, df_dz,    0, 0, df_dxdotz, 0, nx);
                    // factorize
                    blasfeo_dgetrf_rp(nx + nz, nx + nz, df_dxdotz, 0, 0, df_dxdotz, 0, 0,
                                                                                ipiv_one_stage);
                }

                // permute rhs
                blasfeo_dvecpe(nx + nz, ipiv_one_stage, rG, 0);
                // backsolve
                blasfeo_dtrsv_lnu(nx + nz, df_dxdotz, 0, 0, rG, 0, rG, 0);
                blasfeo_dtrsv_unn(nx + nz, df_dxdotz, 0, 0, rG, 0, rG, 0);

                blasfeo_daxpy(nx, -1.0, rG, 0, xtdot, 0, xtdot, 0);
                blasfeo_dveccp(nx, rG, 0, xtdot, 0);
                blasfeo_unpack_dvec(nz, rG, nx, mem->z, 1);
                for (int jj = 0; jj < nz; jj++)
                {
                    out->zn[jj] -= mem->z[jj];
                }
            }
        }

        if (opts->sens_algebraic)
        {
            /* implicit function theorem to get S_alg */
            // eval jacobians at interpolated values
            acados_tic(&ws->timer_ad);
            model->impl_ode_jac_x_xdot_u_z->evaluate(
                    model->impl_ode_jac_x_xdot_u_z, ws->impl_ode_type_in, ws->impl_ode_in,
                    ws->impl_ode_jac_x_xdot_u_z_type_out, ws->impl_ode_jac_x_xdot_u_z_out);
            ws->timing_ad += acados_toc(&ws->timer_ad);

            // set up df_dxdotz
            blasfeo_dgecp(nx + nz, nx, df_dxdot, 0, 0, df_dxdotz, 0, 0);
            blasfeo_dgecp(nx + nz, nz, df_dz,    0, 0, df_dxdotz, 0, nx);
            // set up right hand side dk0_dxu
            blasfeo_dgecp(nx + nz, nx, df_dx, 0, 0, dk0_dxu, 0, 0);
            blasfeo_dgecp(nx + nz, nu, df_du, 0, 0, dk0_dxu, 0, nx);

            // solve linear system
            acados_tic(&ws->timer_la);
            blasfeo_dgetrf_rp(nx + nz, nx + nz, df_dxdotz, 0, 0, df_dxdotz, 0, 0, ipiv_one_stage);
            blasfeo_drowpe(nx + nz, ipiv_one_stage, dk0_dxu);
            blasfeo_dtrsm_llnu(nx + nz, nx + nu, 1.0, df_dxdotz, 0, 0,
                            dk0_dxu, 0, 0, dk0_dxu, 0, 0);
            blasfeo_dtrsm_lunn(nx + nz, nx + nu, 1.0, df_dxdotz, 0, 0,
                            dk0_dxu, 0, 0, dk0_dxu, 0, 0);
            ws->timing_la += acados_toc(&ws->timer_la);

            // solution has different sign
            blasfeo_dgesc(nx + nz, nx + nu, -1.0, dk0_dxu, 0, 0);

            // extract output
            blasfeo_unpack_dmat(nz, nx + nu, dk0_dxu, nx, 0, S_algebraic, nz);
        } // if sens_algebraic
    } // if exact_z_output
    out->info->LAtime += ws->timing_la;
    out->info->ADtime += ws->timing_ad;
}


/************************************************
 * integrator
 ************************************************/

// Note: These macros should be defined _only_ in this translation unit, so they are not prefixed.
//       We #undef them at the end of the file in case someone is doing the silly thing of #include-ing source.
#define UNPACK_DIMS_IRK(dims, opts) \
int nx = dims->nx;\
int nu = dims->nu;\
int nz = dims->nz;\
int np = dims->np;\
int nf_p = opts->sens_forw_p ? np : 0;\
int ns = opts->ns;\
int nK = (nx + nz) * ns;\

void sim_irk_initialize_functions(sim_irk_dims *dims, sim_opts *opts, sim_in *in, sim_out *out, sim_irk_memory *mem, sim_irk_workspace *ws, irk_model *model)
{
    // SET FUNCTION IN- & OUTPUT TYPES
    // INPUT: impl_ode
    ws->impl_ode_type_in[0] = BLASFEO_DVEC;      // xt
    ws->impl_ode_type_in[1] = BLASFEO_DVEC_ARGS; // k_i
    ws->impl_ode_type_in[2] = COLMAJ;            // u
    ws->impl_ode_type_in[3] = BLASFEO_DVEC_ARGS; // z_i
    ws->impl_ode_type_in[4] = COLMAJ;            // t

    ws->impl_ode_in[0] = ws->xt;                // 1st input is always xt
    ws->impl_ode_in[1] = &ws->impl_ode_xdot_in; // 2nd input is part of K[ss],
                                                // always update impl_ode_xdot_in
    ws->impl_ode_in[2] = in->u;                 // 3rd input is u (always)
    ws->impl_ode_in[3] = &ws->impl_ode_z_in;    // 4th input is part of Z[ss]
    ws->impl_ode_in[4] = &ws->t_current;        // current time

    // OUTPUT:
    // impl_ode_fun
    ws->impl_ode_fun_type_out[0] = BLASFEO_DVEC_ARGS;
    ws->impl_ode_fun_out[0] = &ws->impl_ode_res_out;
    ws->impl_ode_res_out.x = ws->rG;

    // impl_ode_fun_jac_x_xdot_z

    ws->impl_ode_fun_jac_x_xdot_z_type_out[0] = BLASFEO_DVEC_ARGS;
    ws->impl_ode_fun_jac_x_xdot_z_out[0] = &ws->impl_ode_res_out;
    ws->impl_ode_fun_jac_x_xdot_z_type_out[1] = BLASFEO_DMAT;
    ws->impl_ode_fun_jac_x_xdot_z_out[1] = &ws->df_dx;
    ws->impl_ode_fun_jac_x_xdot_z_type_out[2] = BLASFEO_DMAT;
    ws->impl_ode_fun_jac_x_xdot_z_out[2] = &ws->df_dxdot;
    ws->impl_ode_fun_jac_x_xdot_z_type_out[3] = BLASFEO_DMAT;
    ws->impl_ode_fun_jac_x_xdot_z_out[3] = &ws->df_dz;

    // impl_ode_jac_x_xdot_u_z
    ws->impl_ode_jac_x_xdot_u_z_type_out[0] = BLASFEO_DMAT;
    ws->impl_ode_jac_x_xdot_u_z_out[0] = &ws->df_dx;
    ws->impl_ode_jac_x_xdot_u_z_type_out[1] = BLASFEO_DMAT;
    ws->impl_ode_jac_x_xdot_u_z_out[1] = &ws->df_dxdot;
    ws->impl_ode_jac_x_xdot_u_z_type_out[2] = BLASFEO_DMAT;
    ws->impl_ode_jac_x_xdot_u_z_out[2] = &ws->df_du;
    ws->impl_ode_jac_x_xdot_u_z_type_out[3] = BLASFEO_DMAT;
    ws->impl_ode_jac_x_xdot_u_z_out[3] = &ws->df_dz;

    // impl_dae_jac_p
    ws->impl_dae_jac_p_type_out[0] = BLASFEO_DMAT;
    ws->impl_dae_jac_p_out[0] = ws->df_dp;

    // impl_ode_hess
    // INPUT: impl_ode_hess
    ws->impl_ode_hess_type_in[0] = BLASFEO_DVEC;            // xt
    ws->impl_ode_hess_in[0] = ws->xt;                       // 1st input is always xt
    ws->impl_ode_hess_type_in[1] = BLASFEO_DVEC_ARGS;       // k_i
    ws->impl_ode_hess_in[1] =  &ws->impl_ode_xdot_in;       // 2nd input is part of K[ss]
    ws->impl_ode_hess_type_in[2] = COLMAJ;                  // u
    ws->impl_ode_hess_in[2] = in->u;                        // 3rd input is u (always)
    ws->impl_ode_hess_type_in[3] = BLASFEO_DVEC_ARGS;       // z_i
    ws->impl_ode_hess_in[3] = &ws->impl_ode_z_in;           // 4th input is part of Z[ss]
    ws->impl_ode_hess_type_in[4] = BLASFEO_DVEC_ARGS;       // lambdaK component, direction
    ws->impl_ode_hess_in[4] = &ws->impl_ode_hess_lambda_in; // 5th input is part of lambdaK[ss]
    ws->impl_ode_hess_type_in[5] = COLMAJ;                  // t
    ws->impl_ode_hess_in[5] = &ws->t_current;               // current time

    // OUTPUT
    ws->impl_ode_hess_type_out[0] = BLASFEO_DMAT;
    ws->impl_ode_hess_out[0] = &ws->f_hess;
}

void sim_irk_initialize(sim_irk_dims *dims, sim_opts *opts, sim_in *in, sim_out *out, sim_irk_memory *mem, sim_irk_workspace *ws, irk_model *model)
{
    UNPACK_DIMS_IRK(dims, opts);

    /* Initialize functions */
    sim_irk_initialize_functions(dims, opts, in, out, mem, ws, model);

    /* Initialize & Pack */
    // initialize times
    out->info->LAtime = 0.0;
    out->info->ADtime = 0.0;

    if (nf_p > 0)
    {
        if (model->impl_dae_jac_p == NULL)
        {
            printf("sim IRK: impl_dae_jac_p is not provided but sens_forw_p=true.\n");
            exit(1);
        }
        blasfeo_dgese(nx, np, 0.0, mem->S_p, 0, 0);
    }

    blasfeo_dvecse(nK, 0.0, ws->lambdaK, 0);
    if (opts->sens_hess){
        blasfeo_dgese(nx + nu, nx + nu, 0.0, &ws->Hess, 0, 0);
    }
    blasfeo_pack_dvec(nx, in->x, 1, ws->xn, 0);
    blasfeo_pack_dmat(nx, nx + nu, in->S_forw, nx, ws->S_forw, 0, 0);
    blasfeo_pack_dvec(nx + nu, in->S_adj, 1, ws->lambda, 0); // TODO set to zero u-part ???
    // initialize integration variables
    for (int i = 0; i < ns; ++i)
    {
        // state derivatives
        blasfeo_pack_dvec(nx, mem->xdot, 1, ws->K, nx*i);
        // algebraic variables
        blasfeo_pack_dvec(nz, mem->z, 1, ws->K, nx*ns + i*nz);
    }
    // printf("sim_irk: K initialization\n");
    // blasfeo_print_exp_dvec(nK, K, 0);
    // exit(1);

    // initialize_cost_model
    ocp_nlp_cost_capsule *cost_capsule = mem->cost_capsule;
    if (opts->cost_computation)
    {
        ocp_nlp_cost_config *cost_config = cost_capsule->config;
        struct blasfeo_dvec *cost_grad = cost_config->memory_get(cost_capsule->memory, "grad");
        double *cost_fun = cost_config->memory_get(cost_capsule->memory, "fun");

        // initialize cost_fun, cost_grad, cost_hess
        blasfeo_dvecse(nx+nu, 0.0, cost_grad, 0);
        blasfeo_dgese(nx+nu, nx+nu, 0.0, mem->cost_hess, 0, 0);
        *cost_fun = 0.0;
        if (nz > 0)
        {
            printf("\nIRK cost_computation not implemented for nz>0!\n\n");
            exit(1);
        }
    }

    // initialize step pointers
    ws->dK_dxu_ss = ws->dK_dxu;
    ws->dG_dK_ss = ws->dG_dK;
    ws->dG_dxu_ss = ws->dG_dxu;
    ws->ipiv_ss = ws->ipiv;
    ws->S_forw_ss = ws->S_forw;
}

// void sim_irk_newton_step(sim_irk_dims *dims, sim_opts *opts, sim_in *in, sim_out *out, sim_irk_memory *mem, sim_irk_workspace *ws, irk_model *model, int ss)

void sim_irk_forward_step(sim_irk_dims *dims, sim_opts *opts, sim_in *in, sim_out *out, sim_irk_memory *mem, sim_irk_workspace *ws, irk_model *model, int ss)
{
    UNPACK_DIMS_IRK(dims, opts);
    int num_steps = opts->num_steps;
    double step = in->T / num_steps;

    // Cost integration
    ocp_nlp_cost_capsule *cost_capsule = mem->cost_capsule;

    // Stagewise pointers
    double a;
    // TODO(@anton) maybe worth it to have this live in ws?
    struct blasfeo_dmat *dG_dK_ss;
    struct blasfeo_dmat *dG_dxu_ss;
    struct blasfeo_dmat *dK_dxu_ss;
    struct blasfeo_dmat *S_forw_ss;
    int *ipiv_ss;

    // decide whether results from forward sensitivity propagation are stored,
    // or if memory has to be reused --> set pointers accordingly
    if (opts->sens_hess){
        ws->dK_dxu_ss = ws->dK_dxu+ss;
        ws->dG_dK_ss = ws->dG_dK+ss;
        ws->dG_dxu_ss = ws->dG_dxu+ss;
        ws->ipiv_ss = ws->ipiv+(ss*nK);
        ws->S_forw_ss = ws->S_forw+(ss+1);
        // copy current S_forw into S_forw_ss
        blasfeo_dgecp(nx, nx + nu, ws->S_forw+ss, 0, 0, ws->S_forw_ss, 0, 0);


        // copy last jacobian factorization into dG_dK_ss
        if (ss > 0 && opts->jac_reuse) {
            blasfeo_dgecp(nK, nK, ws->dG_dK+(ss-1), 0, 0, ws->dG_dK_ss, 0, 0);
            for (int ii = 0; ii < nK; ii++) {
                ws->ipiv_ss[ii] = ws->ipiv[nK*(ss-1) + ii];
            }
        }
    }

    if ( opts->sens_adj || opts->sens_hess )  // store current xn
        blasfeo_dveccp(nx, ws->xn, 0, ws->xn_traj+ss, 0);

    // do newton iters
    for (int iter = 0; iter < opts->newton_iter; iter++)
    {
        if ((opts->jac_reuse && (ss == 0) && (iter == 0)) || (!opts->jac_reuse))
        {
            // if new jacobian gets computed, initialize dG_dK_ss with zeros
            blasfeo_dgese(nK, nK, 0.0, ws->dG_dK_ss, 0, 0);
        }

        for (int ii = 0; ii < ns; ii++)
        {   // ii-th row of tableau
            // take x(n); copy a strvec into a strvec
            blasfeo_dveccp(nx, ws->xn, 0, ws->xt, 0);
            ws->t_current = in->t0 + ss * step + opts->c_vec[ii] * step;

            for (int jj = 0; jj < ns; jj++)
            {   // jj-th col of tableau
                // TODO(oj): precompute A_mat * step;
                a = opts->A_mat[ii + ns * jj] * step;
                // xt = xt + T_int * a[i,j]*K_j
                blasfeo_daxpy(nx, a, ws->K, jj * nx, ws->xt, 0, ws->xt, 0);
            }
            ws->impl_ode_xdot_in.xi = ii * nx;  // use k_i of K = (k_1,..., k_{ns},z_1,..., z_{ns})
            ws->impl_ode_z_in.xi = ns * nx + ii * nz;
            // use z_i of K = (k_1,..., k_{ns},z_1,..., z_{ns})
            ws->impl_ode_res_out.xi = ii * (nx + nz);  // store output in this position of rG

            // compute the residual of implicit ode at time t_ii
            if ((opts->jac_reuse && (ss == 0) && (iter == 0)) || (!opts->jac_reuse))
            {   // evaluate the ode function & jacobian w.r.t. x, xdot;
                // &  compute jacobian dG_dK_ss;
                acados_tic(&ws->timer_ad);
                model->impl_ode_fun_jac_x_xdot_z->evaluate(
                                                           model->impl_ode_fun_jac_x_xdot_z, ws->impl_ode_type_in, ws->impl_ode_in,
                                                           ws->impl_ode_fun_jac_x_xdot_z_type_out, ws->impl_ode_fun_jac_x_xdot_z_out);
                ws->timing_ad += acados_toc(&ws->timer_ad);

                // compute the blocks of dG_dK_ss
                for (int jj = 0; jj < ns; jj++)
                {  // compute the block (ii,jj)th block of dG_dK_ss
                    a = opts->A_mat[ii + ns * jj] * step;
                    blasfeo_dgead(nx + nz, nx, a, &ws->df_dx, 0, 0,
                                  ws->dG_dK_ss, ii * (nx + nz), jj * nx);
                    if (jj == ii)
                    {
                        blasfeo_dgead(nx + nz, nx, 1, &ws->df_dxdot, 0, 0,
                                      ws->dG_dK_ss, ii * (nx + nz), jj * nx);
                        blasfeo_dgead(nx + nz, nz, 1, &ws->df_dz,    0, 0,
                                      ws->dG_dK_ss, ii * (nx + nz), (nx * ns) + jj * nz);
                    }
                }  // end jj
            }
            else // only eval function (without jacobian)
            {
                acados_tic(&ws->timer_ad);
                model->impl_ode_fun->evaluate(model->impl_ode_fun, ws->impl_ode_type_in,
                                              ws->impl_ode_in, ws->impl_ode_fun_type_out,
                                              ws->impl_ode_fun_out);
                ws->timing_ad += acados_toc(&ws->timer_ad);
            }
        }  // end ii

        acados_tic(&ws->timer_la);
        // DGETRF computes an LU factorization of a general M-by-N matrix A
        // using partial pivoting with row interchanges.
        // printf("dG_dK_ss = (IRK) \n");
        // blasfeo_print_exp_dmat((nz+nx) *ns, (nz+nx) *ns, dG_dK_ss, 0, 0);
        if ((opts->jac_reuse && (ss == 0) && (iter == 0)) || (!opts->jac_reuse))
        {
            blasfeo_dgetrf_rp(nK, nK, ws->dG_dK_ss, 0, 0, ws->dG_dK_ss, 0, 0, ws->ipiv_ss);
        }

        // permute also the r.h.s
        blasfeo_dvecpe(nK, ws->ipiv_ss, ws->rG, 0);

        // solve dG_dK_ss * y = rG, dG_dK_ss on the (l)eft, (l)ower-trian, (n)o-trans
        // (u)nit trian
        blasfeo_dtrsv_lnu(nK, ws->dG_dK_ss, 0, 0, ws->rG, 0, ws->rG, 0);

        // solve dG_dK_ss * x = rG, dG_dK_ss on the (l)eft, (u)pper-trian, (n)o-trans
        // (n)o unit trian , and store x in rG
        blasfeo_dtrsv_unn(nK, ws->dG_dK_ss, 0, 0, ws->rG, 0, ws->rG, 0);

        ws->timing_la += acados_toc(&ws->timer_la);

        // scale and add a generic strmat into a generic strmat // K = K - rG, where rG is
        // [DeltaK, DeltaZ]
        blasfeo_daxpy(nK, -1.0, ws->rG, 0, ws->K, 0, ws->K, 0);

        // check early termination based on tolerance
        if (opts->newton_tol > 0)
        {
            blasfeo_dvecnrm_inf(nK, ws->rG, 0, &a);
            if (a < opts->newton_tol)
            {
                break;
            }
        }
    } // end newton_iter

    // save k vectors
    if ( opts->sens_adj || opts->sens_hess )
    {
        blasfeo_dveccp(nK, ws->K, 0, ws->K_traj+ss, 0);
    }

    // evaluate forward sensitivities
    if ( opts->sens_forw || opts->sens_hess || opts->sens_forw_p )
    {
        blasfeo_dgese(nK, nK, 0.0, ws->dG_dK_ss, 0, 0);
        // initialize dG_dK_ss with zeros
        // evaluate dG_dK_ss(xn,Kn)
        for (int ii = 0; ii < ns; ii++)
        {
            ws->impl_ode_xdot_in.xi = ii * nx;  // use k_i of K = (k_1,..., k_{ns},z_1,..., z_{ns})
            ws->impl_ode_z_in.xi    = ns * nx + ii * nz;
            // use z_i of K = (k_1,..., k_{ns},z_1,..., z_{ns})
            blasfeo_dveccp(nx, ws->xn, 0, ws->xt, 0);

            for (int jj = 0; jj < ns; jj++)
            {
                a = opts->A_mat[ii + ns * jj] * step;
                // xt = xt + T_int * a[i,j]*K_j
                blasfeo_daxpy(nx, a, ws->K, jj * nx, ws->xt, 0, ws->xt, 0);
            }
            ws->t_current = in->t0 + ss * step + opts->c_vec[ii] * step;

            acados_tic(&ws->timer_ad);
            model->impl_ode_jac_x_xdot_u_z->evaluate(
                                                     model->impl_ode_jac_x_xdot_u_z, ws->impl_ode_type_in, ws->impl_ode_in,
                                                     ws->impl_ode_jac_x_xdot_u_z_type_out, ws->impl_ode_jac_x_xdot_u_z_out);
            ws->timing_ad += acados_toc(&ws->timer_ad);

            blasfeo_dgecp(nx + nz, nx, &ws->df_dx, 0, 0, ws->dG_dxu_ss, ii * (nx + nz), 0);
            blasfeo_dgecp(nx + nz, nu, &ws->df_du, 0, 0, ws->dG_dxu_ss, ii * (nx + nz), nx);

            // compute the blocks of dG_dK_ss
            for (int jj = 0; jj < ns; jj++)
            {  // compute the block (ii,jj)th block of dG_dK_ss
                a = opts->A_mat[ii + ns * jj] * step;
                blasfeo_dgead(nx + nz, nx, a, &ws->df_dx, 0, 0,
                              ws->dG_dK_ss, ii * (nx + nz), jj * nx);
                if (jj == ii)
                {
                    blasfeo_dgead(nx + nz, nx, 1, &ws->df_dxdot, 0, 0,
                                  ws->dG_dK_ss, ii * (nx + nz), jj * nx);
                    blasfeo_dgead(nx + nz, nz, 1, &ws->df_dz,    0, 0,
                                  ws->dG_dK_ss, ii * (nx + nz), (nx * ns) + jj * nz);
                }
            }  // end jj
        }  // end ii

        // factorize dG_dK_ss
        acados_tic(&ws->timer_la);
        blasfeo_dgetrf_rp(nK, nK, ws->dG_dK_ss, 0, 0, ws->dG_dK_ss, 0, 0, ws->ipiv_ss);
        ws->timing_la += acados_toc(&ws->timer_la);

        if (opts->sens_forw || opts->sens_hess)
        {
            // obtain dK_dxu
            // set up right hand side
            if (in->identity_seed && ss == 0) // omit matrix multiplication for identity seed
                blasfeo_dgecp(nK, nx + nu, ws->dG_dxu_ss, 0, 0, ws->dK_dxu_ss, 0, 0);
            else
            {
                // dK_dw = 0 * dK_dw + 1 * dG_dx * S_forw_old
                blasfeo_dgemm_nn(nK, nx + nu, nx, 1.0, ws->dG_dxu_ss, 0, 0, ws->S_forw_ss, 0,
                                 0, 0.0, ws->dK_dxu_ss, 0, 0, ws->dK_dxu_ss, 0, 0);
                // printf("dG_dxu = \n");
                // blasfeo_print_exp_dmat(nx + nz, nx+nu, dG_dxu_ss, 0, 0);
                // dK_du = dK_du + 1 * dG_du
                blasfeo_dgead(nK, nu, 1.0, ws->dG_dxu_ss, 0, nx, ws->dK_dxu_ss, 0, nx);
            }
            // solve linear system
            acados_tic(&ws->timer_la);
            blasfeo_drowpe(nK, ws->ipiv_ss, ws->dK_dxu_ss);
            blasfeo_dtrsm_llnu(nK, nx + nu, 1.0, ws->dG_dK_ss, 0, 0, ws->dK_dxu_ss, 0, 0, ws->dK_dxu_ss, 0, 0);
            blasfeo_dtrsm_lunn(nK, nx + nu, 1.0, ws->dG_dK_ss, 0, 0, ws->dK_dxu_ss, 0, 0, ws->dK_dxu_ss, 0, 0);
            ws->timing_la += acados_toc(&ws->timer_la);
        }

        // printf("dK_dxu (solved) = (IRK, ss = %d) \n", ss);
        // blasfeo_print_exp_dmat(nK, nx + nu, dK_dxu_ss, 0, 0);

        if (nf_p > 0)
        {
            // dK_dp = dG_dx * S_p
            blasfeo_dgemm_nn(nK, dims->np, nx, 1.0, ws->dG_dxu_ss, 0, 0, mem->S_p, 0, 0, 0.0, ws->dK_dp, 0, 0, ws->dK_dp, 0, 0);

            for (int ii = 0; ii < ns; ii++)
            {
                // stage state: xt = xn + h * sum_j A_ij * k_j
                blasfeo_dveccp(nx, ws->xn, 0, ws->xt, 0);
                for (int jj = 0; jj < ns; jj++)
                {
                    a = opts->A_mat[ii + ns * jj] * step;
                    blasfeo_daxpy(nx, a, ws->K, jj * nx, ws->xt, 0, ws->xt, 0);
                }
                ws->impl_ode_xdot_in.xi = ii * nx;
                ws->impl_ode_z_in.xi    = ns * nx + ii * nz;
                ws->t_current = in->t0 + ss * step + opts->c_vec[ii] * step;

                // eval df/dp
                model->impl_dae_jac_p->evaluate(model->impl_dae_jac_p, ws->impl_ode_type_in, ws->impl_ode_in,
                                                ws->impl_dae_jac_p_type_out, ws->impl_dae_jac_p_out);

                // dK_dp += df/dp
                blasfeo_dgead(nx + nz, dims->np, 1.0, ws->df_dp, 0, 0, ws->dK_dp, ii * (nx + nz), 0);
            }

            // solve linear system
            blasfeo_drowpe(nK, ws->ipiv_ss, ws->dK_dp);
            blasfeo_dtrsm_llnu(nK, np, 1.0, ws->dG_dK_ss, 0, 0, ws->dK_dp, 0, 0, ws->dK_dp, 0, 0);
            blasfeo_dtrsm_lunn(nK, np, 1.0, ws->dG_dK_ss, 0, 0, ws->dK_dp, 0, 0, ws->dK_dp, 0, 0);

            // update S_p
            for (int jj = 0; jj < ns; jj++)
                blasfeo_dgead(nx, dims->np, -step * opts->b_vec[jj], ws->dK_dp, jj * nx, 0, mem->S_p, 0, 0);
        }

        if (opts->cost_computation)
        {
            ocp_nlp_cost_config *cost_config = cost_capsule->config;
            for (int ii = 0; ii < ns; ii++)
            {
                ws->impl_ode_z_in.xi = ns * nx + ii * nz;

                ws->t_current = in->t0 + ss * step + opts->c_vec[ii] * step;
                // compute x at stage (xt) and sensitivity (S_forw_stage)
                blasfeo_dveccp(nx, ws->xn, 0, ws->xt, 0);
                blasfeo_dgecp(nx, nx+nu, ws->S_forw_ss, 0, 0, ws->S_forw_stage, 0, 0);
                for (int jj = 0; jj < ns; jj++)
                {
                    a = opts->A_mat[ii + ns * jj] * step;
                    // xt = xt + T_int * a[i,j]*K_j
                    blasfeo_daxpy(nx, a, ws->K, jj * nx, ws->xt, 0, ws->xt, 0);
                    // NOTE(oj): dK_dxu_ss is actually -dK_dxu_ss, because alpha = -1.0 was not supported by blasfeo initially
                    blasfeo_dgead(nx, nx+nu, -a, ws->dK_dxu_ss, jj*nx, 0, ws->S_forw_stage, 0, 0);
                }
                cost_config->add_integrator_stage_cost_grad_hess(cost_capsule, ws->xt, in->u, &ws->impl_ode_z_in, ws->S_forw_stage, ws->t_current, opts->b_vec[ii]/num_steps, mem->cost_hess);
            }
        }

        // update forward sensitivity
        // NOTE(oj): dK_dxu_ss is actually -dK_dxu_ss, because alpha = -1.0
        // was not supported by blasfeos backsolve initially.
        if (opts->sens_forw || opts->sens_hess)
        {
            for (int jj = 0; jj < ns; jj++)
                blasfeo_dgead(nx, nx + nu, -step * opts->b_vec[jj], ws->dK_dxu_ss, jj * nx, 0,
                              ws->S_forw_ss, 0, 0);
        }
    }  // end if sens_forw || sens_hess || sens_forw_p
    else if (opts->cost_computation) // Cost computation without sensitivities
    {
        ocp_nlp_cost_config *cost_config = cost_capsule->config;
        for (int ii = 0; ii < ns; ii++)
        {
            ws->impl_ode_z_in.xi = ns * nx + ii * nz;

            ws->t_current = in->t0 + ss * step + opts->c_vec[ii] * step;
            // compute x at stage (xt)
            blasfeo_dveccp(nx, ws->xn, 0, ws->xt, 0);
            for (int jj = 0; jj < ns; jj++)
            {
                a = opts->A_mat[ii + ns * jj] * step;
                // xt = xt + T_int * a[i,j]*K_j
                blasfeo_daxpy(nx, a, ws->K, jj * nx, ws->xt, 0, ws->xt, 0);
            }

            cost_config->add_integrator_stage_cost(cost_capsule, ws->xt, in->u, &ws->impl_ode_z_in, ws->t_current, opts->b_vec[ii]/num_steps);
        }
    } // end cost_computation without sens


    // obtain x(n+1)
    for (int ii = 0; ii < ns; ii++){
        // xn += b_i * k_i
        blasfeo_daxpy(nx, step * opts->b_vec[ii], ws->K, ii * nx, ws->xn, 0, ws->xn, 0);
    }

    // algebraic variables output and corresponding sensitivity propagation
    if (ss == 0 && nz > 0)
    {
        sim_irk_compute_z_and_algebraic_sens(dims, opts, in, out, mem, ws, model);
    }

    if (ss == num_steps-1)
    {
        // store last xdot, z values for next initialization
        blasfeo_unpack_dvec(nx, ws->K, (ns-1) * nx, mem->xdot, 1);
        blasfeo_unpack_dvec(nz, ws->K, (ns-1) * nz + ns*nx, mem->z, 1);
    }
}

void sim_irk_forward_sweep(sim_irk_dims *dims, sim_opts *opts, sim_in *in, sim_out *out, sim_irk_memory *mem, sim_irk_workspace *ws, irk_model *model)
{
    UNPACK_DIMS_IRK(dims,opts);
    int num_steps = opts->num_steps;

    // set input for forward sweep
    ws->impl_ode_xdot_in.x = ws->K;
    ws->impl_ode_z_in.x = ws->K;

    // Integrate for num_steps steps
    for (int ss = 0; ss < num_steps; ss++)
    {
        sim_irk_forward_step(dims, opts, in, out, mem, ws, model, ss);
    }

    // extract results from forward sweep to output
    blasfeo_unpack_dvec(nx, ws->xn, 0, out->xn, 1);

    // Extract forward sensitivities
    if ( opts->sens_forw || opts->sens_hess )
        blasfeo_unpack_dmat(nx, nx + nu, ws->S_forw+(opts->sens_hess ? num_steps : 0), 0, 0, out->S_forw, nx);
}

void sim_irk_eval_and_factorize(sim_irk_dims *dims, sim_opts *opts, sim_in *in, sim_irk_workspace *ws, irk_model *model, int ss)
{
    UNPACK_DIMS_IRK(dims,opts);
    int num_steps = opts->num_steps;
    double step = in->T / num_steps;
    double a;

    blasfeo_dgese(nK, nK, 0.0, ws->dG_dK_ss, 0, 0);   // initialize dG_dK_ss with zeros
    /* evaluate function at stage i, build corresponding blocks of dG_dxu, dG_dK_ss */
    for (int ii = 0; ii < ns; ii++)
    {
        // set up input for impl_ode
        ws->impl_ode_xdot_in.xi = ii * nx;
        // use k_i of K = (k_1,..., k_{ns},z_1,..., z_{ns})
        ws->impl_ode_z_in.xi    = ns * nx + ii * nz;
        // use z_i of K = (k_1,..., k_{ns},z_1,..., z_{ns})
        ws->t_current = in->t0 + ss * step + opts->c_vec[ii] * step;

        // build stage value
        blasfeo_dveccp(nx, ws->xn_traj+ss, 0, ws->xt, 0);
        for (int jj = 0; jj < ns; jj++)
        {
            a = opts->A_mat[ii + ns * jj] * step;
            blasfeo_daxpy(nx, a, ws->K_traj+ss, jj * nx, ws->xt, 0, ws->xt, 0);
        }
        // set up input for impl_ode jacobians
        acados_tic(&ws->timer_ad);
        model->impl_ode_jac_x_xdot_u_z->evaluate(
                                                 model->impl_ode_jac_x_xdot_u_z, ws->impl_ode_type_in, ws->impl_ode_in,
                                                 ws->impl_ode_jac_x_xdot_u_z_type_out, ws->impl_ode_jac_x_xdot_u_z_out);
        ws->timing_ad += acados_toc(&ws->timer_ad);

        // build dG_dxu_ss
        blasfeo_dgecp(nx + nz, nx, &ws->df_dx, 0, 0, ws->dG_dxu_ss, ii * (nx + nz), 0);
        blasfeo_dgecp(nx + nz, nu, &ws->df_du, 0, 0, ws->dG_dxu_ss, ii * (nx + nz), nx);

        // build dG_dK_ss
        for (int jj = 0; jj < ns; jj++)
        {  // compute the block (ii,jj)th block of dG_dK_ss
            a = opts->A_mat[ii + ns * jj] * step;
            blasfeo_dgead(nx + nz, nx, a, &ws->df_dx, 0, 0,
                          ws->dG_dK_ss, ii * (nx + nz), jj * nx);
            if (jj == ii)
            {
                blasfeo_dgead(nx + nz, nx, 1.0, &ws->df_dxdot, 0, 0,
                              ws->dG_dK_ss, ii * (nx + nz), jj * nx);
                blasfeo_dgead(nx + nz, nz, 1.0, &ws->df_dz,    0, 0,
                              ws->dG_dK_ss, ii * (nx + nz), (nx * ns) + jj * nz);
            }
        }  // end jj
    }  // end ii

    // factorize dG_dK_ss - already done in forw if hessian is active
    acados_tic(&ws->timer_la);
    blasfeo_dgetrf_rp(nK, nK, ws->dG_dK_ss, 0, 0, ws->dG_dK_ss, 0, 0, ws->ipiv_ss);
    ws->timing_la += acados_toc(&ws->timer_la);
}

void sim_irk_propagate_hessian(sim_irk_dims *dims, sim_opts *opts, sim_in *in, sim_irk_workspace *ws, irk_model *model, int ss)
{
    UNPACK_DIMS_IRK(dims,opts);
    int num_steps = opts->num_steps;
    double step = in->T / num_steps;
    double a;

    ws->impl_ode_hess_lambda_in.x = ws->lambdaK;

    // evaluate second order derivatives and update Hessian
    // - HESSIAN PROPAGATION
    for (int ii = 0; ii < ns; ii++)
    {
        blasfeo_dgecp(nx, nx+nu, ws->S_forw_ss, 0, 0, &ws->dxkzu_dw0, 0, 0);
        // printf("dxkzu_dw0 = (IRK, ss = %d) \n", ss);
        // blasfeo_print_exp_dmat(2 * nx + nu + nz, nx + nu, dxkzu_dw0, 0, 0);
        // build stage value, and dxii_dw0
        blasfeo_dveccp(nx, ws->xn_traj+ss, 0, ws->xt, 0);
        for (int jj = 0; jj < ns; jj++)
        {
            a = opts->A_mat[ii + ns * jj] * step;
            blasfeo_daxpy(nx, a, ws->K_traj+ss, jj * nx, ws->xt, 0, ws->xt, 0);
            // dxii_dw0 += a * dkjj_dxu
            blasfeo_dgead(nx, nx + nu, -a, ws->dK_dxu_ss, jj * nx, 0, &ws->dxkzu_dw0, 0, 0);
        }
        // dk_dw0
        blasfeo_dgecpsc(nx, nx+nu, -1.0, ws->dK_dxu_ss, ii*nx, 0, &ws->dxkzu_dw0, nx, 0);
        // dz_dw0
        blasfeo_dgecpsc(nz, nx+nu, -1.0, ws->dK_dxu_ss, ns*nx+ii*nz, 0, &ws->dxkzu_dw0, 2*nx, 0);
        // du_dw0
        // TODO exploit the fact that this is [0, I] !!!
        blasfeo_dgese(nu, nx+nu, 0.0, &ws->dxkzu_dw0, 2*nx+nz, 0);
        blasfeo_ddiare(nu, 1.0, &ws->dxkzu_dw0, 2*nx+nz, nx);

        // printf("dxkzu_dw0 = (IRK, ss = %d) \n", ss);
        // blasfeo_print_exp_dmat(2 * nx + nu + nz, nx + nu, dxkzu_dw0, 0, 0);
        // printf("xt = (IRK, ss = %d) \n", ss);
        // blasfeo_print_exp_dvec(nx, xt, 0);
        // printf("xdot, z in = (IRK, ss = %d) \n", ss);
        // blasfeo_print_exp_dvec(nx, &K_traj[ss], 0);
        // blasfeo_print_exp_dvec(nz, &K_traj[ss], nx);
        // printf("lambda in = (IRK, ss = %d) \n", ss);
        // blasfeo_print_exp_dvec(nx + nz, lambdaK, 0);
        // set up input for impl_ode_hess
        ws->impl_ode_xdot_in.xi = ii * nx;
        // use k_i of K = (k_1,..., k_{ns},z_1,..., z_{ns})
        ws->impl_ode_z_in.xi    = ns * nx + ii * nz;
        // use z_i of K = (k_1,..., k_{ns},z_1,..., z_{ns})
        ws->impl_ode_hess_lambda_in.xi = ii * (nx + nz);

        ws->t_current = in->t0 + ss * step + opts->c_vec[ii] * step;

        // eval hessian function at stage ii
        // printf("dxkzu_dw0 = (IRK, ss = %d) \n", ss);
        // blasfeo_print_exp_dmat(2 * nx + nu + nz, nx + nu, dxkzu_dw0, 0, 0);
        acados_tic(&ws->timer_ad);

        model->impl_ode_hess->evaluate(model->impl_ode_hess, ws->impl_ode_hess_type_in,
                                       ws->impl_ode_hess_in, ws->impl_ode_hess_type_out, ws->impl_ode_hess_out);

        ws->timing_ad += acados_toc(&ws->timer_ad);
        // exploit that du_dw0 is [0, I]
        blasfeo_dgemm_nn(2*nx+nz+nu, nx+nu, 2*nx+nz, 1.0, &ws->f_hess, 0, 0, &ws->dxkzu_dw0, 0, 0, 0.0, &ws->tmp_dxkzu_dw0, 0, 0, &ws->tmp_dxkzu_dw0, 0, 0);
        blasfeo_dgead(2*nx+nz+nu, nu, 1.0, &ws->f_hess, 0, 2*nx+nz, &ws->tmp_dxkzu_dw0, 0, nx);
        blasfeo_dsyrk_ut(nx+nu, 2*nx+nz, 1.0, &ws->dxkzu_dw0, 0, 0, &ws->tmp_dxkzu_dw0, 0, 0, 1.0, &ws->Hess, 0, 0, &ws->Hess, 0, 0);
        blasfeo_dgead(nu, nx+nu, 1.0, &ws->tmp_dxkzu_dw0, 2*nx+nz, 0, &ws->Hess, nx, 0);
    }  // end for ii
}

void sim_irk_backward_step(sim_irk_dims *dims, sim_opts *opts, sim_in *in, sim_out *out, sim_irk_memory *mem, sim_irk_workspace *ws, irk_model *model, int ss)
{
    UNPACK_DIMS_IRK(dims,opts);
    int num_steps = opts->num_steps;
    double step = in->T / num_steps;

    // Stagewise pointers
    // TODO(@anton) maybe worth it to have this live in ws?
    if (opts->sens_hess){
        ws->dK_dxu_ss = ws->dK_dxu+ss;
        ws->dG_dK_ss = ws->dG_dK+ss;
        ws->dG_dxu_ss = ws->dG_dxu+ss;
        ws->ipiv_ss = ws->ipiv+(ss*nK);
        ws->S_forw_ss = ws->S_forw+ss;
        // lambdaK_ss = &lambdaK[ss];
        // lambda_ss_old = &lambda[ss+1];
        // lambda_ss = &lambda[ss];  // use something like this if lambda should be stored
    }
    ws->impl_ode_xdot_in.x = ws->K_traj+ss; // use K values of step ss
    ws->impl_ode_z_in.x = ws->K_traj+ss;    // use Z values of step ss

    /* evaluate impl_ode_jac_x_xdot_u_z -- build dG_dxu_ss, dG_dK_ss
       & factorize dG_dK_ss  */
    if ( !opts->sens_hess )
    {
        // TODO(@anton) a better name here!
        sim_irk_eval_and_factorize(dims, opts, in, ws, model, ss);
    }

    // update adjoint sensitivities: lambdaK
    // set up right hand side in vector lambdaK
    blasfeo_dvecse(nK, 0.0, ws->lambdaK, 0);
    for (int jj = 0; jj < ns; jj++)
        blasfeo_dveccpsc(nx, -step * opts->b_vec[jj], ws->lambda, 0, ws->lambdaK, jj * nx);
    // lambdaK_jj = -step b_jj * lambda_x

    // obtain lambdaK by solving linear system lambdaK <- (dG_dK)^(-T) lambdaK;
    acados_tic(&ws->timer_la);
    // dG_dK_ss - already factorized
    // solve linear system
    blasfeo_dtrsv_utn(nK, ws->dG_dK_ss, 0, 0, ws->lambdaK, 0, ws->lambdaK, 0);
    blasfeo_dtrsv_ltu(nK, ws->dG_dK_ss, 0, 0, ws->lambdaK, 0, ws->lambdaK, 0);
    blasfeo_dvecpei(nK, ws->ipiv_ss, ws->lambdaK, 0);
    ws->timing_la += acados_toc(&ws->timer_la);

    // update adjoint sensitivities lambda
    // lambda = 1 * lambda + 1 * dG_dxu_ss' * lambdaK
    blasfeo_dgemv_t(nK, nx + nu, 1.0, ws->dG_dxu_ss, 0, 0, ws->lambdaK, 0, 1.0, ws->lambda,
                    0, ws->lambda, 0);
    // Symmetric Hessian Propagation
    if ( opts->sens_hess )
    {
        sim_irk_propagate_hessian(dims, opts, in, ws, model, ss);
    } // end if ( opts->sens_hess )
}

void sim_irk_backward_sweep(sim_irk_dims *dims, sim_opts *opts, sim_in *in, sim_out *out, sim_irk_memory *mem, sim_irk_workspace *ws, irk_model *model)
{
    UNPACK_DIMS_IRK(dims,opts);

    for (int ss = opts->num_steps - 1; ss > -1; ss--)
    {
        sim_irk_backward_step(dims, opts, in, out, mem, ws, model, ss);
    }

    // Extract data
    blasfeo_unpack_dvec(nx + nu, ws->lambda, 0, out->S_adj, 1);
    if ( opts->sens_hess )
    {
        blasfeo_dtrtr_u(nu+nx, &ws->Hess, 0, 0, &ws->Hess, 0, 0);
        blasfeo_unpack_dmat(nx+nu, nx+nu, &ws->Hess, 0, 0, out->S_hess, nx + nu);
    }
}

int sim_irk(void *config_, sim_in *in, sim_out *out, void *opts_, void *mem_, void *work_)
{
    acados_timer timer;
    acados_tic(&timer);

    // Get variables from workspace, etc;
    // cast pointers
    sim_config *config = config_;
    sim_opts *opts = opts_;

    if ( opts->ns != opts->tableau_size )
    {
        printf("Error in sim_irk: the Butcher tableau size does not match ns");
        exit(1);
    }
    void *dims_ = in->dims;

    sim_irk_dims *dims = (sim_irk_dims *) dims_;
    sim_irk_workspace *ws =
        (sim_irk_workspace *) sim_irk_workspace_cast(config, dims, opts, work_);

    ws->timing_ad = 0.0;
    ws->timing_la = 0.0;

    sim_irk_memory *mem = (sim_irk_memory *) mem_;

    irk_model *model = in->model;

    if (model->impl_ode_fun == 0)
    {
        printf("sim IRK: impl_ode_fun is not provided. Exiting.\n");
        exit(1);
    }

    /************************************************
     * Initialize
     *       - Setup external function interfaces
     *       - Initialize workspace & memory
     ************************************************/
    sim_irk_initialize(dims, opts, in, out, mem, ws, model);

    // TODO(dimitris, FreyJo): implement NF (number of forward sensis) properly, instead of nx+nu?

    /************************************************
     * Forward Sweep
     *       - (simulation & forward sensitivities)
     ************************************************/
    sim_irk_forward_sweep(dims, opts, in, out, mem, ws, model);

    /*****************************************************************************
     * Backward Sweep
     *       - (adjoint sensitivities & hessian propagation)
     *       - hessian via symmetric forward-backward sweep
     *                    (see Algorithm 2 from Quirynen2016)
     *       - Quirynen2016: Symmetric Hessian propagation for lifted collocation integrators in direct optimal control
     *******************************************************************************/
    if ( opts->sens_adj  || opts->sens_hess )
    {
        sim_irk_backward_sweep(dims, opts, in, out, mem, ws, model);
    }

    out->info->CPUtime = acados_toc(&timer);
    // note: this is the time for factorization and solving the linear systems
    out->info->LAtime += ws->timing_la;
    out->info->ADtime += ws->timing_ad;

    mem->time_sim = out->info->CPUtime;
    mem->time_ad = out->info->ADtime;
    mem->time_la = out->info->LAtime;

    return ACADOS_SUCCESS;
}



void sim_irk_config_initialize_default(void *config_)
{
    sim_config *config = config_;

    config->evaluate = &sim_irk;
    config->precompute = &sim_irk_precompute;
    config->opts_calculate_size = &sim_irk_opts_calculate_size;
    config->opts_assign = &sim_irk_opts_assign;
    config->opts_initialize_default = &sim_irk_opts_initialize_default;
    config->opts_update = &sim_irk_opts_update;
    config->opts_set = &sim_irk_opts_set;
    config->opts_get = &sim_irk_opts_get;
    config->memory_calculate_size = &sim_irk_memory_calculate_size;
    config->memory_assign = &sim_irk_memory_assign;
    config->memory_set = &sim_irk_memory_set;
    config->memory_set_to_zero = &sim_irk_memory_set_to_zero;
    config->memory_get = &sim_irk_memory_get;
    config->workspace_calculate_size = &sim_irk_workspace_calculate_size;
    config->get_external_fun_workspace_requirement = &sim_irk_get_external_fun_workspace_requirement;
    config->set_external_fun_workspaces = &sim_irk_set_external_fun_workspaces;
    config->model_calculate_size = &sim_irk_model_calculate_size;
    config->model_assign = &sim_irk_model_assign;
    config->model_set = &sim_irk_model_set;
    config->dims_calculate_size = &sim_irk_dims_calculate_size;
    config->dims_assign = &sim_irk_dims_assign;
    config->dims_set = &sim_irk_dims_set;
    config->dims_get = &sim_irk_dims_get;
    return;
}

// #undef local macros
#undef UNPACK_DIMS_IRK
