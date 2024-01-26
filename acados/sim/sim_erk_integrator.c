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


// standard
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// acados
#include "acados/sim/sim_common.h"
#include "acados/sim/sim_collocation_utils.h"
#include "acados/sim/sim_erk_integrator.h"
#include "acados/utils/mem.h"

/************************************************
 * dims
 ************************************************/

acados_size_t sim_erk_dims_calculate_size()
{
    acados_size_t size = sizeof(sim_erk_dims);

    return size;
}



void *sim_erk_dims_assign(void *config_, void *raw_memory)
{
    char *c_ptr = raw_memory;

    sim_erk_dims *dims = (sim_erk_dims *) c_ptr;
    c_ptr += sizeof(sim_erk_dims);

    dims->nx = 0;
    dims->nu = 0;
    dims->nz = 0;

    assert((char *) raw_memory + sim_erk_dims_calculate_size() >= c_ptr);

    return dims;
}



void sim_erk_dims_set(void *config_, void *dims_, const char *field, const int *value)
{
    sim_erk_dims *dims = (sim_erk_dims *) dims_;

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
        if (*value != 0)
        {
            printf("\nerror: nz != 0\n");
            printf("algebraic variables not supported by ERK module\n");
            exit(1);
        }
    }
    else
    {
        printf("\nerror: sim_erk_dims_set: dim type not available: %s\n", field);
        exit(1);
    }
}



void sim_erk_dims_get(void *config_, void *dims_, const char *field, int *value)
{
    sim_erk_dims *dims = (sim_erk_dims *) dims_;

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
        *value = 0;
    }
    else
    {
        printf("\nerror: sim_erk_dims_get: dim type not available: %s\n", field);
        exit(1);
    }
}



/************************************************
 * model
 ************************************************/

acados_size_t sim_erk_model_calculate_size(void *config, void *dims)
{
    acados_size_t size = 0;

    size += sizeof(erk_model);

    return size;
}



void *sim_erk_model_assign(void *config, void *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    erk_model *model = (erk_model *) c_ptr;
    c_ptr += sizeof(erk_model);

    return model;
}



int sim_erk_model_set(void *model_, const char *field, void *value)
{
    erk_model *model = model_;

    if (!strcmp(field, "expl_ode_fun"))
    {
        model->expl_ode_fun = value;
    }
    else if (!strcmp(field, "expl_vde_for") || !strcmp(field, "expl_vde_forw"))
    {
        model->expl_vde_for = value;
    }
    else if (!strcmp(field, "expl_vde_adj"))
    {
        model->expl_vde_adj = value;
    }
    else if (!strcmp(field, "expl_ode_hes") || !strcmp(field, "expl_ode_hess"))
    {
        model->expl_ode_hes = value;
    }
    else
    {
        printf("\nerror: sim_erk_model_set: wrong field: %s\n", field);
        exit(1);
    }

    return ACADOS_SUCCESS;
}



/************************************************
 * opts
 ************************************************/

acados_size_t sim_erk_opts_calculate_size(void *config_, void *dims)
{
    int ns_max = NS_MAX;

    acados_size_t size = sizeof(sim_opts);

    size += ns_max * ns_max * sizeof(double);  // A_mat
    size += ns_max * sizeof(double);           // b_vec
    size += ns_max * sizeof(double);           // c_vec

    make_int_multiple_of(8, &size);
    size += 1 * 8;

    return size;
}



void *sim_erk_opts_assign(void *config_, void *dims, void *raw_memory)
{
    int ns_max = NS_MAX;

    char *c_ptr = (char *) raw_memory;

    sim_opts *opts = (sim_opts *) c_ptr;
    c_ptr += sizeof(sim_opts);

    align_char_to(8, &c_ptr);

    assign_and_advance_double(ns_max * ns_max, &opts->A_mat, &c_ptr);
    assign_and_advance_double(ns_max, &opts->b_vec, &c_ptr);
    assign_and_advance_double(ns_max, &opts->c_vec, &c_ptr);

    assert((char *) raw_memory + sim_erk_opts_calculate_size(config_, dims) >= c_ptr);

    opts->newton_iter = 0;
    // opts->scheme = NULL;
    opts->jac_reuse = false;

    return (void *) opts;
}



void sim_erk_opts_set(void *config_, void *opts_, const char *field, void *value)
{
    sim_opts *opts = (sim_opts *) opts_;
    sim_opts_set_(opts, field, value);
}



void sim_erk_opts_get(void *config_, void *opts_, const char *field, void *value)
{
    sim_opts *opts = (sim_opts *) opts_;
    sim_opts_get_(config_, opts, field, value);
}


void sim_erk_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    sim_opts *opts = opts_;
    sim_erk_dims *dims = (sim_erk_dims *) dims_;

    opts->ns = 4;  // ERK 4
    int ns = opts->ns;

    assert((ns == 1 || ns == 2 || ns == 3 || ns == 4) && "only number of stages = {1,2,3,4} implemented!");

    // set tableau size
    opts->tableau_size = opts->ns;

    double *A = opts->A_mat;
    double *b = opts->b_vec;
    double *c = opts->c_vec;

    get_explicit_butcher_tableau(ns, A, b, c);

    opts->num_steps = 1;
    opts->num_forw_sens = dims->nx + dims->nu;
    opts->sens_forw = true;
    opts->sens_adj = false;
    opts->sens_hess = false;
    opts->cost_computation = false;

    opts->output_z = false;
    opts->sens_algebraic = false;
}



void sim_erk_opts_update(void *config_, void *dims, void *opts_)
{
    sim_opts *opts = opts_;

    int ns = opts->ns;

    opts->tableau_size = opts->ns;

    assert((ns == 1 || ns == 2 || ns == 3 || ns == 4) && "only number of stages = {1,2,3,4} implemented!");

    assert(ns <= NS_MAX && "ns > NS_MAX!");

    // set tableau size
    opts->tableau_size = opts->ns;

    double *A = opts->A_mat;
    double *b = opts->b_vec;
    double *c = opts->c_vec;

    get_explicit_butcher_tableau(ns, A, b, c);

    return;
}



/************************************************
 * memory
 ************************************************/

acados_size_t sim_erk_memory_calculate_size(void *config, void *dims, void *opts_)
{
    acados_size_t size = sizeof(sim_erk_memory);

    return size;
}



void *sim_erk_memory_assign(void *config, void *dims, void *opts_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    sim_erk_memory *mem = (sim_erk_memory *) c_ptr;
    c_ptr += sizeof(sim_erk_memory);

    return mem;
}



int sim_erk_memory_set(void *config_, void *dims_, void *mem_, const char *field, void *value)
{
    printf("sim_erk_memory_set field %s is not supported! \n", field);
    exit(1);
}



int sim_erk_memory_set_to_zero(void *config_, void * dims_, void *opts_, void *mem_, const char *field)
{
    int status = ACADOS_SUCCESS;

    if (!strcmp(field, "guesses"))
    {
        // no guesses/initialization in ERK
    }
    else
    {
        printf("sim_erk_memory_set_to_zero field %s is not supported! \n", field);
        exit(1);
    }

    return status;
}



void sim_erk_memory_get(void *config_, void *dims_, void *mem_, const char *field, void *value)
{
    sim_erk_memory *mem = mem_;

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
        printf("sim_erk_memory_get field %s is not supported! \n", field);
        exit(1);
    }
}



/************************************************
 * workspace
 ************************************************/

acados_size_t sim_erk_workspace_calculate_size(void *config_, void *dims_, void *opts_)
{
    sim_opts *opts = opts_;
    sim_erk_dims *dims = (sim_erk_dims *) dims_;

    int ns = opts->ns;

    int nx = dims->nx;
    int nu = dims->nu;
    int nf = opts->num_forw_sens;

    int nX = nx * (1 + nf);  // (nx) for ODE and (nf*nx) for VDE
    int nhess = (nf + 1) * nf / 2;
    int num_steps = opts->num_steps;  // number of steps

    acados_size_t size = sizeof(sim_erk_workspace);

    size += (nX + nu) * sizeof(double);  // rhs_forw_in

    if (opts->sens_adj | opts->sens_hess)
    {
        size += num_steps * ns * nX * sizeof(double);   // K_traj
        size += (num_steps + 1) * nX * sizeof(double);  // out_forw_traj
    }
    else
    {
        size += ns * nX * sizeof(double);  // K_traj
        size += nX * sizeof(double);       // out_forw_traj
    }

    if (opts->sens_hess) // && opts->sens_adj)
    {
        size += (nX + nx + nu) * sizeof(double);          // rhs_adj_in
        size += (nx + nu + nhess) * sizeof(double);       // out_adj_tmp
        size += ns * (nx + nu + nhess) * sizeof(double);  // adj_traj
    }
    else if (opts->sens_adj)
    {
        size += (nx * 2 + nu) * sizeof(double);   // rhs_adj_in
        size += (nx + nu) * sizeof(double);       // out_adj_tmp
        size += ns * (nx + nu) * sizeof(double);  // adj_traj
    }

    make_int_multiple_of(8, &size);
    size += 1 * 8;

    return size;
}



static void *sim_erk_cast_workspace(void *config_, void *dims_, void *opts_, void *raw_memory, void *mem_)
{
    sim_opts *opts = opts_;
    sim_erk_dims *dims = (sim_erk_dims *) dims_;
    sim_erk_memory *mem = mem_;

    int ns = opts->ns;

    int nx = dims->nx;
    int nu = dims->nu;
    int nf = opts->num_forw_sens;

    int nX = nx * (1 + nf);  // (nx) for ODE and (nf*nx) for VDE
    int nhess = (nf + 1) * nf / 2;
    int num_steps = opts->num_steps;  // number of steps

    char *c_ptr = (char *) raw_memory;

    sim_erk_workspace *work = (sim_erk_workspace *) c_ptr;
    c_ptr += sizeof(sim_erk_workspace);

    align_char_to(8, &c_ptr);

    // switch to d_ptr
    double *d_ptr = (double *) c_ptr;

    //assign_and_advance_double(nX + nu, &workspace->rhs_forw_in, &c_ptr);
    work->rhs_forw_in = d_ptr;
    d_ptr += (nX+nu);

    if (opts->sens_adj | opts->sens_hess)
    {
        //
        //assign_and_advance_double(ns * num_steps * nX, &workspace->K_traj, &c_ptr);
        work->K_traj = d_ptr;
        d_ptr += ns*num_steps*nX;
        //assign_and_advance_double((num_steps + 1) * nX, &workspace->out_forw_traj, &c_ptr);
        work->out_forw_traj = d_ptr;
        d_ptr += (num_steps+1)*nX;
    }
    else
    {
        //assign_and_advance_double(ns * nX, &workspace->K_traj, &c_ptr);
        work->K_traj = d_ptr;
        d_ptr += ns*nX;
        //assign_and_advance_double(nX, &workspace->out_forw_traj, &c_ptr);
        work->out_forw_traj = d_ptr;
        d_ptr += nX;
    }

    if (opts->sens_hess) // && opts->sens_adj)
    {
        //assign_and_advance_double(nx + nX + nu, &workspace->rhs_adj_in, &c_ptr);
        work->rhs_adj_in = d_ptr;
        d_ptr += nx+nX+nu;
        //assign_and_advance_double(nx + nu + nhess, &workspace->out_adj_tmp, &c_ptr);
        work->out_adj_tmp = d_ptr;
        d_ptr += nx+nu+nhess;
        //assign_and_advance_double(ns * (nx + nu + nhess), &workspace->adj_traj, &c_ptr);
        work->adj_traj = d_ptr;
        d_ptr += ns*(nx+nu+nhess);
    }
    else if (opts->sens_adj)
    {
        //assign_and_advance_double((nx * 2 + nu), &workspace->rhs_adj_in, &c_ptr);
        work->rhs_adj_in = d_ptr;
        d_ptr += 2*nx+nu;
        //assign_and_advance_double(nx + nu, &workspace->out_adj_tmp, &c_ptr);
        work->out_adj_tmp = d_ptr;
        d_ptr += nu+nx;
        //assign_and_advance_double(ns * (nx + nu), &workspace->adj_traj, &c_ptr);
        work->adj_traj = d_ptr;
        d_ptr += ns*(nu+nx);
    }

    // update c_ptr
    c_ptr = (char *) d_ptr;

    assert((char *) raw_memory + mem->workspace_size >= c_ptr);

    return (void *) work;
}

/************************************************
 * functions
 ************************************************/

int sim_erk_precompute(void *config_, sim_in *in, sim_out *out, void *opts_, void *mem_,
                       void *work_)
{
    sim_erk_memory *mem = mem_;
    mem->workspace_size = sim_erk_workspace_calculate_size(config_, in->dims, opts_);
    return ACADOS_SUCCESS;
}



int sim_erk(void *config_, sim_in *in, sim_out *out, void *opts_, void *mem_, void *work_)
{
    acados_timer timer, timer_ad;
    acados_tic(&timer);

    sim_config *config = config_;
    sim_opts *opts = opts_;

    sim_erk_memory *mem = mem_;

    if ( opts->ns != opts->tableau_size )
    {
        printf("Error in sim_erk: the Butcher tableau size does not match ns\n");
        exit(1);
    }
    int ns = opts->ns;

    void *dims_ = in->dims;
    sim_erk_dims *dims = (sim_erk_dims *) dims_;

    sim_erk_workspace *work = sim_erk_cast_workspace(config, dims, opts, work_, mem_);

    int i, j, s, istep;
    double a = 0, b = 0;  // temp values of A_mat and b_vec
    int nx = dims->nx;
    int nu = dims->nu;
    int nz = dims->nz;

    // assert - only use supported features
    if (nz != 0)
    {
        printf("sim_erk: nz should be zero - DAEs are not supported by the ERK integrator\n");
        exit(1);
    }
    // if (opts->output_z)
    // {
    //     printf("sim_erk: opts->output_z should be false - DAEs are not supported for the ERK integrator\n");
    //     exit(1);
    // }
    if (opts->sens_algebraic)
    {
        printf("sim_erk: opts->sens_algebraic should be false - DAEs are not supported for the ERK integrator\n");
        exit(1);
    }

    int nf = opts->num_forw_sens;
    if (!opts->sens_forw) nf = 0;

    int nhess = (nf + 1) * nf / 2;
    int nX = nx + nx * nf;

    double *x = in->x;
    double *u = in->u;
    double *S_forw_in = in->S_forw;
    int num_steps = opts->num_steps;
    double step = in->T / num_steps;

    double *S_adj_in = in->S_adj;

    double *A_mat = opts->A_mat;
    double *b_vec = opts->b_vec;
    //    double *c_vec = opts->c_vec;

    double *K_traj = work->K_traj;
    double *forw_traj = work->out_forw_traj;
    double *rhs_forw_in = work->rhs_forw_in;

    double *adj_tmp = work->out_adj_tmp;
    double *adj_traj = work->adj_traj;
    double *rhs_adj_in = work->rhs_adj_in;

    double *xn = out->xn;
    double *S_forw_out = out->S_forw;
    double *S_adj_out = out->S_adj;
    double *S_hess_out = out->S_hess;

    ext_fun_arg_t ext_fun_type_in[5];
    void *ext_fun_in[5];
    ext_fun_arg_t ext_fun_type_out[3];
    void *ext_fun_out[3];

    ext_fun_arg_t expl_vde_type_in[4];
    void *expl_vde_in[4];
    ext_fun_arg_t expl_vde_type_out[3];
    void *expl_vde_out[3];

    int nx_squared_plus_nx = nx * nx + nx;
    int nx_times_nu = nx * nu;

    if (opts->sens_forw)
    {  // simulation + forward sensitivities
        expl_vde_type_in[0] = COLMAJ;
        expl_vde_in[0] = rhs_forw_in;  // x: nx
        expl_vde_type_in[1] = COLMAJ;
        expl_vde_in[1] = rhs_forw_in + nx;  // Sx: nx*nx
        expl_vde_type_in[2] = COLMAJ;
        expl_vde_in[2] = rhs_forw_in + nx_squared_plus_nx;  // Su: nx*nu
        expl_vde_type_in[3] = COLMAJ;
        expl_vde_in[3] = rhs_forw_in + nx_squared_plus_nx + nx_times_nu;  // u: nu

        expl_vde_type_out[0] = COLMAJ;
        expl_vde_type_out[1] = COLMAJ;
        expl_vde_type_out[2] = COLMAJ;
    }
    else
    {
        expl_vde_type_in[0] = COLMAJ;
        expl_vde_in[0] = rhs_forw_in;  // x: nx
        expl_vde_type_in[1] = COLMAJ;
        expl_vde_in[1] = rhs_forw_in + nx;  // u: nu

        expl_vde_type_out[0] = COLMAJ;
    }

    erk_model *model = in->model;

    double timing_ad = 0.0;

    /************************************************
     * forward sweep
     ************************************************/

    // initialize integrator variables
    for (i = 0; i < nx; i++) forw_traj[i] = x[i];  // x0
    if (opts->sens_forw)
    {
        for (i = 0; i < nx * nf; i++) forw_traj[nx + i] = S_forw_in[i];  // sensitivities
    }
    for (i = 0; i < nu; i++) rhs_forw_in[nX + i] = u[i];  // controls

    for (istep = 0; istep < num_steps; istep++)
    {
        if (opts->sens_adj | opts->sens_hess)
        {
            K_traj = work->K_traj + istep * ns * nX;
            forw_traj = work->out_forw_traj + (istep + 1) * nX;
            for (i = 0; i < nX; i++)
                forw_traj[i] = forw_traj[i - nX];
        }

        for (s = 0; s < ns; s++)
        {
            for (i = 0; i < nX; i++)
                rhs_forw_in[i] = forw_traj[i];
            for (j = 0; j < s; j++)
            {
                a = A_mat[j * ns + s];
                if (a != 0)
                {
                    a *= step;
                    for (i = 0; i < nX; i++)
                        rhs_forw_in[i] += a * K_traj[j * nX + i];
                }
            }

            acados_tic(&timer_ad);
            if (opts->sens_forw)
            {  // simulation + forward sensitivities
                // forward VDE evaluation
                expl_vde_out[0] = K_traj + s * nX;  // fun: nx
                expl_vde_out[1] = K_traj + s * nX + nx;  // Sx: nx*nx
                expl_vde_out[2] = K_traj + s * nX + nx_squared_plus_nx;  // Su: nx*nu
                model->expl_vde_for->evaluate(model->expl_vde_for, expl_vde_type_in, expl_vde_in,
                                              expl_vde_type_out, expl_vde_out);
            }
            else
            {  // simulation only
                if (model->expl_ode_fun == 0)
                {
                    printf("sim ERK: expl_ode_fun is not provided. Exiting.\n");
                    exit(1);
                }
                expl_vde_out[0] = K_traj + s * nX;  // fun: nx
                model->expl_ode_fun->evaluate(model->expl_ode_fun, expl_vde_type_in, expl_vde_in,
                                              expl_vde_type_out, expl_vde_out);  // ODE evaluation
            }
            timing_ad += acados_toc(&timer_ad);
        }
        for (s = 0; s < ns; s++)
        {
            b = step * b_vec[s];
            for (i = 0; i < nX; i++) forw_traj[i] += b * K_traj[s * nX + i];  // ERK step
        }
    }

    // store trajectory
    for (i = 0; i < nx; i++) xn[i] = forw_traj[i];
    // store forward sensitivities
    if (opts->sens_forw)
    {
        for (i = 0; i < nx * nf; i++) S_forw_out[i] = forw_traj[nx + i];
    }

    /************************************************
     * adjoint sweep
     ************************************************/
    if (opts->sens_adj | opts->sens_hess)
    {
        // initialize integrator variables
        for (i = 0; i < nx; i++)
            adj_tmp[i] = S_adj_in[i];
        // TODO read from S_adj also u-part ???
        for (i = 0; i < nu; i++)
            adj_tmp[nx + i] = 0.0;

        int nForw = nx;
        int nAdj = nx + nu;
        if (opts->sens_hess)
        {
            nForw = nX;
            nAdj = nx + nu + nhess;
            for (i = 0; i < nhess; i++)
                adj_tmp[nx + nu + i] = 0.0;
        }

        for (i = 0; i < nu; i++)
            rhs_adj_in[nForw + nx + i] = u[i];

        for (istep = num_steps - 1; istep >= 0; istep--)
        {

            K_traj = work->K_traj + istep * ns * nX;
            forw_traj = work->out_forw_traj + istep*nX;

            for (s = ns - 1; s >= 0; s--)
            {

                // forward variables:
                for (i = 0; i < nForw; i++)
                    rhs_adj_in[i] = forw_traj[i];  // extract x trajectory

                for (j = 0; j < s; j++)
                {
                    a = A_mat[j*ns + s];
                    if (a != 0)
                    {
                        a *= step;
                        for (i = 0; i < nForw; i++)
                            rhs_adj_in[i] += a * K_traj[j*nX + i];
                    }  // plus k traj
                }

                // adjoint variables:
                b = step * b_vec[s];

                for (i = 0; i < nx; i++) // why nx ???
                    rhs_adj_in[nForw + i] = b * adj_tmp[i];

                for (j = s + 1; j < ns; j++)
                {
                    a = A_mat[s*ns + j];
                    if (a != 0)
                    {
                        a *= step;
                        for (i = 0; i < nx; i++) // why nx ???
                            rhs_adj_in[nForw + i] += a * adj_traj[j*nAdj + i];
                    }
                }

                acados_tic(&timer_ad);
                if (!opts->sens_hess)
                {
                    ext_fun_type_in[0] = COLMAJ;
                    ext_fun_in[0] = rhs_adj_in + 0;  // x: nx

                    ext_fun_type_in[1] = COLMAJ;
                    ext_fun_in[1] = rhs_adj_in + nx;  // lam: nx
                    ext_fun_type_in[2] = COLMAJ;
                    ext_fun_in[2] = rhs_adj_in + nx + nx;  // u: nu

                    ext_fun_type_out[0] = COLMAJ;
                    ext_fun_out[0] = adj_traj + s * nAdj + 0;  // adj: nx+nu

                    if (model->expl_vde_adj == 0)
                    {
                        printf("sim ERK: expl_vde_adj is not provided. Exiting.\n");
                        exit(1);
                    }
                    // adjoint VDE evaluation
                    model->expl_vde_adj->evaluate(model->expl_vde_adj, ext_fun_type_in, ext_fun_in,
                            ext_fun_type_out, ext_fun_out);
                }
                else
                {
                    ext_fun_type_in[0] = COLMAJ;
                    ext_fun_in[0] = rhs_adj_in + 0;  // x: nx
                    ext_fun_type_in[1] = COLMAJ;
                    ext_fun_in[1] = rhs_adj_in + nx;  // Sx: nx*nx
                    ext_fun_type_in[2] = COLMAJ;
                    ext_fun_in[2] = rhs_adj_in + nx_squared_plus_nx;  // Su: nx*nu
                    ext_fun_type_in[3] = COLMAJ;
                    ext_fun_in[3] = rhs_adj_in + nx_squared_plus_nx + nx_times_nu;  // lam: nx
                    ext_fun_type_in[4] = COLMAJ;
                    ext_fun_in[4] = rhs_adj_in + nx_squared_plus_nx + nx_times_nu + nx;  // u: nu

                    ext_fun_type_out[0] = COLMAJ;
                    ext_fun_out[0] = adj_traj + s * nAdj + 0;  // adj: nx+nu
                    ext_fun_type_out[1] = COLMAJ;
                    ext_fun_out[1] = adj_traj + s * nAdj + nx + nu;  // hess: (nx+nu)*(nx+nu)

                    model->expl_ode_hes->evaluate(model->expl_ode_hes, ext_fun_type_in, ext_fun_in,
                            ext_fun_type_out, ext_fun_out);
                }
                timing_ad += acados_toc(&timer_ad);
            }

            for (s = 0; s < ns; s++)
                for (i = 0; i < nAdj; i++)
                    adj_tmp[i] += adj_traj[s*nAdj + i];  // ERK step

        }

        // store adjoint sensitivities
        for (i = 0; i < nx + nu; i++)
            S_adj_out[i] = adj_tmp[i];
        // store hessian
        if (opts->sens_hess)
        {
            // former line for tridiagonal export was
            //            for (i = 0; i < nhess; i++) S_hess_out[i] = adj_tmp[nx + nu + i];
            int count_upper = 0;
            for (int j = 0; j < nx + nu; j++)
            {
                for (int i = j; i < nx + nu; i++)
                {
                    S_hess_out[i + (nf) *j] = adj_tmp[nx + nu + count_upper];
                    // copy to upper part
                    S_hess_out[j + (nf) *i] = adj_tmp[nx + nu + count_upper];
                    count_upper++;
                }
            }
        }
    }

    // store timings
    out->info->CPUtime = acados_toc(&timer);
    out->info->LAtime = 0.0;
    out->info->ADtime = timing_ad;

    mem->time_sim = out->info->CPUtime;
    mem->time_ad = out->info->ADtime;
    mem->time_la = out->info->LAtime;

    return 0;
}


void sim_erk_config_initialize_default(void *config_)
{
    sim_config *config = config_;

    config->opts_calculate_size = &sim_erk_opts_calculate_size;
    config->opts_assign = &sim_erk_opts_assign;
    config->opts_initialize_default = &sim_erk_opts_initialize_default;
    config->opts_update = &sim_erk_opts_update;
    config->opts_set = &sim_erk_opts_set;
    config->opts_get = &sim_erk_opts_get;
    config->memory_calculate_size = &sim_erk_memory_calculate_size;
    config->memory_assign = &sim_erk_memory_assign;
    config->memory_set = &sim_erk_memory_set;
    config->memory_set_to_zero = &sim_erk_memory_set_to_zero;
    config->memory_get = &sim_erk_memory_get;
    config->workspace_calculate_size = &sim_erk_workspace_calculate_size;
    config->model_calculate_size = &sim_erk_model_calculate_size;
    config->model_assign = &sim_erk_model_assign;
    config->model_set = &sim_erk_model_set;
    config->evaluate = &sim_erk;
    config->precompute = &sim_erk_precompute;
    config->config_initialize_default = &sim_erk_config_initialize_default;
    config->dims_calculate_size = &sim_erk_dims_calculate_size;
    config->dims_assign = &sim_erk_dims_assign;
    config->dims_set = &sim_erk_dims_set;
    config->dims_get = &sim_erk_dims_get;
    return;
}
