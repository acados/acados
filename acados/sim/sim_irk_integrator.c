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
// acados
#include "acados/utils/mem.h"
#include "acados/utils/print.h"
#include "acados/utils/math.h"

#include "acados/sim/sim_common.h"

#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_blas.h"


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
    dims->ny = 0;

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
    else if (!strcmp(field, "ny"))
    {
        dims->ny = *value;
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

    irk_model *data = (irk_model *) c_ptr;
    c_ptr += sizeof(irk_model);

    assert((char *) raw_memory + sim_irk_model_calculate_size(config, dims) >= c_ptr);

    return data;
}



int sim_irk_model_set(void *model_, const char *field, void *value)
{
    irk_model *model = model_;

    if (!strcmp(field, "impl_ode_fun") || !strcmp(field, "impl_dae_fun"))
    {
        model->impl_ode_fun = value;
    }
    else if (!strcmp(field, "impl_ode_fun_jac_x_xdot") || !strcmp(field, "impl_dae_fun_jac_x_xdot"))
    {
        // TODO(oj): remove this case and fix dependencies
        model->impl_ode_fun_jac_x_xdot_z = value;
    }
    else if (!strcmp(field, "impl_ode_fun_jac_x_xdot_z") || !strcmp(field, "impl_dae_fun_jac_x_xdot_z"))
    {
        model->impl_ode_fun_jac_x_xdot_z = value;
    }
    else if (!strcmp(field, "impl_ode_jac_x_xdot_u") || !strcmp(field, "impl_dae_jac_x_xdot_u"))
    {
        // TODO(oj): remove this and update with z everywhere
        model->impl_ode_jac_x_xdot_u_z = value;
    }
    else if (!strcmp(field, "impl_ode_jac_x_xdot_u_z") || !strcmp(field, "impl_dae_jac_x_xdot_u_z"))
    {
        model->impl_ode_jac_x_xdot_u_z = value;
    }
    else if (!strcmp(field, "impl_ode_hes") || !strcmp(field, "impl_ode_hess") || !strcmp(field, "impl_dae_hess"))
    {
        model->impl_ode_hess = value;
    }
    else if (!strcmp(field, "nls_y_fun_jac") )
    {
        model->nls_y_fun_jac = value;
    }
    else if (!strcmp(field, "nls_y_fun") )
    {
        model->nls_y_fun = value;
    }
    else if (!strcmp(field, "conl_cost_fun_jac_hess") )
    {
        model->conl_cost_fun_jac_hess = value;
    }
    else if (!strcmp(field, "conl_cost_fun") )
    {
        model->conl_cost_fun = value;
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
    else if (!strcmp(field, "cost_fun"))
    {
        mem->cost_fun = value;
    }
    else if (!strcmp(field, "cost_grad"))
    {
        mem->cost_grad = value;
    }
    else if (!strcmp(field, "W_chol"))
    {
        mem->W_chol = value;
    }
    else if (!strcmp(field, "y_ref"))
    {
        mem->y_ref = value;
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



int sim_irk_memory_set_to_zero(void *config_, void * dims_, void *opts_, void *mem_, const char *field)
{
    sim_config *config = config_;
    sim_irk_memory *mem = (sim_irk_memory *) mem_;

    int status = ACADOS_SUCCESS;

    if (!strcmp(field, "guesses"))
    {
        int nx, nz;
        config->dims_get(config_, dims_, "nz", &nz);
        config->dims_get(config_, dims_, "nx", &nx);
        for (int ii=0; ii < nz; ii++)
            mem->z[ii] = 0.0;
        for (int ii=0; ii < nx; ii++)
            mem->xdot[ii] = 0.0;
    }
    else
    {
        printf("sim_irk_memory_set: field %s is not supported! \n", field);
        exit(1);
    }

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
    int ny = dims->ny;

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

    if (opts->cost_computation)
    {
        size += 4 * sizeof(struct blasfeo_dmat);  // J_y_tilde, tmp_nux_ny, S_forw_stage, tmp_nux_ny2
        size += 2 * sizeof(struct blasfeo_dvec);  // tmp_ny, nls_res
        if (opts->cost_type == CONVEX_OVER_NONLINEAR)
        {
            size += 3 * sizeof(struct blasfeo_dmat); // tmp_nv_ny, W, Jt_z
        }
    }

    /* blasfeo mem */
    if (opts->cost_computation)
    {
        size += 1 * blasfeo_memsize_dmat(ny, nx+nu);  // J_y_tilde
        size += 2 * blasfeo_memsize_dmat(nx+nu, ny);  // tmp_nux_ny, tmp_nux_ny2
        size += 2 * blasfeo_memsize_dvec(ny);  // tmp_ny, nls_res
        size += 1 * blasfeo_memsize_dmat(nx, nx + nu);  // S_forw_stage
        if (opts->cost_type == CONVEX_OVER_NONLINEAR)
        {
            size += 1 * blasfeo_memsize_dmat(nx + nu, ny);  // tmp_nv_ny
            size += 1 * blasfeo_memsize_dmat(ny, ny);  // W
            size += 1 * blasfeo_memsize_dmat(nz, ny);  // Jt_z
        }
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
    int ny = dims->ny;
    int nK = (nx + nz) * ns;

    int steps = opts->num_steps;

    char *c_ptr = (char *) raw_memory;

    // initial align
    align_char_to(8, &c_ptr);

    sim_irk_workspace *workspace = (sim_irk_workspace *) c_ptr;
    c_ptr += sizeof(sim_irk_workspace);

    if ( opts->sens_adj || opts->sens_hess ){
        assign_and_advance_blasfeo_dvec_structs(steps, &workspace->xn_traj, &c_ptr);
        assign_and_advance_blasfeo_dvec_structs(steps, &workspace->K_traj, &c_ptr);
    }

    assign_and_advance_blasfeo_dvec_structs(1, &workspace->rG, &c_ptr);
    assign_and_advance_blasfeo_dvec_structs(1, &workspace->K, &c_ptr);
    assign_and_advance_blasfeo_dvec_structs(1, &workspace->lambda, &c_ptr);
    assign_and_advance_blasfeo_dvec_structs(1, &workspace->lambdaK, &c_ptr);

    // dG_dxu, dG_dK, dK_dxu, S_forw
    if (!opts->sens_hess){
        assign_and_advance_blasfeo_dmat_structs(1, &workspace->dG_dxu, &c_ptr);
        assign_and_advance_blasfeo_dmat_structs(1, &workspace->dG_dK, &c_ptr);
        assign_and_advance_blasfeo_dmat_structs(1, &workspace->dK_dxu, &c_ptr);
        assign_and_advance_blasfeo_dmat_structs(1, &workspace->S_forw, &c_ptr);
    } else {
        assign_and_advance_blasfeo_dmat_structs(steps, &workspace->dG_dxu, &c_ptr);
        assign_and_advance_blasfeo_dmat_structs(steps, &workspace->dG_dK, &c_ptr);
        assign_and_advance_blasfeo_dmat_structs(steps, &workspace->dK_dxu, &c_ptr);
        assign_and_advance_blasfeo_dmat_structs(steps + 1, &workspace->S_forw, &c_ptr);
    }
    assign_and_advance_blasfeo_dvec_structs(1, &workspace->xt, &c_ptr);
    assign_and_advance_blasfeo_dvec_structs(1, &workspace->xn, &c_ptr);

    if (opts->cost_computation)
    {
        assign_and_advance_blasfeo_dmat_structs(1, &workspace->J_y_tilde, &c_ptr);
        assign_and_advance_blasfeo_dmat_structs(1, &workspace->tmp_nux_ny, &c_ptr);
        assign_and_advance_blasfeo_dmat_structs(1, &workspace->tmp_nux_ny2, &c_ptr);
        assign_and_advance_blasfeo_dmat_structs(1, &workspace->S_forw_stage, &c_ptr);
        if (opts->cost_type == CONVEX_OVER_NONLINEAR)
        {
            assign_and_advance_blasfeo_dmat_structs(1, &workspace->tmp_nv_ny, &c_ptr);
            assign_and_advance_blasfeo_dmat_structs(1, &workspace->W, &c_ptr);
            assign_and_advance_blasfeo_dmat_structs(1, &workspace->Jt_z, &c_ptr);
        }
        assign_and_advance_blasfeo_dvec_structs(1, &workspace->tmp_ny, &c_ptr);
        assign_and_advance_blasfeo_dvec_structs(1, &workspace->nls_res, &c_ptr);
    }

    /* algin c_ptr to 64 blasfeo_dmat_mem has to be assigned directly after that  */
    align_char_to(64, &c_ptr);

    if (opts->cost_computation)
    {
        assign_and_advance_blasfeo_dmat_mem(ny, nx+nu, workspace->J_y_tilde, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(nx+nu, ny, workspace->tmp_nux_ny, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(nx+nu, ny, workspace->tmp_nux_ny2, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(nx, nx+nu, workspace->S_forw_stage, &c_ptr);
        if (opts->cost_type == CONVEX_OVER_NONLINEAR)
        {
            assign_and_advance_blasfeo_dmat_mem(nx+nu, ny, workspace->tmp_nv_ny, &c_ptr);
            assign_and_advance_blasfeo_dmat_mem(ny, ny, workspace->W, &c_ptr);
            assign_and_advance_blasfeo_dmat_mem(nz, ny, workspace->Jt_z, &c_ptr);
        }
    }

    if (!opts->sens_hess){
        assign_and_advance_blasfeo_dmat_mem(nK, nx + nu, workspace->dG_dxu, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(nK, nK,      workspace->dG_dK, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(nK, nx + nu, workspace->dK_dxu, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(nx, nx + nu, workspace->S_forw, &c_ptr);
    }
    else
    {
        for (int ii = 0; ii < steps; ii++) {
            assign_and_advance_blasfeo_dmat_mem(nK, nx + nu, &workspace->dG_dxu[ii], &c_ptr);
            assign_and_advance_blasfeo_dmat_mem(nK, nK, &workspace->dG_dK[ii], &c_ptr);
            assign_and_advance_blasfeo_dmat_mem(nK, nx + nu, &workspace->dK_dxu[ii], &c_ptr);
            assign_and_advance_blasfeo_dmat_mem(nx, nx + nu, &workspace->S_forw[ii], &c_ptr);
        }
        assign_and_advance_blasfeo_dmat_mem(nx, nx + nu, &workspace->S_forw[steps], &c_ptr);

        assign_and_advance_blasfeo_dmat_mem(2*nx+nz+nu, 2*nx+nz+nu, &workspace->f_hess, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(2*nx+nz+nu, nx+nu, &workspace->dxkzu_dw0, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(2*nx+nz+nu, nx+nu, &workspace->tmp_dxkzu_dw0, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(nx + nu, nx + nu, &workspace->Hess, &c_ptr);
    }

    assign_and_advance_blasfeo_dmat_mem(nx + nz, nx, &workspace->df_dx, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx + nz, nx, &workspace->df_dxdot, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx + nz, nu, &workspace->df_du, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx + nz, nz, &workspace->df_dz, &c_ptr);

    if (opts->sens_algebraic && opts->exact_z_output)
    {
        assign_and_advance_blasfeo_dmat_mem(nx + nz, nx + nz, &workspace->df_dxdotz, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(nx + nz, nx + nu, &workspace->dk0_dxu, &c_ptr);
    }

    if (opts->cost_computation)
    {
        assign_and_advance_blasfeo_dvec_mem(ny, workspace->tmp_ny, &c_ptr);
        assign_and_advance_blasfeo_dvec_mem(ny, workspace->nls_res, &c_ptr);
    }

    assign_and_advance_blasfeo_dvec_mem(nK, workspace->rG, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nK, workspace->K, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nx, workspace->xt, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nx, workspace->xn, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nx, &workspace->xtdot, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nx + nu, workspace->lambda, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nK, workspace->lambdaK, &c_ptr);


    if ( opts->sens_adj || opts->sens_hess ){
        for (int i = 0; i < steps; i++)
        {
            assign_and_advance_blasfeo_dvec_mem(nx, &workspace->xn_traj[i], &c_ptr);
            assign_and_advance_blasfeo_dvec_mem(nK, &workspace->K_traj[i], &c_ptr);
        }
    }

    if (opts->sens_algebraic || opts->output_z){
        assign_and_advance_double(ns, &workspace->Z_work, &c_ptr);
        assign_and_advance_int((nx + nz), &workspace->ipiv_one_stage, &c_ptr);
    }

    if (!opts->sens_hess){
        assign_and_advance_int(nK, &workspace->ipiv, &c_ptr);
    } else {
        assign_and_advance_int(steps * nK, &workspace->ipiv, &c_ptr);
    }

    // printf("\npointer moved - size calculated = %d bytes\n", c_ptr- (char*)raw_memory -
    // sim_irk_calculate_workspace_size(dims, opts_));

    assert((char *) raw_memory + sim_irk_workspace_calculate_size(config_, dims, opts_) >= c_ptr);

    return (void *) workspace;
}



int sim_irk_precompute(void *config_, sim_in *in, sim_out *out, void *opts_, void *mem_,
                       void *work_)
{
    return ACADOS_SUCCESS;
}



/************************************************
 * integrator
 ************************************************/

int sim_irk(void *config_, sim_in *in, sim_out *out, void *opts_, void *mem_, void *work_)
{
    acados_timer timer, timer_ad, timer_la;
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
    int ns = opts->ns;

    void *dims_ = in->dims;
    sim_irk_dims *dims = (sim_irk_dims *) dims_;
    sim_irk_workspace *workspace =
        (sim_irk_workspace *) sim_irk_workspace_cast(config, dims, opts, work_);

    sim_irk_memory *mem = (sim_irk_memory *) mem_;

    irk_model *model = in->model;

    int nx = dims->nx;
    int nu = dims->nu;
    int nz = dims->nz;
    int ny = dims->ny;

    int nK = (nx + nz) * ns;

    double *u = in->u;
    double t0 = in->t0;
    double t_current;

    int newton_iter = opts->newton_iter;
    double *A_mat = opts->A_mat;
    double *b_vec = opts->b_vec;
    int num_steps = opts->num_steps;
    double step = in->T / num_steps;

    int *ipiv = workspace->ipiv;
    double *Z_work = workspace->Z_work;

    struct blasfeo_dmat *dG_dK = workspace->dG_dK;
    struct blasfeo_dvec *rG = workspace->rG;
    struct blasfeo_dvec *K = workspace->K;
    struct blasfeo_dmat *dG_dxu = workspace->dG_dxu;
    struct blasfeo_dmat *dK_dxu = workspace->dK_dxu;
    struct blasfeo_dvec *xt = workspace->xt;

    struct blasfeo_dvec *xn = workspace->xn;
    struct blasfeo_dmat *S_forw = workspace->S_forw;

    struct blasfeo_dmat *df_dx = &workspace->df_dx;
    struct blasfeo_dmat *df_dxdot = &workspace->df_dxdot;
    struct blasfeo_dmat *df_du = &workspace->df_du;
    struct blasfeo_dmat *df_dz = &workspace->df_dz;
    struct blasfeo_dmat *f_hess = &workspace->f_hess;
    struct blasfeo_dmat *dxkzu_dw0 = &workspace->dxkzu_dw0;
    struct blasfeo_dmat *tmp_dxkzu_dw0 = &workspace->tmp_dxkzu_dw0;

    struct blasfeo_dmat *df_dxdotz = &workspace->df_dxdotz;
    struct blasfeo_dmat *dk0_dxu = &workspace->dk0_dxu;
    struct blasfeo_dvec *xtdot = &workspace->xtdot;
    int *ipiv_one_stage = workspace->ipiv_one_stage;

    // for adjoint
    struct blasfeo_dvec *lambda = workspace->lambda;
    struct blasfeo_dvec *lambdaK = workspace->lambdaK;
    struct blasfeo_dvec *xn_traj = workspace->xn_traj;
    struct blasfeo_dvec *K_traj = workspace->K_traj;

    // for hessians only
    struct blasfeo_dmat *Hess = &workspace->Hess;

    double *x_out = out->xn;
    double *S_forw_out = out->S_forw;
    double *S_adj_out = out->S_adj;
    double *S_algebraic = out->S_algebraic;

    // for cost propagation only
    struct blasfeo_dvec *cost_grad = mem->cost_grad;
    struct blasfeo_dvec *nls_res = workspace->nls_res;
    struct blasfeo_dvec *tmp_ny = workspace->tmp_ny;

    struct blasfeo_dmat *cost_hess = mem->cost_hess;
    struct blasfeo_dmat *J_y_tilde = workspace->J_y_tilde;
    struct blasfeo_dmat *tmp_nux_ny = workspace->tmp_nux_ny;
    struct blasfeo_dmat *tmp_nux_ny2 = workspace->tmp_nux_ny2;
    // struct blasfeo_dmat *tmp_nx_nu = workspace->tmp_nx_nu;
    struct blasfeo_dmat *S_forw_stage = workspace->S_forw_stage;

    // declare
    double a;
    struct blasfeo_dmat *dG_dK_ss;
    struct blasfeo_dmat *dG_dxu_ss;
    struct blasfeo_dmat *dK_dxu_ss;
    struct blasfeo_dmat *S_forw_ss = S_forw;
    int *ipiv_ss;

    // SET FUNCTION IN- & OUTPUT TYPES
    // INPUT: impl_ode
    ext_fun_arg_t impl_ode_type_in[5];
    void *impl_ode_in[5];

    impl_ode_type_in[0] = BLASFEO_DVEC;       // xt
    impl_ode_type_in[1] = BLASFEO_DVEC_ARGS;  // k_i
    impl_ode_type_in[2] = COLMAJ;             // u
    impl_ode_type_in[3] = BLASFEO_DVEC_ARGS;  // z_i
    impl_ode_type_in[4] = COLMAJ;  // t

    struct blasfeo_dvec_args impl_ode_xdot_in;
    struct blasfeo_dvec_args impl_ode_z_in;

    impl_ode_in[0] = xt;                  // 1st input is always xt
    impl_ode_in[1] = &impl_ode_xdot_in;    // 2nd input is part of K[ss],
                                        // always update impl_ode_xdot_in
    impl_ode_in[2] = u;                 // 3rd input is u (always)
    impl_ode_in[3] = &impl_ode_z_in;    // 4th input is part of Z[ss]
    impl_ode_in[4] = &t_current;    // current time

    // OUTPUT:
    // impl_ode_fun
    ext_fun_arg_t impl_ode_fun_type_out[1];
    void *impl_ode_fun_out[1];
    impl_ode_fun_type_out[0] = BLASFEO_DVEC_ARGS;

    struct blasfeo_dvec_args impl_ode_res_out;
    impl_ode_res_out.x = rG;

    impl_ode_fun_out[0] = &impl_ode_res_out;

    // impl_ode_fun_jac_x_xdot_z
    ext_fun_arg_t impl_ode_fun_jac_x_xdot_z_type_out[4];
    void *impl_ode_fun_jac_x_xdot_z_out[4];
    impl_ode_fun_jac_x_xdot_z_type_out[0] = BLASFEO_DVEC_ARGS;
    impl_ode_fun_jac_x_xdot_z_out[0] = &impl_ode_res_out;
    impl_ode_fun_jac_x_xdot_z_type_out[1] = BLASFEO_DMAT;
    impl_ode_fun_jac_x_xdot_z_out[1] = df_dx;
    impl_ode_fun_jac_x_xdot_z_type_out[2] = BLASFEO_DMAT;
    impl_ode_fun_jac_x_xdot_z_out[2] = df_dxdot;
    impl_ode_fun_jac_x_xdot_z_type_out[3] = BLASFEO_DMAT;
    impl_ode_fun_jac_x_xdot_z_out[3] = df_dz;

    // impl_ode_jac_x_xdot_u_z
    ext_fun_arg_t impl_ode_jac_x_xdot_u_z_type_out[4];
    void *impl_ode_jac_x_xdot_u_z_out[4];
    impl_ode_jac_x_xdot_u_z_type_out[0] = BLASFEO_DMAT;
    impl_ode_jac_x_xdot_u_z_out[0] = df_dx;
    impl_ode_jac_x_xdot_u_z_type_out[1] = BLASFEO_DMAT;
    impl_ode_jac_x_xdot_u_z_out[1] = df_dxdot;
    impl_ode_jac_x_xdot_u_z_type_out[2] = BLASFEO_DMAT;
    impl_ode_jac_x_xdot_u_z_out[2] = df_du;
    impl_ode_jac_x_xdot_u_z_type_out[3] = BLASFEO_DMAT;
    impl_ode_jac_x_xdot_u_z_out[3] = df_dz;

    // impl_ode_hess
    // INPUT: impl_ode_hess
    ext_fun_arg_t impl_ode_hess_type_in[6];
    void *impl_ode_hess_in[6];

    struct blasfeo_dvec_args impl_ode_hess_lambda_in;

    impl_ode_hess_type_in[0] = BLASFEO_DVEC;           // xt
    impl_ode_hess_in[0] = xt;                          // 1st input is always xt
    impl_ode_hess_type_in[1] = BLASFEO_DVEC_ARGS;      // k_i
    impl_ode_hess_in[1] =  &impl_ode_xdot_in;          // 2nd input is part of K[ss]
    impl_ode_hess_type_in[2] = COLMAJ;                 // u
    impl_ode_hess_in[2] = u;                           // 3rd input is u (always)
    impl_ode_hess_type_in[3] = BLASFEO_DVEC_ARGS;      // z_i
    impl_ode_hess_in[3] = &impl_ode_z_in;              // 4th input is part of Z[ss]
    impl_ode_hess_type_in[4] = BLASFEO_DVEC_ARGS;      // lambdaK component, direction
    impl_ode_hess_in[4] = &impl_ode_hess_lambda_in; // 5th input is part of lambdaK[ss]
    impl_ode_hess_type_in[5] = COLMAJ;                 // t
    impl_ode_hess_in[5] = &t_current;                  // current time

    // OUTPUT
    ext_fun_arg_t impl_ode_hess_type_out[1];
    void *impl_ode_hess_out[1];
    impl_ode_hess_type_out[0] = BLASFEO_DMAT;
    impl_ode_hess_out[0] = f_hess;

    /* Initialize & Pack */
    // initialize
    double timing_ad = 0.0;
    double timing_la = 0.0;
    blasfeo_dvecse(nK, 0.0, lambdaK, 0);
    if (opts->sens_hess){
        blasfeo_dgese(nx + nu, nx + nu, 0.0, Hess, 0, 0);
    }

    if (opts->cost_computation)
    {
        // initialize cost_fun, cost_grad, cost_hess
        blasfeo_dvecse(nx+nu, 0.0, cost_grad, 0);
        blasfeo_dgese(nx+nu, nx+nu, 0.0, cost_hess, 0, 0);
        mem->cost_fun[0] = 0.0;
        if (nz > 0)
        {
            printf("\nIRK cost_computation not implemented for nz>0!\n\n");
            exit(1);
        }
    }

    // pack
    blasfeo_pack_dvec(nx, in->x, 1, xn, 0);
    blasfeo_pack_dmat(nx, nx + nu, in->S_forw, nx, S_forw, 0, 0);
    blasfeo_pack_dvec(nx + nu, in->S_adj, 1, lambda, 0); // TODO set to zero u-part ???

    // initialize integration variables
    for (int i = 0; i < ns; ++i)
    {
        // state derivatives
        blasfeo_pack_dvec(nx, mem->xdot, 1, K, nx*i);
        // algebraic variables
        blasfeo_pack_dvec(nz, mem->z, 1, K, nx*ns + i*nz);
    }
    // printf("sim_irk: K initialization\n");
    // blasfeo_print_exp_dvec(nK, K, 0);
    // exit(1);

    // TODO(dimitris, FreyJo): implement NF (number of forward sensis) properly, instead of nx+nu?

    /************************************************
    * Forward Sweep
    *       - (simulation & forward sensitivities)
    ************************************************/
    // set input for forward sweep
    impl_ode_xdot_in.x = K;
    impl_ode_z_in.x = K;

    // start the loop
    for (int ss = 0; ss < num_steps; ss++)
    {

        // decide whether results from forward sensitivity propagation are stored,
        // or if memory has to be reused --> set pointers accordingly
        if (opts->sens_hess){
            dK_dxu_ss = &dK_dxu[ss];
            dG_dK_ss = &dG_dK[ss];
            dG_dxu_ss = &dG_dxu[ss];
            ipiv_ss = &ipiv[ss*nK];
            S_forw_ss = &S_forw[ss+1];
            // copy current S_forw into S_forw_ss
            blasfeo_dgecp(nx, nx + nu, &S_forw[ss], 0, 0, S_forw_ss, 0, 0);


            // copy last jacobian factorization into dG_dK_ss
            if (ss > 0 && opts->jac_reuse) {
                blasfeo_dgecp(nK, nK, &dG_dK[ss-1], 0, 0, dG_dK_ss, 0, 0);
                for (int ii = 0; ii < nK; ii++) {
                    ipiv_ss[ii] = ipiv[nK*(ss-1) + ii];
                }
            }
        }
        else
        {
            dK_dxu_ss = dK_dxu;
            dG_dK_ss = dG_dK;
            dG_dxu_ss = dG_dxu;
            ipiv_ss = ipiv;
            S_forw_ss = S_forw;
        }

        if ( opts->sens_adj || opts->sens_hess )  // store current xn
            blasfeo_dveccp(nx, xn, 0, &xn_traj[ss], 0);

        for (int iter = 0; iter < newton_iter; iter++)
        {
            if ((opts->jac_reuse && (ss == 0) && (iter == 0)) || (!opts->jac_reuse))
            {
                // if new jacobian gets computed, initialize dG_dK_ss with zeros
                blasfeo_dgese(nK, nK, 0.0, dG_dK_ss, 0, 0);
            }

            for (int ii = 0; ii < ns; ii++)
            {  // ii-th row of tableau
                // take x(n); copy a strvec into a strvec
                blasfeo_dveccp(nx, xn, 0, xt, 0);
                t_current = t0 + ss * step + opts->c_vec[ii] * step;

                for (int jj = 0; jj < ns; jj++)
                {  // jj-th col of tableau
                    // TODO(oj): precompute A_mat * step;
                    a = A_mat[ii + ns * jj] * step;
                    // xt = xt + T_int * a[i,j]*K_j
                    blasfeo_daxpy(nx, a, K, jj * nx, xt, 0, xt, 0);
                }
                impl_ode_xdot_in.xi = ii * nx;  // use k_i of K = (k_1,..., k_{ns},z_1,..., z_{ns})
                impl_ode_z_in.xi = ns * nx + ii * nz;
                                // use z_i of K = (k_1,..., k_{ns},z_1,..., z_{ns})
                impl_ode_res_out.xi = ii * (nx + nz);  // store output in this position of rG

                // compute the residual of implicit ode at time t_ii
                if ((opts->jac_reuse && (ss == 0) && (iter == 0)) || (!opts->jac_reuse))
                {   // evaluate the ode function & jacobian w.r.t. x, xdot;
                    // &  compute jacobian dG_dK_ss;
                    acados_tic(&timer_ad);
                    model->impl_ode_fun_jac_x_xdot_z->evaluate(
                        model->impl_ode_fun_jac_x_xdot_z, impl_ode_type_in, impl_ode_in,
                        impl_ode_fun_jac_x_xdot_z_type_out, impl_ode_fun_jac_x_xdot_z_out);
                    timing_ad += acados_toc(&timer_ad);

                    // compute the blocks of dG_dK_ss
                    for (int jj = 0; jj < ns; jj++)
                    {  // compute the block (ii,jj)th block of dG_dK_ss
                        a = A_mat[ii + ns * jj] * step;
                        blasfeo_dgead(nx + nz, nx, a, df_dx, 0, 0,
                                            dG_dK_ss, ii * (nx + nz), jj * nx);
                        if (jj == ii)
                        {
                            blasfeo_dgead(nx + nz, nx, 1, df_dxdot, 0, 0,
                                          dG_dK_ss, ii * (nx + nz), jj * nx);
                            blasfeo_dgead(nx + nz, nz, 1, df_dz,    0, 0,
                                          dG_dK_ss, ii * (nx + nz), (nx * ns) + jj * nz);
                        }
                    }  // end jj
                }
                else // only eval function (without jacobian)
                {
                    if (model->impl_ode_fun == 0)
                    {
                        printf("sim IRK: impl_ode_fun is not provided. Exiting.\n");
                        exit(1);
                    }
                    acados_tic(&timer_ad);
                    model->impl_ode_fun->evaluate(model->impl_ode_fun, impl_ode_type_in,
                                                  impl_ode_in, impl_ode_fun_type_out,
                                                  impl_ode_fun_out);
                    timing_ad += acados_toc(&timer_ad);
                }
            }  // end ii

            acados_tic(&timer_la);
            // DGETRF computes an LU factorization of a general M-by-N matrix A
            // using partial pivoting with row interchanges.
            // printf("dG_dK_ss = (IRK) \n");
            // blasfeo_print_exp_dmat((nz+nx) *ns, (nz+nx) *ns, dG_dK_ss, 0, 0);
            if ((opts->jac_reuse && (ss == 0) && (iter == 0)) || (!opts->jac_reuse))
            {
                blasfeo_dgetrf_rp(nK, nK, dG_dK_ss, 0, 0, dG_dK_ss, 0, 0, ipiv_ss);
            }

            // permute also the r.h.s
            blasfeo_dvecpe(nK, ipiv_ss, rG, 0);

            // solve dG_dK_ss * y = rG, dG_dK_ss on the (l)eft, (l)ower-trian, (n)o-trans
            // (u)nit trian
            blasfeo_dtrsv_lnu(nK, dG_dK_ss, 0, 0, rG, 0, rG, 0);

            // solve dG_dK_ss * x = rG, dG_dK_ss on the (l)eft, (u)pper-trian, (n)o-trans
            // (n)o unit trian , and store x in rG
            blasfeo_dtrsv_unn(nK, dG_dK_ss, 0, 0, rG, 0, rG, 0);

            timing_la += acados_toc(&timer_la);

            // scale and add a generic strmat into a generic strmat // K = K - rG, where rG is
            // [DeltaK, DeltaZ]
            blasfeo_daxpy(nK, -1.0, rG, 0, K, 0, K, 0);

            // check early termination based on tolerance
            if (opts->newton_tol > 0)
            {
                blasfeo_dvecnrm_inf(nK, rG, 0, &a);
                if (a < opts->newton_tol)
                {
                    break;
                }
            }
        } // end newton_iter

        if ( opts->sens_adj || opts->sens_hess )
        {
            blasfeo_dveccp(nK, K, 0, &K_traj[ss], 0);
        }

        // evaluate forward sensitivities
        if ( opts->sens_forw || opts->sens_hess )
        {
            blasfeo_dgese(nK, nK, 0.0, dG_dK_ss, 0, 0);
            // initialize dG_dK_ss with zeros
            // evaluate dG_dK_ss(xn,Kn)
            for (int ii = 0; ii < ns; ii++)
            {
                impl_ode_xdot_in.xi = ii * nx;  // use k_i of K = (k_1,..., k_{ns},z_1,..., z_{ns})
                impl_ode_z_in.xi    = ns * nx + ii * nz;
                                                // use z_i of K = (k_1,..., k_{ns},z_1,..., z_{ns})
                blasfeo_dveccp(nx, xn, 0, xt, 0);

                for (int jj = 0; jj < ns; jj++)
                {
                    a = A_mat[ii + ns * jj] * step;
                    // xt = xt + T_int * a[i,j]*K_j
                    blasfeo_daxpy(nx, a, K, jj * nx, xt, 0, xt, 0);
                }
                t_current = t0 + ss * step + opts->c_vec[ii] * step;

                acados_tic(&timer_ad);
                model->impl_ode_jac_x_xdot_u_z->evaluate(
                    model->impl_ode_jac_x_xdot_u_z, impl_ode_type_in, impl_ode_in,
                    impl_ode_jac_x_xdot_u_z_type_out, impl_ode_jac_x_xdot_u_z_out);
                timing_ad += acados_toc(&timer_ad);

                blasfeo_dgecp(nx + nz, nx, df_dx, 0, 0, dG_dxu_ss, ii * (nx + nz), 0);
                blasfeo_dgecp(nx + nz, nu, df_du, 0, 0, dG_dxu_ss, ii * (nx + nz), nx);

                // compute the blocks of dG_dK_ss
                for (int jj = 0; jj < ns; jj++)
                {  // compute the block (ii,jj)th block of dG_dK_ss
                    a = A_mat[ii + ns * jj] * step;
                    blasfeo_dgead(nx + nz, nx, a, df_dx, 0, 0,
                                        dG_dK_ss, ii * (nx + nz), jj * nx);
                    if (jj == ii)
                    {
                        blasfeo_dgead(nx + nz, nx, 1, df_dxdot, 0, 0,
                                        dG_dK_ss, ii * (nx + nz), jj * nx);
                        blasfeo_dgead(nx + nz, nz, 1, df_dz,    0, 0,
                                        dG_dK_ss, ii * (nx + nz), (nx * ns) + jj * nz);
                    }
                }  // end jj
            }  // end ii

            // factorize dG_dK_ss
            acados_tic(&timer_la);
            blasfeo_dgetrf_rp(nK, nK, dG_dK_ss, 0, 0, dG_dK_ss, 0, 0, ipiv_ss);
            timing_la += acados_toc(&timer_la);

            // obtain dK_dxu
            // set up right hand side
            if (in->identity_seed && ss == 0) // omit matrix multiplication for identity seed
                blasfeo_dgecp(nK, nx + nu, dG_dxu_ss, 0, 0, dK_dxu_ss, 0, 0);
            else
            {
                // dK_dw = 0 * dK_dw + 1 * dG_dx * S_forw_old
                blasfeo_dgemm_nn(nK, nx + nu, nx, 1.0, dG_dxu_ss, 0, 0, S_forw_ss, 0,
                                    0, 0.0, dK_dxu_ss, 0, 0, dK_dxu_ss, 0, 0);
                // printf("dG_dxu = \n");
                // blasfeo_print_exp_dmat(nx + nz, nx+nu, dG_dxu_ss, 0, 0);
                // dK_du = dK_du + 1 * dG_du
                blasfeo_dgead(nK, nu, 1.0, dG_dxu_ss, 0, nx, dK_dxu_ss, 0, nx);
            }
            // solve linear system
            acados_tic(&timer_la);
            blasfeo_drowpe(nK, ipiv_ss, dK_dxu_ss);
            blasfeo_dtrsm_llnu(nK, nx + nu, 1.0, dG_dK_ss, 0, 0, dK_dxu_ss, 0, 0, dK_dxu_ss, 0, 0);
            blasfeo_dtrsm_lunn(nK, nx + nu, 1.0, dG_dK_ss, 0, 0, dK_dxu_ss, 0, 0, dK_dxu_ss, 0, 0);
            timing_la += acados_toc(&timer_la);

            // printf("dK_dxu (solved) = (IRK, ss = %d) \n", ss);
            // blasfeo_print_exp_dmat(nK, nx + nu, dK_dxu_ss, 0, 0);

            if (opts->cost_computation && opts->cost_type == NONLINEAR_LS)
            {
                ext_fun_arg_t nls_y_fun_jac_type_in[4];
                void *nls_y_fun_jac_in[4];
                ext_fun_arg_t nls_y_fun_jac_type_out[3];
                void *nls_y_fun_jac_out[3];

                nls_y_fun_jac_type_in[0] = BLASFEO_DVEC;
                nls_y_fun_jac_in[0] = xt;
                nls_y_fun_jac_type_in[1] = COLMAJ;
                nls_y_fun_jac_in[1] = u;
                nls_y_fun_jac_type_in[2] = BLASFEO_DVEC_ARGS;
                nls_y_fun_jac_in[2] = &impl_ode_z_in;
                nls_y_fun_jac_type_in[3] = COLMAJ;
                nls_y_fun_jac_in[3] = &t_current;

                impl_ode_z_in.x = &K_traj[ss];
                nls_y_fun_jac_type_out[0] = BLASFEO_DVEC;
                nls_y_fun_jac_out[0] = nls_res;  // fun: ny
                nls_y_fun_jac_type_out[1] = BLASFEO_DMAT;
                // nls_y_fun_jac_out[1] = &workspace->tmp_nux_ny;  // jac': (nu+nx) * ny
                nls_y_fun_jac_out[1] = tmp_nux_ny;  // jac': (nu+nx) * ny
                // dy_dux^T
                // TODO: set output in case of nz > 0
                nls_y_fun_jac_type_out[2] = BLASFEO_DMAT;
                nls_y_fun_jac_out[2] = tmp_nux_ny;   // jac_yexpr_z:  ny * nz

                for (int ii = 0; ii < ns; ii++)
                {
                    impl_ode_z_in.xi = ns * nx + ii * nz;

                    t_current = t0 + ss * step + opts->c_vec[ii] * step;
                    // compute x at stage (xt) and sensitivity (S_forw_stage)
                    blasfeo_dveccp(nx, xn, 0, xt, 0);
                    blasfeo_dgecp(nx, nx+nu, S_forw_ss, 0, 0, S_forw_stage, 0, 0);
                    for (int jj = 0; jj < ns; jj++)
                    {
                        a = A_mat[ii + ns * jj] * step;
                        // xt = xt + T_int * a[i,j]*K_j
                        blasfeo_daxpy(nx, a, K, jj * nx, xt, 0, xt, 0);
                        // NOTE(oj): dK_dxu_ss is actually -dK_dxu_ss, because alpha = -1.0 was not supported by blasfeo initially
                        blasfeo_dgead(nx, nx+nu, -a, dK_dxu_ss, jj*nx, 0, S_forw_stage, 0, 0);
                    }

                    model->nls_y_fun_jac->evaluate(model->nls_y_fun_jac, nls_y_fun_jac_type_in, nls_y_fun_jac_in,
                                        nls_y_fun_jac_type_out, nls_y_fun_jac_out);

                    // nls_res = nls_res - y_ref
                    blasfeo_daxpy(ny, -1.0, mem->y_ref, 0, nls_res, 0, nls_res, 0);

                    /* J_y_tilde = tmp_ny_nux * current_forward_sens( in [u,x] form) */
                    // NOTE: tmp_nux_ny = dy_dux^T here
                    // J_y_tilde[:nu, :ny] = tmp_nux_ny^T
                    blasfeo_dgetr(nu, ny, tmp_nux_ny, 0, 0, J_y_tilde, 0, 0);
                    // J_y_tilde[:nu, :ny] += tmp_nux_ny^T[nu:,:] * S_forw_stage[:, nx:]
                    blasfeo_dgemm_tn(ny, nu, nx, 1.0, tmp_nux_ny, nu, 0, S_forw_stage, 0, nx, 1.0, J_y_tilde, 0, 0,
                                    J_y_tilde, 0, 0);

                    // J_y_tilde (x part)
                    blasfeo_dgemm_tn(ny, nx, nx, 1.0, tmp_nux_ny, nu, 0, S_forw_stage, 0, 0, 0.0, J_y_tilde, 0, nu,
                                    J_y_tilde, 0, nu);

                    // transpose
                    blasfeo_dgetr(ny, nx+nu, J_y_tilde, 0, 0, tmp_nux_ny, 0, 0);

                    // tmp_nux_ny2 = W_chol * J_y_tilde (ny * (nx+nu))
                    blasfeo_dtrmm_rlnn(nu+nx, ny, 1.0, mem->W_chol, 0, 0, tmp_nux_ny, 0, 0,
                                       tmp_nux_ny2, 0, 0);

                    // tmp_ny = W_chol * nls_res
                    blasfeo_dtrmv_lnn(ny, mem->W_chol, 0, 0, nls_res, 0, tmp_ny, 0);

                    // cost_grad += b * tmp_ny^T * tmp_ny_nux = b * tmp_ny_nux^T * tmp_ny
                    blasfeo_dgemv_n(nx+nu, ny, b_vec[ii]/num_steps, tmp_nux_ny2, 0, 0, tmp_ny, 0,
                                    1.0, cost_grad, 0, cost_grad, 0);

                    // cost_hess += b * tmp_nux_ny_2 * tmp_nux_ny_2^T
                    // TODO: use syrk (exploit symmetry)
                    blasfeo_dgemm_nt(nx+nu, nx+nu, ny, b_vec[ii]/num_steps, tmp_nux_ny2, 0, 0, tmp_nux_ny2, 0, 0,
                                    1.0, cost_hess, 0, 0, cost_hess, 0, 0);

                    // cost function value
                    // NOTE: slack contribution and scaling done in cost module
                    mem->cost_fun[0] += 0.5 * b_vec[ii]/num_steps * blasfeo_ddot(ny, tmp_ny, 0, tmp_ny, 0);
                } // end ii
            } // end cost propagation NLS COST
            else if (opts->cost_computation && opts->cost_type == CONVEX_OVER_NONLINEAR)
            {
                // prepare function
                ext_fun_arg_t conl_fun_jac_hess_type_in[5];
                void *conl_fun_jac_hess_in[5];
                ext_fun_arg_t conl_fun_jac_hess_type_out[5];
                void *conl_fun_jac_hess_out[5];

                // inputs
                conl_fun_jac_hess_type_in[0] = BLASFEO_DVEC;
                conl_fun_jac_hess_in[0] = xt;
                conl_fun_jac_hess_type_in[1] = COLMAJ;
                conl_fun_jac_hess_in[1] = u;
                conl_fun_jac_hess_type_in[2] = BLASFEO_DVEC_ARGS;
                conl_fun_jac_hess_in[2] = &impl_ode_z_in;
                conl_fun_jac_hess_type_in[3] = BLASFEO_DVEC;
                conl_fun_jac_hess_in[3] = mem->y_ref;
                conl_fun_jac_hess_type_in[4] = COLMAJ;
                conl_fun_jac_hess_in[4] = &t_current;
                impl_ode_z_in.x = &K_traj[ss];

                // outputs // TODO!!
                conl_fun_jac_hess_type_out[0] = COLMAJ;
                conl_fun_jac_hess_out[0] = &a;         // fun: scalar
                conl_fun_jac_hess_type_out[1] = BLASFEO_DVEC;
                conl_fun_jac_hess_out[1] = workspace->tmp_ny;  // grad of outer loss wrt residual, ny
                conl_fun_jac_hess_type_out[2] = BLASFEO_DMAT;
                conl_fun_jac_hess_out[2] = tmp_nux_ny;  // inner Jacobian wrt ux, transposed, (nu+nx) x ny
                conl_fun_jac_hess_type_out[3] = BLASFEO_DMAT;
                conl_fun_jac_hess_out[3] = workspace->Jt_z; // inner Jacobian wrt z, transposed, nz x ny
                conl_fun_jac_hess_type_out[4] = BLASFEO_DMAT;
                conl_fun_jac_hess_out[4] = workspace->W;    // outer hessian: ny x ny

                for (int ii = 0; ii < ns; ii++)
                {
                    impl_ode_z_in.xi = ns * nx + ii * nz;

                    t_current = t0 + ss * step + opts->c_vec[ii] * step;
                    // compute x at stage (xt) and sensitivity (S_forw_stage)
                    blasfeo_dveccp(nx, xn, 0, xt, 0);
                    blasfeo_dgecp(nx, nx+nu, S_forw_ss, 0, 0, S_forw_stage, 0, 0);
                    for (int jj = 0; jj < ns; jj++)
                    {
                        a = A_mat[ii + ns * jj] * step;
                        // xt = xt + T_int * a[i,j]*K_j
                        blasfeo_daxpy(nx, a, K, jj * nx, xt, 0, xt, 0);
                        // NOTE(oj): dK_dxu_ss is actually -dK_dxu_ss, because alpha = -1.0 was not supported by blasfeo initially
                        blasfeo_dgead(nx, nx+nu, -a, dK_dxu_ss, jj*nx, 0, S_forw_stage, 0, 0);
                    }

                    // evaluate external function
                    model->conl_cost_fun_jac_hess->evaluate(model->conl_cost_fun_jac_hess, conl_fun_jac_hess_type_in,
                                                conl_fun_jac_hess_in, conl_fun_jac_hess_type_out, conl_fun_jac_hess_out);
                    // hessian of outer loss function
                    blasfeo_dpotrf_l(ny, workspace->W, 0, 0, mem->W_chol, 0, 0);

                    if (nz > 0) // TODO: test this!
                    {
                        // // Jt_ux_tilde = workspace->tmp_nux_ny + dzdux_tran*Jt_z
                        // blasfeo_dgemm_nn(nu + nx, ny, nz, 1.0, memory->dzdux_tran, 0, 0,
                        //         &workspace->Jt_z, 0, 0, 1.0, &workspace->tmp_nux_ny, 0, 0, &Jt_ux_tilde, 0, 0);

                        // // cost_grad += b * Jt_ux_tilde * tmp_ny
                        // blasfeo_dgemv_n(nu+nx, ny, b_vec[ii]/num_steps, &Jt_ux_tilde, 0, 0, tmp_ny, 0,
                        //                 1.0, cost_grad, 0, cost_grad, 0);

                        // // tmp_nv_ny = Jt_ux_tilde * W_chol
                        // blasfeo_dtrmm_rlnn(nu + nx, ny, 1.0, mem->W_chol, 0, 0,
                        //                 &Jt_ux_tilde, 0, 0, worksace->tmp_nv_ny, 0, 0);
                    }
                    else
                    {
                        /* J_y_tilde = tmp_ny_nux * current_forward_sens( in [u,x] form) */
                        // NOTE: tmp_nux_ny = dy_dux^T here
                        // J_y_tilde[:nu, :ny] = tmp_nux_ny^T
                        blasfeo_dgetr(nu, ny, tmp_nux_ny, 0, 0, J_y_tilde, 0, 0);
                        // J_y_tilde[:nu, :ny] += tmp_nux_ny^T[nu:,:] * S_forw_stage[:, nx:]
                        blasfeo_dgemm_tn(ny, nu, nx, 1.0, tmp_nux_ny, nu, 0, S_forw_stage, 0, nx, 1.0, J_y_tilde, 0, 0,
                                        J_y_tilde, 0, 0);

                        // J_y_tilde (x part)
                        blasfeo_dgemm_tn(ny, nx, nx, 1.0, tmp_nux_ny, nu, 0, S_forw_stage, 0, 0,
                                    0.0, J_y_tilde, 0, nu, J_y_tilde, 0, nu);

                        // transpose
                        blasfeo_dgetr(ny, nx+nu, J_y_tilde, 0, 0, tmp_nux_ny, 0, 0);

                        // tmp_nux_ny2 = W_chol * J_y_tilde (ny * (nx+nu))
                        blasfeo_dtrmm_rlnn(nu+nx, ny, 1.0, mem->W_chol, 0, 0, tmp_nux_ny, 0, 0,
                                        tmp_nux_ny2, 0, 0);

                        // cost_grad += b * J_y_tilde^T * tmp_ny
                        blasfeo_dgemv_t(ny, nx+nu, b_vec[ii]/num_steps, J_y_tilde, 0, 0, tmp_ny, 0,
                                        1.0, cost_grad, 0, cost_grad, 0);
                    }
                    // cost_hess += b * tmp_nux_ny2 * tmp_nux_ny2^T
                    blasfeo_dsyrk_ln(nu+nx, ny, b_vec[ii]/num_steps, tmp_nux_ny2, 0, 0, tmp_nux_ny2, 0, 0,
                            1.0, cost_hess, 0, 0, cost_hess, 0, 0);
                    // cost function value
                    // NOTE: slack contribution and scaling done in cost module
                    mem->cost_fun[0] += b_vec[ii]/num_steps * a;
                }
            }

            // update forward sensitivity
            // NOTE(oj): dK_dxu_ss is actually -dK_dxu_ss, because alpha = -1.0
            // was not supported by blasfeos backsolve initially.
            for (int jj = 0; jj < ns; jj++)
                blasfeo_dgead(nx, nx + nu, -step * b_vec[jj], dK_dxu_ss, jj * nx, 0,
                                                     S_forw_ss, 0, 0);
        }  // end if sens_forw || sens_hess
        // Cost computation without sensitivities
        else if (opts->cost_computation && opts->cost_type == NONLINEAR_LS)
        {
            ext_fun_arg_t nls_y_fun_type_in[4];
            void *nls_y_fun_in[4];
            ext_fun_arg_t nls_y_fun_type_out[1];
            void *nls_y_fun_out[1];

            nls_y_fun_type_in[0] = BLASFEO_DVEC;
            nls_y_fun_in[0] = xt;
            nls_y_fun_type_in[1] = COLMAJ;
            nls_y_fun_in[1] = u;
            nls_y_fun_type_in[2] = BLASFEO_DVEC_ARGS;
            nls_y_fun_in[2] = &impl_ode_z_in;
            nls_y_fun_type_in[3] = COLMAJ;
            nls_y_fun_in[3] = &t_current;
            impl_ode_z_in.x = &K_traj[ss];

            nls_y_fun_type_out[0] = BLASFEO_DVEC;
            nls_y_fun_out[0] = nls_res;  // fun: ny
            for (int ii = 0; ii < ns; ii++)
            {
                impl_ode_z_in.xi = ns * nx + ii * nz;

                t_current = t0 + ss * step + opts->c_vec[ii] * step;
                // compute x at stage (xt)
                blasfeo_dveccp(nx, xn, 0, xt, 0);
                for (int jj = 0; jj < ns; jj++)
                {
                    a = A_mat[ii + ns * jj] * step;
                    // xt = xt + T_int * a[i,j]*K_j
                    blasfeo_daxpy(nx, a, K, jj * nx, xt, 0, xt, 0);
                }

                model->nls_y_fun->evaluate(model->nls_y_fun, nls_y_fun_type_in,
                                nls_y_fun_in, nls_y_fun_type_out, nls_y_fun_out);

                // nls_res = nls_res - y_ref
                blasfeo_daxpy(ny, -1.0, mem->y_ref, 0, nls_res, 0, nls_res, 0);

                // tmp_ny = W_chol * nls_res
                blasfeo_dtrmv_lnn(ny, mem->W_chol, 0, 0, nls_res, 0, tmp_ny, 0);

                // cost function value
                // NOTE: slack contribution and scaling done in cost module
                mem->cost_fun[0] += 0.5 * b_vec[ii]/num_steps * blasfeo_ddot(ny, tmp_ny, 0, tmp_ny, 0);
            }
        } // end NLS cost_computation without sens
        else if (opts->cost_computation && opts->cost_type == CONVEX_OVER_NONLINEAR)
        {
            /* specify input types and pointers for external cost function */
            ext_fun_arg_t ext_fun_type_in[5];
            void *ext_fun_in[5];
            ext_fun_arg_t ext_fun_type_out[1];
            void *ext_fun_out[1];

            // INPUT
            ext_fun_type_in[0] = BLASFEO_DVEC;
            ext_fun_in[0] = xt;
            ext_fun_type_in[1] = COLMAJ;
            ext_fun_in[1] = u;
            ext_fun_type_in[2] = BLASFEO_DVEC;
            ext_fun_in[2] = &impl_ode_z_in;
            impl_ode_z_in.x = &K_traj[ss];
            ext_fun_type_in[3] = BLASFEO_DVEC;
            ext_fun_in[3] = mem->y_ref;
            ext_fun_type_in[4] = COLMAJ;
            ext_fun_in[4] = &t_current;

            // OUTPUT
            ext_fun_type_out[0] = COLMAJ;
            ext_fun_out[0] = &a;  // function: scalar

            for (int ii = 0; ii < ns; ii++)
            {
                impl_ode_z_in.xi = ns * nx + ii * nz;

                t_current = t0 + ss * step + opts->c_vec[ii] * step;
                // compute x at stage (xt)
                blasfeo_dveccp(nx, xn, 0, xt, 0);
                for (int jj = 0; jj < ns; jj++)
                {
                    a = A_mat[ii + ns * jj] * step;
                    // xt = xt + T_int * a[i,j]*K_j
                    blasfeo_daxpy(nx, a, K, jj * nx, xt, 0, xt, 0);
                }

                model->conl_cost_fun->evaluate(model->conl_cost_fun, ext_fun_type_in, ext_fun_in,
                                   ext_fun_type_out, ext_fun_out);

                // cost function value
                // NOTE: slack contribution and scaling done in cost module
                mem->cost_fun[0] += b_vec[ii]/num_steps * a;
            }
        } // end NLS cost_computation without sens


        // obtain x(n+1)
        for (int ii = 0; ii < ns; ii++){
            // xn += b_i * k_i
            blasfeo_daxpy(nx, step * b_vec[ii], K, ii * nx, xn, 0, xn, 0);
        }

        // algebraic variables output and corresponding sensitivity propagation
        if (ss == 0)
        {
            // generate z output
            if ((opts->output_z || opts->sens_algebraic) && nz > 0)
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

            if (opts->sens_algebraic && nz > 0 && !opts->exact_z_output)  // generate S_algebraic
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
                        // printf("\ndz[ii=%d]_dxu[jj=%d] = %e\n", ii, jj, interpolated_value);
                        // blasfeo_pack_dvec(1, &interpolated_value, 1, xtdot, ii);
                    }
                }
            }  // end if sens_algebraic


            if (opts->exact_z_output)
            {
                // printf("\nz output before newton:\n");
                // d_print_mat(1, nz, out->zn, 1);

                // set input for impl_ode
                impl_ode_type_in[0] = COLMAJ;
                impl_ode_type_in[1] = BLASFEO_DVEC;
                impl_ode_type_in[3] = COLMAJ;
                impl_ode_in[0] = in->x;  // 1st input is always xn
                impl_ode_in[1] = xtdot;
                impl_ode_in[3] = &out->zn[0];

                impl_ode_res_out.xi = 0;

                if (opts->output_z || opts->sens_algebraic)
                {
                    // initial guess for xdot0
                    for (int ii = 0; ii < nx; ii++)
                    {
                        double interpolated_value;
                        for (int jj = 0; jj < ns; jj++)
                        {
                            Z_work[jj] = blasfeo_dvecex1(K, nx * jj + ii);
                                // copy values of k_ii in first step, into Z_work
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
                            t_current = t0;
                            // eval jacobians at interpolated values
                            acados_tic(&timer_ad);
                            model->impl_ode_fun_jac_x_xdot_z->evaluate(
                                model->impl_ode_fun_jac_x_xdot_z, impl_ode_type_in, impl_ode_in,
                                impl_ode_fun_jac_x_xdot_z_type_out, impl_ode_fun_jac_x_xdot_z_out);
                            timing_ad += acados_toc(&timer_ad);

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
                // printf("\nz output after newton:\n");
                // d_print_mat(1, nz, out->zn, 1);

                if (opts->sens_algebraic)
                {
                    /* implicit function theorem to get S_alg */
                    // eval jacobians at interpolated values
                    acados_tic(&timer_ad);
                    model->impl_ode_jac_x_xdot_u_z->evaluate(
                            model->impl_ode_jac_x_xdot_u_z, impl_ode_type_in, impl_ode_in,
                            impl_ode_jac_x_xdot_u_z_type_out, impl_ode_jac_x_xdot_u_z_out);
                    timing_ad += acados_toc(&timer_ad);

                    // set up df_dxdotz
                    blasfeo_dgecp(nx + nz, nx, df_dxdot, 0, 0, df_dxdotz, 0, 0);
                    blasfeo_dgecp(nx + nz, nz, df_dz,    0, 0, df_dxdotz, 0, nx);
                    // set up right hand side dk0_dxu
                    blasfeo_dgecp(nx + nz, nx, df_dx, 0, 0, dk0_dxu, 0, 0);
                    blasfeo_dgecp(nx + nz, nu, df_du, 0, 0, dk0_dxu, 0, nx);

                    // solve linear system
                    acados_tic(&timer_la);
                    blasfeo_dgetrf_rp(nx + nz, nx + nz, df_dxdotz, 0, 0, df_dxdotz, 0, 0,
                                                                                ipiv_one_stage);
                    blasfeo_drowpe(nx + nz, ipiv_one_stage, dk0_dxu);
                    blasfeo_dtrsm_llnu(nx + nz, nx + nu, 1.0, df_dxdotz, 0, 0,
                                    dk0_dxu, 0, 0, dk0_dxu, 0, 0);
                    blasfeo_dtrsm_lunn(nx + nz, nx + nu, 1.0, df_dxdotz, 0, 0,
                                    dk0_dxu, 0, 0, dk0_dxu, 0, 0);
                    timing_la += acados_toc(&timer_la);

                    // solution has different sign
                    blasfeo_dgesc(nx + nz, nx + nu, -1.0, dk0_dxu, 0, 0);

                    // extract output
                    blasfeo_unpack_dmat(nz, nx + nu, dk0_dxu, nx, 0, S_algebraic, nz);
                } // if sens_algebraic
                // Reset impl_ode inputs
                impl_ode_type_in[0] = BLASFEO_DVEC;       // xt
                impl_ode_type_in[1] = BLASFEO_DVEC_ARGS;  // k_i
                impl_ode_type_in[3] = BLASFEO_DVEC_ARGS;  // z_i
                impl_ode_in[0] = xt;  // 1st input is always xt
                impl_ode_in[1] = &impl_ode_xdot_in;
                impl_ode_in[3] = &impl_ode_z_in;     // 4th input is part of Z[ss]
            } // if exact_z_output
        }  //  end if (ss == 0)
        if (ss == num_steps-1)
        {
            // store last xdot, z values for next initialization
            blasfeo_unpack_dvec(nx, K, (ns-1) * nx, mem->xdot, 1);
            blasfeo_unpack_dvec(nz, K, (ns-1) * nz + ns*nx, mem->z, 1);
        }
    }  // end step loop (ss)

    // extract results from forward sweep to output
    blasfeo_unpack_dvec(nx, xn, 0, x_out, 1);

    if  ( opts->sens_forw || opts->sens_hess )
        blasfeo_unpack_dmat(nx, nx + nu, S_forw_ss, 0, 0, S_forw_out, nx);

/*****************************************************************************
* Backward Sweep
*       - (adjoint sensitivities & hessian propagation)
*       - hessian via symmetric forward-backward sweep
*                    (see Algorithm 2 from Quirynen2016)
*******************************************************************************/
    if ( opts->sens_adj  || opts->sens_hess )
    {
        for (int ss = num_steps - 1; ss > -1; ss--)
        {
            if (opts->sens_hess){
                dK_dxu_ss = &dK_dxu[ss];
                dG_dK_ss = &dG_dK[ss];
                dG_dxu_ss = &dG_dxu[ss];
                ipiv_ss = &ipiv[ss*nK];
                S_forw_ss = &S_forw[ss];
                // lambdaK_ss = &lambdaK[ss];
                // lambda_ss_old = &lambda[ss+1];
                // lambda_ss = &lambda[ss];  // use something like this if lambda should be stored
            }
            else
            {
                dK_dxu_ss = dK_dxu;
                dG_dK_ss = dG_dK;
                dG_dxu_ss = dG_dxu;
                ipiv_ss = ipiv;
            }
            impl_ode_xdot_in.x = &K_traj[ss];              // use K values of step ss
            impl_ode_z_in.x = &K_traj[ss];                 // use Z values of step ss

        /* evaluate impl_ode_jac_x_xdot_u_z -- build dG_dxu_ss, dG_dK_ss
                                    & factorize dG_dK_ss  */
            if ( !opts->sens_hess )
            {
                blasfeo_dgese(nK, nK, 0.0, dG_dK_ss, 0, 0);   // initialize dG_dK_ss with zeros
                /* evaluate function at stage i, build corresponding blocks of dG_dxu, dG_dK_ss */
                for (int ii = 0; ii < ns; ii++)
                {
                    // set up input for impl_ode
                    impl_ode_xdot_in.xi = ii * nx;
                    // use k_i of K = (k_1,..., k_{ns},z_1,..., z_{ns})
                    impl_ode_z_in.xi    = ns * nx + ii * nz;
                    // use z_i of K = (k_1,..., k_{ns},z_1,..., z_{ns})
                    t_current = t0 + ss * step + opts->c_vec[ii] * step;

                    // build stage value
                    blasfeo_dveccp(nx, &xn_traj[ss], 0, xt, 0);
                    for (int jj = 0; jj < ns; jj++)
                    {
                        a = A_mat[ii + ns * jj] * step;
                        blasfeo_daxpy(nx, a, &K_traj[ss], jj * nx, xt, 0, xt, 0);
                    }
                    // set up input for impl_ode jacobians
                    acados_tic(&timer_ad);
                    model->impl_ode_jac_x_xdot_u_z->evaluate(
                        model->impl_ode_jac_x_xdot_u_z, impl_ode_type_in, impl_ode_in,
                        impl_ode_jac_x_xdot_u_z_type_out, impl_ode_jac_x_xdot_u_z_out);
                    timing_ad += acados_toc(&timer_ad);

                    // build dG_dxu_ss
                    blasfeo_dgecp(nx + nz, nx, df_dx, 0, 0, dG_dxu_ss, ii * (nx + nz), 0);
                    blasfeo_dgecp(nx + nz, nu, df_du, 0, 0, dG_dxu_ss, ii * (nx + nz), nx);

                    // build dG_dK_ss
                    for (int jj = 0; jj < ns; jj++)
                    {  // compute the block (ii,jj)th block of dG_dK_ss
                        a = A_mat[ii + ns * jj] * step;
                        blasfeo_dgead(nx + nz, nx, a, df_dx, 0, 0,
                                      dG_dK_ss, ii * (nx + nz), jj * nx);
                        if (jj == ii)
                        {
                            blasfeo_dgead(nx + nz, nx, 1.0, df_dxdot, 0, 0,
                                            dG_dK_ss, ii * (nx + nz), jj * nx);
                            blasfeo_dgead(nx + nz, nz, 1.0, df_dz,    0, 0,
                                            dG_dK_ss, ii * (nx + nz), (nx * ns) + jj * nz);
                        }
                    }  // end jj
                }  // end ii

                // factorize dG_dK_ss - already done in forw if hessian is active
                acados_tic(&timer_la);
                blasfeo_dgetrf_rp(nK, nK, dG_dK_ss, 0, 0, dG_dK_ss, 0, 0, ipiv_ss);
                timing_la += acados_toc(&timer_la);

            }  // end if( !opts->sens_hess )

            // update adjoint sensitivities: lambdaK
            // set up right hand side in vector lambdaK
            blasfeo_dvecse(nK, 0.0, lambdaK, 0);
            for (int jj = 0; jj < ns; jj++)
                blasfeo_dveccpsc(nx, -step * b_vec[jj], lambda, 0, lambdaK, jj * nx);
                // lambdaK_jj = -step b_jj * lambda_x

            //  obtain lambdaK by solving linear system lambdaK <- (dG_dK)^(-T) lambdaK;
            acados_tic(&timer_la);
            // dG_dK_ss - already factorized
            // solve linear system
            blasfeo_dtrsv_utn(nK, dG_dK_ss, 0, 0, lambdaK, 0, lambdaK, 0);
            blasfeo_dtrsv_ltu(nK, dG_dK_ss, 0, 0, lambdaK, 0, lambdaK, 0);
            blasfeo_dvecpei(nK, ipiv_ss, lambdaK, 0);
            timing_la += acados_toc(&timer_la);

            // update adjoint sensitivities lambda
            // lambda = 1 * lambda + 1 * dG_dxu_ss' * lambdaK
            blasfeo_dgemv_t(nK, nx + nu, 1.0, dG_dxu_ss, 0, 0, lambdaK, 0, 1.0, lambda,
                             0, lambda, 0);
            // Symmetric Hessian Propagation
            if ( opts->sens_hess )
            {
                impl_ode_hess_lambda_in.x = lambdaK;

                // evaluate second order derivatives and update Hessian
                // - HESSIAN PROPAGATION
                for (int ii = 0; ii < ns; ii++)
                {
                    blasfeo_dgecp(nx, nx+nu, S_forw_ss, 0, 0, dxkzu_dw0, 0, 0);
                    // printf("dxkzu_dw0 = (IRK, ss = %d) \n", ss);
                    // blasfeo_print_exp_dmat(2 * nx + nu + nz, nx + nu, dxkzu_dw0, 0, 0);
                    // build stage value, and dxii_dw0
                    blasfeo_dveccp(nx, &xn_traj[ss], 0, xt, 0);
                    for (int jj = 0; jj < ns; jj++)
                    {
                        a = A_mat[ii + ns * jj] * step;
                        blasfeo_daxpy(nx, a, &K_traj[ss], jj * nx, xt, 0, xt, 0);
                        // dxii_dw0 += a * dkjj_dxu
                        blasfeo_dgead(nx, nx + nu, -a, dK_dxu_ss, jj * nx, 0, dxkzu_dw0, 0, 0);
                    }
                    // dk_dw0
                    blasfeo_dgecpsc(nx, nx+nu, -1.0, dK_dxu_ss, ii*nx, 0, dxkzu_dw0, nx, 0);
                    // dz_dw0
                    blasfeo_dgecpsc(nz, nx+nu, -1.0, dK_dxu_ss, ns*nx+ii*nz, 0, dxkzu_dw0, 2*nx, 0);
                    // du_dw0
                    // TODO exploit the fact that this is [0, I] !!!
                    blasfeo_dgese(nu, nx+nu, 0.0, dxkzu_dw0, 2*nx+nz, 0);
                    blasfeo_ddiare(nu, 1.0, dxkzu_dw0, 2*nx+nz, nx);

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
                    impl_ode_xdot_in.xi = ii * nx;
                    // use k_i of K = (k_1,..., k_{ns},z_1,..., z_{ns})
                    impl_ode_z_in.xi    = ns * nx + ii * nz;
                    // use z_i of K = (k_1,..., k_{ns},z_1,..., z_{ns})
                    impl_ode_hess_lambda_in.xi = ii * (nx + nz);

                    t_current = t0 + ss * step + opts->c_vec[ii] * step;

                    // eval hessian function at stage ii
                    // printf("dxkzu_dw0 = (IRK, ss = %d) \n", ss);
                    // blasfeo_print_exp_dmat(2 * nx + nu + nz, nx + nu, dxkzu_dw0, 0, 0);
                    acados_tic(&timer_ad);

                    model->impl_ode_hess->evaluate(model->impl_ode_hess, impl_ode_hess_type_in,
                            impl_ode_hess_in, impl_ode_hess_type_out, impl_ode_hess_out);

                    timing_ad += acados_toc(&timer_ad);

#if 1
                    // exploit that du_dw0 is [0, I]
                    blasfeo_dgemm_nn(2*nx+nz+nu, nx+nu, 2*nx+nz, 1.0, f_hess, 0, 0, dxkzu_dw0, 0, 0, 0.0, tmp_dxkzu_dw0, 0, 0, tmp_dxkzu_dw0, 0, 0);
                    blasfeo_dgead(2*nx+nz+nu, nu, 1.0, f_hess, 0, 2*nx+nz, tmp_dxkzu_dw0, 0, nx);
                    blasfeo_dsyrk_ut(nx+nu, 2*nx+nz, 1.0, dxkzu_dw0, 0, 0, tmp_dxkzu_dw0, 0, 0, 1.0, Hess, 0, 0, Hess, 0, 0);
                    blasfeo_dgead(nu, nx+nu, 1.0, tmp_dxkzu_dw0, 2*nx+nz, 0, Hess, nx, 0);
#else
                    blasfeo_dgemm_nn(2*nx+nz+nu, nx+nu, 2*nx+nz+nu, 1.0, f_hess, 0, 0, dxkzu_dw0, 0, 0, 0.0, tmp_dxkzu_dw0, 0, 0, tmp_dxkzu_dw0, 0, 0);
                    blasfeo_dsyrk_ut(nx+nu, 2*nx+nz+nu, 1.0, dxkzu_dw0, 0, 0, tmp_dxkzu_dw0, 0, 0, 1.0, Hess, 0, 0, Hess, 0, 0);
#endif
                }  // end for ii
            }  // end if ( opts->sens_hess )
        }  // end for ss
    }  // end if ( opts->sens_adj  || opts->sens_hess )


    // extract output
    if  ( opts->sens_adj  || opts->sens_hess )
        blasfeo_unpack_dvec(nx + nu, lambda, 0, S_adj_out, 1);
    if  ( opts->sens_hess )
    {
        blasfeo_dtrtr_u(nu+nx, Hess, 0, 0, Hess, 0, 0);
        // printf("Hess = (IRK) \n");
        // blasfeo_print_exp_dmat(nx + nu, nx + nu, Hess, 0, 0);
        blasfeo_unpack_dmat(nx+nu, nx+nu, Hess, 0, 0, out->S_hess, nx + nu);
    }

    out->info->CPUtime = acados_toc(&timer);
    // note: this is the time for factorization and solving the linear systems
    out->info->LAtime = timing_la;
    out->info->ADtime = timing_ad;

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
    config->model_calculate_size = &sim_irk_model_calculate_size;
    config->model_assign = &sim_irk_model_assign;
    config->model_set = &sim_irk_model_set;
    config->dims_calculate_size = &sim_irk_dims_calculate_size;
    config->dims_assign = &sim_irk_dims_assign;
    config->dims_set = &sim_irk_dims_set;
    config->dims_get = &sim_irk_dims_get;
    return;
}
