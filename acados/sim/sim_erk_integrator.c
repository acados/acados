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

// standard
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// acados
#include "acados/sim/sim_common.h"
#include "acados/sim/sim_erk_integrator.h"
#include "acados/utils/mem.h"

/************************************************
 * dims
 ************************************************/

int sim_erk_dims_calculate_size()
{
    int size = sizeof(sim_erk_dims);

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

int sim_erk_model_calculate_size(void *config, void *dims)
{
    int size = 0;

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
//    printf("\nsim_erk_model_set\n");
    erk_model *model = model_;

    if (!strcmp(field, "expl_ode_fun"))
    {
//    printf("\nsim_erk_model_set expl_ode_fun\n");
        model->expl_ode_fun = value;
    }
    else if (!strcmp(field, "expl_vde_for"))
    {
        model->expl_vde_for = value;
    }
    else if (!strcmp(field, "expl_vde_adj"))
    {
        model->expl_vde_adj = value;
    }
    else if (!strcmp(field, "expl_ode_hes"))
    {
        model->expl_ode_hes = value;
    }
    else
    {
        printf("\nerror: sim_erk_model_set: wrong field: %s\n", field);
		exit(1);
//        return ACADOS_FAILURE;
    }

    return ACADOS_SUCCESS;
}



/************************************************
 * opts
 ************************************************/

int sim_erk_opts_calculate_size(void *config_, void *dims)
{
    int ns_max = NS_MAX;

    int size = sizeof(sim_opts);

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
    opts->scheme = NULL;
    opts->jac_reuse = false;

    return (void *) opts;
}



int sim_erk_opts_set(void *config_, void *opts_, const char *field, void *value)
{
    sim_opts *opts = (sim_opts *) opts_;
    return sim_opts_set_(opts, field, value);
}



void sim_erk_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    sim_opts *opts = opts_;
    sim_erk_dims *dims = (sim_erk_dims *) dims_;

    opts->ns = 4;  // ERK 4
    int ns = opts->ns;

    assert((ns == 1 || ns == 2 || ns == 4) && "only number of stages = {1,2,4} implemented!");

    // set tableau size
    opts->tableau_size = opts->ns;

    double *A = opts->A_mat;
    double *b = opts->b_vec;
    double *c = opts->c_vec;

    switch (ns)
    {
        case 1:
        {
            // A
            A[0 + ns * 0] = 0.0;
            // b
            b[0] = 1.0;
            // c
            c[0] = 0.0;
            break;
        }
        case 2:
        {
            // A
            A[0 + ns * 0] = 0.0;
            A[0 + ns * 1] = 0.0;
            A[1 + ns * 0] = 0.5;
            A[1 + ns * 1] = 0.0;
            // b
            b[0] = 0.0;
            b[1] = 1.0;
            // c
            c[0] = 0.0;
            c[1] = 0.5;
            break;
        }
        case 4:
        {
            // A
            A[0 + ns * 0] = 0.0;
            A[0 + ns * 1] = 0.0;
            A[0 + ns * 2] = 0.0;
            A[0 + ns * 3] = 0.0;
            A[1 + ns * 0] = 0.5;
            A[1 + ns * 1] = 0.0;
            A[1 + ns * 2] = 0.0;
            A[1 + ns * 3] = 0.0;
            A[2 + ns * 0] = 0.0;
            A[2 + ns * 1] = 0.5;
            A[2 + ns * 2] = 0.0;
            A[2 + ns * 3] = 0.0;
            A[3 + ns * 0] = 0.0;
            A[3 + ns * 1] = 0.0;
            A[3 + ns * 2] = 1.0;
            A[3 + ns * 3] = 0.0;
            // b
            b[0] = 1.0 / 6.0;
            b[1] = 1.0 / 3.0;
            b[2] = 1.0 / 3.0;
            b[3] = 1.0 / 6.0;
            // c
            c[0] = 0.0;
            c[1] = 0.5;
            c[2] = 0.5;
            c[3] = 1.0;
            break;
        }
        default:
        {
            // impossible
            assert((ns == 1 || ns == 2 || ns == 4) &&
                   "only number of stages = {1,2,4} implemented!");
        }
    }

    opts->num_steps = 1;
    opts->num_forw_sens = dims->nx + dims->nu;
    opts->sens_forw = true;
    opts->sens_adj = false;
    opts->sens_hess = false;

    opts->output_z = false;
    opts->sens_algebraic = false;
}



void sim_erk_opts_update(void *config_, void *dims, void *opts_)
{
    sim_opts *opts = opts_;

    int ns = opts->ns;

    opts->tableau_size = opts->ns;

    assert((ns == 1 || ns == 2 || ns == 4) && "only number of stages = {1,2,4} implemented!");

    assert(ns <= NS_MAX && "ns > NS_MAX!");

    // set tableau size
    opts->tableau_size = opts->ns;

    double *A = opts->A_mat;
    double *b = opts->b_vec;
    double *c = opts->c_vec;

    switch (ns)
    {
        case 1:
        {
            // A
            A[0 + ns * 0] = 0.0;
            // b
            b[0] = 1.0;
            // c
            c[0] = 0.0;
            break;
        }
        case 2:
        {
            // A
            A[0 + ns * 0] = 0.0;
            A[0 + ns * 1] = 0.0;
            A[1 + ns * 0] = 0.5;
            A[1 + ns * 1] = 0.0;
            // b
            b[0] = 0.0;
            b[1] = 1.0;
            // c
            c[0] = 0.0;
            c[1] = 0.5;
            break;
        }
        case 4:
        {
            // A
            A[0 + ns * 0] = 0.0;
            A[0 + ns * 1] = 0.0;
            A[0 + ns * 2] = 0.0;
            A[0 + ns * 3] = 0.0;
            A[1 + ns * 0] = 0.5;
            A[1 + ns * 1] = 0.0;
            A[1 + ns * 2] = 0.0;
            A[1 + ns * 3] = 0.0;
            A[2 + ns * 0] = 0.0;
            A[2 + ns * 1] = 0.5;
            A[2 + ns * 2] = 0.0;
            A[2 + ns * 3] = 0.0;
            A[3 + ns * 0] = 0.0;
            A[3 + ns * 1] = 0.0;
            A[3 + ns * 2] = 1.0;
            A[3 + ns * 3] = 0.0;
            // b
            b[0] = 1.0 / 6.0;
            b[1] = 1.0 / 3.0;
            b[2] = 1.0 / 3.0;
            b[3] = 1.0 / 6.0;
            // c
            c[0] = 0.0;
            c[1] = 0.5;
            c[2] = 0.5;
            c[3] = 1.0;
            break;
        }
        default:
        {
            // impossible
            assert((ns == 1 || ns == 2 || ns == 4) &&
                   "only number of stages = {1,2,4} implemented!");
        }
    }

    return;
}



/************************************************
 * memory
 ************************************************/

int sim_erk_memory_calculate_size(void *config, void *dims, void *opts_) { return 0; }
void *sim_erk_memory_assign(void *config, void *dims, void *opts_, void *raw_memory)
{
    return NULL;
}

/************************************************
 * workspace
 ************************************************/

int sim_erk_workspace_calculate_size(void *config_, void *dims_, void *opts_)
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

    int size = sizeof(sim_erk_workspace);

    size += (nX + nu) * sizeof(double);  // rhs_forw_in

    if (opts->sens_adj)
    {
        size += num_steps * ns * nX * sizeof(double);   // K_traj
        size += (num_steps + 1) * nX * sizeof(double);  // out_forw_traj
    }
    else
    {
        size += ns * nX * sizeof(double);  // K_traj
        size += nX * sizeof(double);       // out_forw_traj
    }

    if (opts->sens_hess && opts->sens_adj)
    {
        size += (nx + nX + nu) * sizeof(double);          // rhs_adj_in
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

static void *sim_erk_cast_workspace(void *config_, void *dims_, void *opts_, void *raw_memory)
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

    char *c_ptr = (char *) raw_memory;

    sim_erk_workspace *workspace = (sim_erk_workspace *) c_ptr;
    c_ptr += sizeof(sim_erk_workspace);

    align_char_to(8, &c_ptr);

    assign_and_advance_double(nX + nu, &workspace->rhs_forw_in, &c_ptr);

    if (opts->sens_adj)
    {
        assign_and_advance_double(ns * num_steps * nX, &workspace->K_traj, &c_ptr);
        assign_and_advance_double((num_steps + 1) * nX, &workspace->out_forw_traj, &c_ptr);
    }
    else
    {
        assign_and_advance_double(ns * nX, &workspace->K_traj, &c_ptr);
        assign_and_advance_double(nX, &workspace->out_forw_traj, &c_ptr);
    }

    if (opts->sens_hess && opts->sens_adj)
    {
        assign_and_advance_double(nx + nX + nu, &workspace->rhs_adj_in, &c_ptr);
        assign_and_advance_double(nx + nu + nhess, &workspace->out_adj_tmp, &c_ptr);
        assign_and_advance_double(ns * (nx + nu + nhess), &workspace->adj_traj, &c_ptr);
    }
    else if (opts->sens_adj)
    {
        assign_and_advance_double((nx * 2 + nu), &workspace->rhs_adj_in, &c_ptr);
        assign_and_advance_double(nx + nu, &workspace->out_adj_tmp, &c_ptr);
        assign_and_advance_double(ns * (nx + nu), &workspace->adj_traj, &c_ptr);
    }

    assert((char *) raw_memory + sim_erk_workspace_calculate_size(config_, dims, opts_) >= c_ptr);

    return (void *) workspace;
}

/************************************************
 * functions
 ************************************************/

int sim_erk_precompute(void *config_, sim_in *in, sim_out *out, void *opts_, void *mem_,
                       void *work_)
{
    return ACADOS_SUCCESS;
}

int sim_erk(void *config_, sim_in *in, sim_out *out, void *opts_, void *mem_, void *work_)
{
    sim_config *config = config_;
    sim_opts *opts = opts_;

    if ( opts->ns != opts->tableau_size )
    {
        printf("Error in sim_erk: the Butcher tableau size does not match ns");
        return ACADOS_FAILURE;
    }
    int ns = opts->ns;

    void *dims_ = in->dims;
    sim_erk_dims *dims = (sim_erk_dims *) dims_;

    sim_erk_workspace *workspace =
        (sim_erk_workspace *) sim_erk_cast_workspace(config, dims, opts, work_);

    int i, j, s, istep;
    double a = 0, b = 0;  // temp values of A_mat and b_vec
    int nx = dims->nx;
    int nu = dims->nu;
    int nz = dims->nz;

    // assert - only use supported features
    if (nz != 0)
    {
        printf("nz should be zero - DAEs are not supported by the ERK integrator");
        return ACADOS_FAILURE;
    }
    if (opts->output_z)
    {
        printf("opts->output_z should be false - DAEs are not supported for the ERK integrator");
        return ACADOS_FAILURE;
    }
    if (opts->sens_algebraic)
    {
        printf("opts->sens_algebraic should be false - DAEs are not supported for the ERK integrator");
        return ACADOS_FAILURE;
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

    double *K_traj = workspace->K_traj;
    double *forw_traj = workspace->out_forw_traj;
    double *rhs_forw_in = workspace->rhs_forw_in;

    double *adj_tmp = workspace->out_adj_tmp;
    double *adj_traj = workspace->adj_traj;
    double *rhs_adj_in = workspace->rhs_adj_in;

    double *xn = out->xn;
    double *S_forw_out = out->S_forw;
    double *S_adj_out = out->S_adj;
    double *S_hess_out = out->S_hess;

    ext_fun_arg_t ext_fun_type_in[5];
    void *ext_fun_in[5];  // XXX large enough ?
    ext_fun_arg_t ext_fun_type_out[5];
    void *ext_fun_out[5];  // XXX large enough ?

    erk_model *model = in->model;

    acados_timer timer, timer_ad;
    double timing_ad = 0.0;

    acados_tic(&timer);

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
        if (opts->sens_adj)
        {
            K_traj = workspace->K_traj + istep * ns * nX;
            forw_traj = workspace->out_forw_traj + (istep + 1) * nX;
            for (i = 0; i < nX; i++) forw_traj[i] = forw_traj[i - nX];
        }

        for (s = 0; s < ns; s++)
        {
            for (i = 0; i < nX; i++) rhs_forw_in[i] = forw_traj[i];
            for (j = 0; j < s; j++)
            {
                a = A_mat[j * ns + s];
                if (a != 0)
                {
                    a *= step;
                    for (i = 0; i < nX; i++) rhs_forw_in[i] += a * K_traj[j * nX + i];
                }
            }

            acados_tic(&timer_ad);
            if (opts->sens_forw)
            {  // simulation + forward sensitivities
                ext_fun_type_in[0] = COLMAJ;
                ext_fun_in[0] = rhs_forw_in + 0;  // x: nx
                ext_fun_type_in[1] = COLMAJ;
                ext_fun_in[1] = rhs_forw_in + nx;  // Sx: nx*nx
                ext_fun_type_in[2] = COLMAJ;
                ext_fun_in[2] = rhs_forw_in + nx + nx * nx;  // Su: nx*nu
                ext_fun_type_in[3] = COLMAJ;
                ext_fun_in[3] = rhs_forw_in + nx + nx * nx + nx * nu;  // u: nu

                ext_fun_type_out[0] = COLMAJ;
                ext_fun_out[0] = K_traj + s * nX + 0;  // fun: nx
                ext_fun_type_out[1] = COLMAJ;
                ext_fun_out[1] = K_traj + s * nX + nx;  // Sx: nx*nx
                ext_fun_type_out[2] = COLMAJ;
                ext_fun_out[2] = K_traj + s * nX + nx + nx * nx;  // Su: nx*nu

                // forward VDE evaluation
                model->expl_vde_for->evaluate(model->expl_vde_for, ext_fun_type_in, ext_fun_in,
                                              ext_fun_type_out, ext_fun_out);
            }
            else
            {  // simulation only
                ext_fun_type_in[0] = COLMAJ;
                ext_fun_in[0] = rhs_forw_in + 0;  // x: nx
                ext_fun_type_in[1] = COLMAJ;
                ext_fun_in[1] = rhs_forw_in + nx;  // u: nu

                ext_fun_type_out[0] = COLMAJ;
                ext_fun_out[0] = K_traj + s * nX + 0;  // fun: nx

                model->expl_ode_fun->evaluate(model->expl_ode_fun, ext_fun_type_in, ext_fun_in,
                                              ext_fun_type_out, ext_fun_out);  // ODE evaluation
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
    if (opts->sens_adj || opts->sens_hess)
    {
        // initialize integrator variables
        for (i = 0; i < nx; i++) adj_tmp[i] = S_adj_in[i];
        for (i = 0; i < nu; i++) adj_tmp[nx + i] = 0.0;

        int nForw = nx;
        int nAdj = nx + nu;
        if (opts->sens_hess)
        {
            nForw = nX;
            nAdj = nx + nu + nhess;
            for (i = 0; i < nhess; i++) adj_tmp[nx + nu + i] = 0.0;
        }

        //        printf("\nnFOrw=%d nAdj=%d\n", nForw, nAdj);

        for (i = 0; i < nu; i++) rhs_adj_in[nForw + nx + i] = u[i];

        for (istep = num_steps - 1; istep > -1; istep--)
        {
            K_traj = workspace->K_traj + istep * ns * nX;
            forw_traj = workspace->out_forw_traj + (istep + 1) * nX;
            for (s = ns - 1; s > -1; s--)
            {
                // forward variables:
                for (i = 0; i < nForw; i++) rhs_adj_in[i] = forw_traj[i];  // extract x trajectory
                for (j = 0; j < s; j++)
                {
                    a = A_mat[j * ns + s];
                    if (a != 0)
                    {
                        a *= step;
                        for (i = 0; i < nForw; i++) rhs_adj_in[i] += a * K_traj[j * nX + i];
                    }  // plus k traj
                }
                // adjoint variables:
                b = step * b_vec[s];
                for (i = 0; i < nx; i++) rhs_adj_in[nForw + i] = b * adj_tmp[i];
                for (j = s + 1; j < ns; j++)
                {
                    a = A_mat[s * ns + j];
                    if (a != 0)
                    {
                        a *= step;
                        for (i = 0; i < nx; i++)
                            rhs_adj_in[nForw + i] += a * adj_traj[j * nAdj + i];
                    }
                }
                // TODO(oj): fix this whole file or write from scratch, not really readable :/
                acados_tic(&timer_ad);
                if (!opts->sens_hess)
                {
                    ext_fun_type_in[0] = COLMAJ;
                    ext_fun_in[0] = rhs_adj_in + 0;  // x: nx
                    ext_fun_type_in[1] = COLMAJ;
                    ext_fun_in[1] = rhs_adj_in + nx;  // Sx: nx*nx
                    ext_fun_type_in[2] = COLMAJ;
                    ext_fun_in[2] = rhs_adj_in + nx + nx * nx;  // Su: nx*nu
                    ext_fun_type_in[3] = COLMAJ;
                    ext_fun_in[3] = rhs_adj_in + nx + nx * nx + nx * nu;  // lam: nx
                    ext_fun_type_in[4] = COLMAJ;
                    ext_fun_in[4] = rhs_adj_in + nx + nx * nx + nx * nu + nx;  // u: nu

                    ext_fun_type_out[0] = COLMAJ;
                    ext_fun_out[0] = adj_traj + s * nAdj + 0;  // adj: nx+nu
                    ext_fun_type_out[1] = COLMAJ;
                    ext_fun_out[1] = adj_traj + s * nAdj + nx + nu;  // hess: (nx+nu)*(nx+nu)

                    model->expl_vde_adj->evaluate(model->expl_vde_adj, ext_fun_type_in, ext_fun_in,
                                                  ext_fun_type_out,
                                                  ext_fun_out);  // adjoint VDE evaluation
                }
                else
                {
                    ext_fun_type_in[0] = COLMAJ;
                    ext_fun_in[0] = rhs_adj_in + 0;  // x: nx
                    ext_fun_type_in[1] = COLMAJ;
                    ext_fun_in[1] = rhs_adj_in + nx;  // Sx: nx*nx
                    ext_fun_type_in[2] = COLMAJ;
                    ext_fun_in[2] = rhs_adj_in + nx + nx * nx;  // Su: nx*nu
                    ext_fun_type_in[3] = COLMAJ;
                    ext_fun_in[3] = rhs_adj_in + nx + nx * nx + nx * nu;  // lam: nx
                    ext_fun_type_in[4] = COLMAJ;
                    ext_fun_in[4] = rhs_adj_in + nx + nx * nx + nx * nu + nx;  // u: nu

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
                for (i = 0; i < nAdj; i++) adj_tmp[i] += adj_traj[s * nAdj + i];  // ERK step
        }

        // store adjoint sensitivities
        for (i = 0; i < nx + nu; i++) S_adj_out[i] = adj_tmp[i];
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
                    S_hess_out[j + (nf) *i] = adj_tmp[nx + nu + count_upper];
                    // copy to upper part
                    count_upper++;
                }
            }
        }
    }

    // store timings
    out->info->CPUtime = acados_toc(&timer);
    out->info->LAtime = 0.0;
    out->info->ADtime = timing_ad;

    // return
    return 0;  // success
}

void sim_erk_config_initialize_default(void *config_)
{
    sim_config *config = config_;

    config->opts_calculate_size = &sim_erk_opts_calculate_size;
    config->opts_assign = &sim_erk_opts_assign;
    config->opts_initialize_default = &sim_erk_opts_initialize_default;
    config->opts_update = &sim_erk_opts_update;
    config->opts_set = &sim_erk_opts_set;
    config->memory_calculate_size = &sim_erk_memory_calculate_size;
    config->memory_assign = &sim_erk_memory_assign;
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
