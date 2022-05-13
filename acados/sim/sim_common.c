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


// standard
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// acados
#include "acados/sim/sim_common.h"
#include "acados/utils/mem.h"
#include "acados/utils/print.h"



/************************************************
 * config
 ************************************************/

acados_size_t sim_config_calculate_size()
{
    acados_size_t size = 0;

    size += sizeof(sim_config);

    return size;
}



sim_config *sim_config_assign(void *raw_memory)
{
    char *c_ptr = raw_memory;

    sim_config *config = (sim_config *) c_ptr;
    c_ptr += sizeof(sim_config);

    assert((char *) raw_memory + sim_config_calculate_size() >= c_ptr);

    return config;
}



/************************************************
 * in
 ************************************************/

acados_size_t sim_in_calculate_size(void *config_, void *dims)
{
    sim_config *config = config_;

    acados_size_t size = sizeof(sim_in);

    int nx, nu, nz;

    config->dims_get(config_, dims, "nx", &nx);
    config->dims_get(config_, dims, "nu", &nu);
    config->dims_get(config_, dims, "nz", &nz);

    size += nx * sizeof(double);              // x
    size += nu * sizeof(double);              // u
    size += nx * (nx + nu) * sizeof(double);  // S_forw (max dimension)
    size += (nx + nu) * sizeof(double);       // S_adj

    size += config->model_calculate_size(config, dims);

    make_int_multiple_of(8, &size);
    size += 2 * 8;

    return size;
}



sim_in *sim_in_assign(void *config_, void *dims, void *raw_memory)
{
    sim_config *config = config_;

    char *c_ptr = (char *) raw_memory;

    // assign structure itself
    sim_in *in = (sim_in *) c_ptr;
    c_ptr += sizeof(sim_in);

    // align
    align_char_to(8, &c_ptr);

    // assign substructures
    in->model = config->model_assign(config, dims, c_ptr);
    c_ptr += config->model_calculate_size(config, dims);

    // set pointers and dimensions, defaults
    in->dims = dims;
    int nx, nu, nz;
    config->dims_get(config_, dims, "nx", &nx);
    config->dims_get(config_, dims, "nu", &nu);
    config->dims_get(config_, dims, "nz", &nz);
    int NF = nx + nu;
    in->identity_seed = false;

    // align
    align_char_to(8, &c_ptr);

    // assign doubles
    assign_and_advance_double(nx, &in->x, &c_ptr);
    assign_and_advance_double(nu, &in->u, &c_ptr);
    assign_and_advance_double(nx * NF, &in->S_forw, &c_ptr);
    assign_and_advance_double(NF, &in->S_adj, &c_ptr);

    assert((char *) raw_memory + sim_in_calculate_size(config_, dims) >= c_ptr);

    return in;
}



int sim_in_set_(void *config_, void *dims_, sim_in *in, const char *field, void *value)
{
    sim_config *config = config_;

    int status = ACADOS_SUCCESS;

    if (!strcmp(field, "T"))
    {
        double *T = value;
        in->T = T[0];
    }
    else if (!strcmp(field, "x"))
    {
        int nx;
        config->dims_get(config_, dims_, "nx", &nx);
        double *x = value;
        for (int ii=0; ii < nx; ii++)
            in->x[ii] = x[ii];
    }
    else if (!strcmp(field, "u"))
    {
        int nu;
        config->dims_get(config_, dims_, "nu", &nu);
        double *u = value;
        for (int ii=0; ii < nu; ii++)
            in->u[ii] = u[ii];
    }
    else if (!strcmp(field, "Sx"))
    {
        // note: this assumes nf = nu+nx !!!
        int nx;
        config->dims_get(config_, dims_, "nx", &nx);
        double *Sx = value;
        for (int ii=0; ii < nx*nx; ii++)
            in->S_forw[ii] = Sx[ii];
    }
    else if (!strcmp(field, "Su"))
    {
        // note: this assumes nf = nu+nx !!!
        int nx, nu;
        config->dims_get(config_, dims_, "nx", &nx);
        config->dims_get(config_, dims_, "nu", &nu);
        double *Su = value;
        for (int ii=0; ii < nx*nu; ii++)
            in->S_forw[nx*nx+ii] = Su[ii];
    }
    else if (!strcmp(field, "S_forw"))
    {
        // NOTE: this assumes nf = nu+nx !!!
        int nx, nu;
        config->dims_get(config_, dims_, "nx", &nx);
        config->dims_get(config_, dims_, "nu", &nu);
        double *S_forw = value;
        for (int ii=0; ii < nx*(nu+nx); ii++)
            in->S_forw[ii] = S_forw[ii];
    }
    else if (!strcmp(field, "S_adj"))
    {
        // NOTE: this assumes nf = nu+nx !!!
        int nx, nu;
        config->dims_get(config_, dims_, "nx", &nx);
        config->dims_get(config_, dims_, "nu", &nu);
        double *S_adj = value;
        for (int ii=0; ii < nx+nu; ii++)
            in->S_adj[ii] = S_adj[ii];
    }
    else if (!strcmp(field, "seed_adj"))
    {
        // NOTE: this assumes nf = nu+nx !!!
        // NOTE: this correctly initialized the u-part to 0, unless the above S_adj which copies it from outside
        int nx, nu;
        config->dims_get(config_, dims_, "nx", &nx);
        config->dims_get(config_, dims_, "nu", &nu);
        double *seed_adj = value;
        for (int ii=0; ii < nx; ii++)
            in->S_adj[ii] = seed_adj[ii];
        for (int ii=0; ii < nu; ii++)
            in->S_adj[nx+ii] = 0;
    }
    else
    {
        status = config->model_set(in->model, field, value);
    }

    return status;
}



/************************************************
 * out
 ************************************************/

acados_size_t sim_out_calculate_size(void *config_, void *dims)
{
    sim_config *config = config_;

    acados_size_t size = sizeof(sim_out);

    int nx, nu, nz;
    config->dims_get(config_, dims, "nx", &nx);
    config->dims_get(config_, dims, "nu", &nu);
    config->dims_get(config_, dims, "nz", &nz);

    int NF = nx + nu;
    size += sizeof(sim_info);

    size += nx * sizeof(double);                // xn
    size += nx * NF * sizeof(double);           // S_forw
    size += (nx + nu) * sizeof(double);         // S_adj
    size += (NF * NF) * sizeof(double);         // S_hess

    size += nz * sizeof(double);                //  zn
    size += nz * NF * sizeof(double);           // S_algebraic

    size += NF * sizeof(double);                // grad

    make_int_multiple_of(8, &size);
    size += 1 * 8;

    return size;
}



sim_out *sim_out_assign(void *config_, void *dims, void *raw_memory)
{
    sim_config *config = config_;

    char *c_ptr = (char *) raw_memory;

    int nx, nu, nz;
    config->dims_get(config_, dims, "nx", &nx);
    config->dims_get(config_, dims, "nu", &nu);
    config->dims_get(config_, dims, "nz", &nz);

    int NF = nx + nu;

    sim_out *out = (sim_out *) c_ptr;
    c_ptr += sizeof(sim_out);

    out->info = (sim_info *) c_ptr;
    c_ptr += sizeof(sim_info);

    align_char_to(8, &c_ptr);

    assign_and_advance_double(nx, &out->xn, &c_ptr);
    assign_and_advance_double(nx * NF, &out->S_forw, &c_ptr);
    assign_and_advance_double(nx + nu, &out->S_adj, &c_ptr);
    assign_and_advance_double(NF * NF, &out->S_hess, &c_ptr);
    assign_and_advance_double(NF, &out->grad, &c_ptr);

    assign_and_advance_double(nz, &out->zn, &c_ptr);
    assign_and_advance_double(nz * NF, &out->S_algebraic, &c_ptr);

    assert((char *) raw_memory + sim_out_calculate_size(config_, dims) >= c_ptr);

    return out;
}



int sim_out_get_(void *config_, void *dims_, sim_out *out, const char *field, void *value)
{
    sim_config *config = config_;

    int status = ACADOS_SUCCESS;

    if (!strcmp(field, "xn") || !strcmp(field, "x"))
    {
        int nx;
        config->dims_get(config_, dims_, "nx", &nx);
        double *xn = value;
        for (int ii=0; ii < nx; ii++)
            xn[ii] = out->xn[ii];
    }
    else if (!strcmp(field, "zn") || !strcmp(field, "z"))
    {
        int nz;
        config->dims_get(config_, dims_, "nz", &nz);
        double *zn = value;
        for (int ii=0; ii < nz; ii++)
            zn[ii] = out->zn[ii];
    }
    else if (!strcmp(field, "S_forw"))
    {
        // note: this assumes nf = nu+nx !!!
        int nx, nu;
        config->dims_get(config_, dims_, "nx", &nx);
        config->dims_get(config_, dims_, "nu", &nu);
        int ii;
        double *S_forw = value;
        for (ii=0; ii < nx*(nu+nx); ii++)
            S_forw[ii] = out->S_forw[ii];
    }
    else if (!strcmp(field, "Sx"))
    {
        // note: this assumes nf = nu+nx !!!
        int nx;
        config->dims_get(config_, dims_, "nx", &nx);
        int ii;
        double *Sx = value;
        for (ii=0; ii < nx*nx; ii++)
            Sx[ii] = out->S_forw[ii];
    }
    else if (!strcmp(field, "Su"))
    {
        // note: this assumes nf = nu+nx !!!
        int nx, nu;
        config->dims_get(config_, dims_, "nx", &nx);
        config->dims_get(config_, dims_, "nu", &nu);
        int ii;
        double *Su = value;
        for (ii=0; ii < nx*nu; ii++)
            Su[ii] = out->S_forw[nx*nx+ii];
    }
    else if (!strcmp(field, "S_adj"))
    {
        // note: this assumes nf = nu+nx !!!
        int nx, nu;
        config->dims_get(config_, dims_, "nx", &nx);
        config->dims_get(config_, dims_, "nu", &nu);
        int ii;
        double *S_adj = value;
        for (ii=0; ii < nu+nx; ii++)
            S_adj[ii] = out->S_adj[ii];
    }
    else if (!strcmp(field, "S_hess"))
    {
        // note: this assumes nf = nu+nx !!!
        int nx, nu;
        config->dims_get(config_, dims_, "nx", &nx);
        config->dims_get(config_, dims_, "nu", &nu);
        int ii;
        double *S_hess = value;
        for (ii=0; ii < (nu+nx)*(nu+nx); ii++)
            S_hess[ii] = out->S_hess[ii];
    }
    else if (!strcmp(field, "S_algebraic"))
    {
        // note: this assumes nf = nu+nx !!!
        int nx, nu, nz;
        config->dims_get(config_, dims_, "nz", &nz);
        config->dims_get(config_, dims_, "nx", &nx);
        config->dims_get(config_, dims_, "nu", &nu);
        double *S_algebraic = value;
        for (int ii=0; ii < nz*(nu+nx); ii++)
            S_algebraic[ii] = out->S_algebraic[ii];
    }
    else if (!strcmp(field, "CPUtime") || !strcmp(field, "time_tot"))
    {
        double *time = value;
        *time = out->info->CPUtime;
    }
    else if (!strcmp(field, "ADtime") || !strcmp(field, "time_ad"))
    {
        double *time = value;
        *time = out->info->ADtime;
    }
    else if (!strcmp(field, "LAtime") || !strcmp(field, "time_la"))
    {
        double *time = value;
        *time = out->info->LAtime;
    }
    else
    {
        printf("sim_out_get_: field %s not supported \n", field);
        exit(1);
    }

    return status;
}



/************************************************
* sim_opts
************************************************/

void sim_opts_set_(sim_opts *opts, const char *field, void *value)
{

    if (!strcmp(field, "ns") ||!strcmp(field, "num_stages"))
    {
        int *ns = (int *) value;
        opts->ns = *ns;
    }
    else if (!strcmp(field, "num_steps"))
    {
        int *num_steps = (int *) value;
        opts->num_steps = *num_steps;
    }
    else if (!strcmp(field, "newton_iter"))
    {
        int *newton_iter = (int *) value;
        opts->newton_iter = *newton_iter;
    }
    else if (!strcmp(field, "jac_reuse"))
    {
        bool *jac_reuse = (bool *) value;
        opts->jac_reuse = *jac_reuse;
    }
    else if (!strcmp(field, "sens_forw"))
    {
        bool *sens_forw = (bool *) value;
        opts->sens_forw = *sens_forw;
    }
    else if (!strcmp(field, "sens_adj"))
    {
        bool *sens_adj = (bool *) value;
        opts->sens_adj = *sens_adj;
    }
    else if (!strcmp(field, "sens_hess"))
    {
        bool *sens_hess = (bool *) value;
        opts->sens_hess = *sens_hess;
    }
    else if (!strcmp(field, "output_z"))
    {
        bool *output_z = (bool *) value;
        opts->output_z = *output_z;
    }
    else if (!strcmp(field, "exact_z_output"))
    {
        bool *exact_z_output = (bool *) value;
        opts->exact_z_output = *exact_z_output;
    }
    else if (!strcmp(field, "sens_algebraic"))
    {
        bool *sens_algebraic = (bool *) value;
        opts->sens_algebraic = *sens_algebraic;
    }
    else if (!strcmp(field, "collocation_type"))
    {
        sim_collocation_type *collocation_type = (sim_collocation_type *) value;
        opts->collocation_type = *collocation_type;
    }
    else
    {
        printf("\nerror: field %s not available in sim_opts_set\n", field);
        exit(1);
    }

    return;
}


void sim_opts_get_(sim_config *config, sim_opts *opts, const char *field, void *value)
{

    if (!strcmp(field, "sens_forw"))
    {
        bool *sens_forw = value;
        *sens_forw = opts->sens_forw;
    }
    else if (!strcmp(field, "sens_adj"))
    {
        bool *sens_adj = value;
        *sens_adj = opts->sens_adj;
    }
    else if (!strcmp(field, "sens_hess"))
    {
        bool *sens_hess = value;
        *sens_hess = opts->sens_hess;
    }
    else
    {
        printf("sim_opts_get: field %s not supported \n", field);
        exit(1);
    }

    return;
}
