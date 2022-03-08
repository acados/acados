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

// system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
// acados
#include "acados_c/sim_interface.h"

// mex
#include "mex.h"
#include "mex_macros.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    // sizeof(long long) == sizeof(void *) = 64 !!!
    int nx, nu, ii;

    long long *l_ptr;
    char *c_ptr;
    char fun_name[50] = "sim_create";
    char buffer [300]; // for error messages

    /* RHS */
    const mxArray *matlab_model = prhs[0];
    const mxArray *matlab_opts = prhs[1];

    /* LHS */
    // field names of output struct
    char *fieldnames[7];
    fieldnames[0] = (char*)mxMalloc(50);
    fieldnames[1] = (char*)mxMalloc(50);
    fieldnames[2] = (char*)mxMalloc(50);
    fieldnames[3] = (char*)mxMalloc(50);
    fieldnames[4] = (char*)mxMalloc(50);
    fieldnames[5] = (char*)mxMalloc(50);
    fieldnames[6] = (char*)mxMalloc(50);

    memcpy(fieldnames[0],"config",sizeof("config"));
    memcpy(fieldnames[1],"dims",sizeof("dims"));
    memcpy(fieldnames[2],"opts",sizeof("opts"));
    memcpy(fieldnames[3],"in",sizeof("in"));
    memcpy(fieldnames[4],"out",sizeof("out"));
    memcpy(fieldnames[5],"solver",sizeof("solver"));
    memcpy(fieldnames[6],"method",sizeof("method"));
    // create output struct
    plhs[0] = mxCreateStructMatrix(1, 1, 7, (const char **) fieldnames);

    mxFree( fieldnames[0] );
    mxFree( fieldnames[1] );
    mxFree( fieldnames[2] );
    mxFree( fieldnames[3] );
    mxFree( fieldnames[4] );
    mxFree( fieldnames[5] );
    mxFree( fieldnames[6] );


    /* plan & config */
    char *method = mxArrayToString( mxGetField( matlab_opts, 0, "method" ) );

    sim_solver_plan_t plan;

    if (!strcmp(method, "erk"))
    {
        plan.sim_solver = ERK;
    }
    else if (!strcmp(method, "irk"))
    {
        plan.sim_solver = IRK;
    }
    else if (!strcmp(method, "irk_gnsf"))
    {
        plan.sim_solver = GNSF;
    }
    else
    {
        MEX_FIELD_VALUE_NOT_SUPPORTED_SUGGEST(fun_name, "method", method, "erk, irk, irk_gnsf");
    }

    sim_config *config = sim_config_create(plan);

    /* dims */
    void *dims = sim_dims_create(config);

    if (mxGetField( matlab_model, 0, "dim_nx" )!=NULL)
    {
        nx = mxGetScalar( mxGetField( matlab_model, 0, "dim_nx" ) );
        sim_dims_set(config, dims, "nx", &nx);
    }
    else
    {
        MEX_MISSING_ARGUMENT(fun_name, "dim_nx");
    }
    
    if (mxGetField( matlab_model, 0, "dim_nu" )!=NULL)
    {
        nu = mxGetScalar( mxGetField( matlab_model, 0, "dim_nu" ) );
        sim_dims_set(config, dims, "nu", &nu);
    }
    else
    {
        nu = 0;
    }

    if (mxGetField( matlab_model, 0, "dim_nz" )!=NULL)
    {
        int nz = mxGetScalar( mxGetField( matlab_model, 0, "dim_nz" ) );
        sim_dims_set(config, dims, "nz", &nz);
    }

    if (!strcmp(method, "irk_gnsf"))
    {
        if (mxGetField( matlab_model, 0, "dim_gnsf_nx1" )!=NULL)
        {
            int gnsf_nx1 = mxGetScalar( mxGetField( matlab_model, 0, "dim_gnsf_nx1" ) );
            sim_dims_set(config, dims, "gnsf_nx1", &gnsf_nx1);
        }
        if (mxGetField( matlab_model, 0, "dim_gnsf_nz1" )!=NULL)
        {
            int gnsf_nz1 = mxGetScalar( mxGetField( matlab_model, 0, "dim_gnsf_nz1" ) );
            sim_dims_set(config, dims, "gnsf_nz1", &gnsf_nz1);
        }
        if (mxGetField( matlab_model, 0, "dim_gnsf_nuhat" )!=NULL)
        {
            int gnsf_nuhat = mxGetScalar( mxGetField( matlab_model, 0, "dim_gnsf_nuhat" ) );
            sim_dims_set(config, dims, "gnsf_nuhat", &gnsf_nuhat);
        }
        if (mxGetField( matlab_model, 0, "dim_gnsf_ny" )!=NULL)
        {
            int gnsf_ny = mxGetScalar( mxGetField( matlab_model, 0, "dim_gnsf_ny" ) );
            sim_dims_set(config, dims, "gnsf_ny", &gnsf_ny);
        }
        if (mxGetField( matlab_model, 0, "dim_gnsf_nout" )!=NULL)
        {
            int gnsf_nout = mxGetScalar( mxGetField( matlab_model, 0, "dim_gnsf_nout" ) );
            sim_dims_set(config, dims, "gnsf_nout", &gnsf_nout);
        }
    }


    /* opts */
    sim_opts *opts = sim_opts_create(config, dims);


    int num_stages = mxGetScalar( mxGetField( matlab_opts, 0, "num_stages" ) );
    int num_steps = mxGetScalar( mxGetField( matlab_opts, 0, "num_steps" ) );
    int newton_iter = mxGetScalar( mxGetField( matlab_opts, 0, "newton_iter" ) );

    // convert options to bool
    bool sens_forw = false;
    c_ptr = mxArrayToString( mxGetField( matlab_opts, 0, "sens_forw" ) );
    if (!strcmp(c_ptr, "true"))
    {
        sens_forw = true;
    }

    bool sens_adj = false;
    c_ptr = mxArrayToString( mxGetField( matlab_opts, 0, "sens_adj" ) );
    if (!strcmp(c_ptr, "true"))
    {
        sens_adj = true;
    }

    bool sens_hess = false;
    c_ptr = mxArrayToString( mxGetField( matlab_opts, 0, "sens_hess" ) );
    if (!strcmp(c_ptr, "true"))
    {
        sens_hess = true;
    }

    bool sens_algebraic = false;
    c_ptr = mxArrayToString( mxGetField( matlab_opts, 0, "sens_algebraic" ) );
    if (!strcmp(c_ptr, "true"))
    {
        sens_algebraic = true;
    }

    bool jac_reuse = false;
    c_ptr = mxArrayToString( mxGetField( matlab_opts, 0, "jac_reuse" ) );
    if (!strcmp(c_ptr, "true"))
    {
        jac_reuse = true;
    }

    char *collocation_type = mxArrayToString( mxGetField( matlab_opts, 0, "collocation_type" ) );
    sim_collocation_type collo_type;
    if (!strcmp(collocation_type, "gauss_legendre"))
    {
        collo_type = GAUSS_LEGENDRE;
    }
    else if (!strcmp(collocation_type, "gauss_radau_iia"))
    {
        collo_type = GAUSS_RADAU_IIA;
    }
    else
    {
        MEX_FIELD_VALUE_NOT_SUPPORTED_SUGGEST(fun_name, "collocation_type", collocation_type, "gauss_legendre, gauss_radau_iia");
    }

    sim_opts_set(config, opts, "num_stages", &num_stages);
    sim_opts_set(config, opts, "num_steps", &num_steps);
    sim_opts_set(config, opts, "newton_iter", &newton_iter);
    sim_opts_set(config, opts, "sens_forw", &sens_forw);
    sim_opts_set(config, opts, "sens_adj", &sens_adj);
    sim_opts_set(config, opts, "sens_hess", &sens_hess);
    sim_opts_set(config, opts, "sens_algebraic", &sens_algebraic);
    sim_opts_set(config, opts, "jac_reuse", &jac_reuse);
    sim_opts_set(config, opts, "collocation_type", &collo_type);



    /* in */
    sim_in *in = sim_in_create(config, dims);

    if (mxGetField( matlab_model, 0, "T" )!=NULL)
    {
        double T = mxGetScalar( mxGetField( matlab_model, 0, "T" ) );
        sim_in_set(config, dims, in, "T", &T);
    }
    else
    {
        MEX_MISSING_ARGUMENT(fun_name, "T");
    }

    if (mxGetField( matlab_model, 0, "x" )!=NULL)
    {
        double *x = mxGetPr( mxGetField( matlab_model, 0, "x" ) );
        sim_in_set(config, dims, in, "x", x);
    }
    if (mxGetField( matlab_model, 0, "u" )!=NULL)
    {
        double *u = mxGetPr( mxGetField( matlab_model, 0, "u" ) );
        sim_in_set(config, dims, in, "u", u);
    }
    if (mxGetField( matlab_model, 0, "seed_adj" )!=NULL)
    {
        double *seed_adj = mxGetPr( mxGetField( matlab_model, 0, "seed_adj" ) );
        sim_in_set(config, dims, in, "S_adj", seed_adj);
    }

    // unit forward seed
    if ((sens_forw==true) || (sens_hess==true))
    {
        double *Sx = calloc(nx*nx, sizeof(double));
        for(ii=0; ii<nx; ii++)
            Sx[ii*(nx+1)] = 1.0;
        double *Su = calloc(nx*nu, sizeof(double));
//        d_print_mat(nx, nx, Sx, nx);
//        d_print_mat(nx, nu, Su, nx);
        sim_in_set(config, dims, in, "Sx", Sx);
        sim_in_set(config, dims, in, "Su", Su);
        free(Sx);
        free(Su);
    }


    /* out */
    sim_out *out = sim_out_create(config, dims);

    /* solver */
    sim_solver *solver = sim_solver_create(config, dims, opts);

    /* populate output struct */
    // config
    mxArray *config_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(config_mat);
    l_ptr[0] = (long long) config;
    mxSetField(plhs[0], 0, "config", config_mat);

    // dims
    mxArray *dims_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(dims_mat);
    l_ptr[0] = (long long) dims;
    mxSetField(plhs[0], 0, "dims", dims_mat);

    // opts
    mxArray *opts_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(opts_mat);
    l_ptr[0] = (long long) opts;
    mxSetField(plhs[0], 0, "opts", opts_mat);

    // in
    mxArray *in_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(in_mat);
    l_ptr[0] = (long long) in;
    mxSetField(plhs[0], 0, "in", in_mat);

    // out
    mxArray *out_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(out_mat);
    l_ptr[0] = (long long) out;
    mxSetField(plhs[0], 0, "out", out_mat);

    // solver
    mxArray *solver_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(solver_mat);
    l_ptr[0] = (long long) solver;
    mxSetField(plhs[0], 0, "solver", solver_mat);

    // method
    mxArray *method_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(method_mat);
    l_ptr[0] = plan.sim_solver;
    mxSetField(plhs[0], 0, "method", method_mat);

    return;

}
