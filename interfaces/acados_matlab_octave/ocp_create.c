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
#include "acados_c/ocp_nlp_interface.h"
// mex
#include "mex.h"
#include "mex_macros.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//mexPrintf("\nstart ocp create\n");

    long long *l_ptr;
    int *i_ptr;
    int *tmp_idx;
    double *d_ptr;
    char *c_ptr;

    int N;
    char fun_name[50] = "ocp_create";
    char buffer[500]; // for error messages
    char matlab_field_name[100];

    double acados_inf = 1e8;

    /* RHS */
    const mxArray *matlab_model = prhs[0];
    const mxArray *matlab_opts = prhs[1];

    /* LHS */
    // field names of output struct
    char *fieldnames[8];
    fieldnames[0] = (char*) mxMalloc(50);
    fieldnames[1] = (char*) mxMalloc(50);
    fieldnames[2] = (char*) mxMalloc(50);
    fieldnames[3] = (char*) mxMalloc(50);
    fieldnames[4] = (char*) mxMalloc(50);
    fieldnames[5] = (char*) mxMalloc(50);
    fieldnames[6] = (char*) mxMalloc(50);
    fieldnames[7] = (char*) mxMalloc(50);

    memcpy(fieldnames[0],"config",sizeof("config"));
    memcpy(fieldnames[1],"dims",sizeof("dims"));
    memcpy(fieldnames[2],"opts",sizeof("opts"));
    memcpy(fieldnames[3],"in",sizeof("in"));
    memcpy(fieldnames[4],"out",sizeof("out"));
    memcpy(fieldnames[5],"solver",sizeof("solver"));
    memcpy(fieldnames[6],"sens_out",sizeof("sens_out"));
    memcpy(fieldnames[7],"plan",sizeof("plan"));

    // create output struct
    plhs[0] = mxCreateStructMatrix(1, 1, 8, (const char **) fieldnames);

    mxFree( fieldnames[0] );
    mxFree( fieldnames[1] );
    mxFree( fieldnames[2] );
    mxFree( fieldnames[3] );
    mxFree( fieldnames[4] );
    mxFree( fieldnames[5] );
    mxFree( fieldnames[6] );
    mxFree( fieldnames[7] );


    /* plan & config */
    // param_scheme_N
    if (mxGetField( matlab_opts, 0, "param_scheme_N" )!=NULL)
    {
        N = mxGetScalar( mxGetField( matlab_opts, 0, "param_scheme_N" ) );
    }
    else
    {
       MEX_MISSING_ARGUMENT(fun_name, "param_scheme_N");
    }
    ocp_nlp_plan *plan = ocp_nlp_plan_create(N);


    // nlp solver
    char *nlp_solver = mxArrayToString( mxGetField( matlab_opts, 0, "nlp_solver" ) );
    if (!strcmp(nlp_solver, "sqp"))
    {
        plan->nlp_solver = SQP;
    }
    else if (!strcmp(nlp_solver, "sqp_rti"))
    {
        plan->nlp_solver = SQP_RTI;
    }
    else
    {
        MEX_FIELD_VALUE_NOT_SUPPORTED_SUGGEST(fun_name, "nlp_solver", nlp_solver, "sqp, sqp_rti");
    }

    // cost type
    char *cost_type;
    ocp_nlp_cost_t cost_type_enum;
    if (mxGetField( matlab_model, 0, "cost_type" )!=NULL)
    {
        cost_type = mxArrayToString( mxGetField( matlab_model, 0, "cost_type" ) );
        if (!strcmp(cost_type, "linear_ls"))
        {
            cost_type_enum = LINEAR_LS;
        }
        else if (!strcmp(cost_type, "nonlinear_ls"))
        {
            cost_type_enum = NONLINEAR_LS;
        }
        else if (!strcmp(cost_type, "ext_cost"))
        {
            cost_type_enum = EXTERNAL;
        }
        else
        {
            MEX_FIELD_VALUE_NOT_SUPPORTED_SUGGEST(fun_name, "cost_type", cost_type, "linear_ls, nonlinear_ls, ext_cost");
        }
    }
    else
    {
        MEX_MISSING_ARGUMENT(fun_name, "cost_type");
    }

    for (int ii=0; ii<N; ii++)
    {
        plan->nlp_cost[ii] = cost_type_enum;
    }

    // cost e_type
    char *cost_type_e;
    ocp_nlp_cost_t cost_type_e_enum;
    if (mxGetField( matlab_model, 0, "cost_type_e" )!=NULL)
    {
        cost_type_e = mxArrayToString( mxGetField( matlab_model, 0, "cost_type_e" ) );
        if (!strcmp(cost_type_e, "linear_ls"))
        {
            cost_type_e_enum = LINEAR_LS;
        }
        else if (!strcmp(cost_type_e, "nonlinear_ls"))
        {
            cost_type_e_enum = NONLINEAR_LS;
        }
        else if (!strcmp(cost_type_e, "ext_cost"))
        {
            cost_type_e_enum = EXTERNAL;
        }
        else
        {
            MEX_FIELD_VALUE_NOT_SUPPORTED_SUGGEST(fun_name, "cost_type_e", cost_type_e, "linear_ls, nonlinear_ls, ext_cost");
        }
    }
    else
    {
        MEX_MISSING_ARGUMENT(fun_name, "cost_type_e");
    }
    plan->nlp_cost[N] = cost_type_e_enum;


    // dynamics: dyn_type, sim_method
    char *dyn_type;
    ocp_nlp_dynamics_t dyn_type_enum;
    if (mxGetField( matlab_model, 0, "dyn_type" )!=NULL)
    {
        dyn_type = mxArrayToString( mxGetField( matlab_model, 0, "dyn_type" ) );
        if (!strcmp(dyn_type, "explicit"))
        {
            dyn_type_enum = CONTINUOUS_MODEL;
        }
        else if (!strcmp(dyn_type, "implicit"))
        {
            dyn_type_enum = CONTINUOUS_MODEL;
        }
        else if (!strcmp(dyn_type, "discrete"))
        {
            dyn_type_enum = DISCRETE_MODEL;
        }
        else
        {
            MEX_FIELD_VALUE_NOT_SUPPORTED_SUGGEST(fun_name, "dyn_type", dyn_type, "explicit, implicit, discrete");
        }
    }
    else
    {
        MEX_MISSING_ARGUMENT(fun_name, "dyn_type");
    }
    for (int ii=0; ii<N; ii++)
    {
        plan->nlp_dynamics[ii] = dyn_type_enum;
    }

    char *sim_method;
    sim_solver_t sim_method_enum = INVALID_SIM_SOLVER;

    if (strcmp(dyn_type, "discrete")!=0) // discrete dont need integrator
    {
        if (mxGetField( matlab_opts, 0, "sim_method" )!=NULL)
        {
            sim_method = mxArrayToString( mxGetField( matlab_opts, 0, "sim_method" ) );
        }
        else
        {
            MEX_MISSING_ARGUMENT(fun_name, "sim_method");
        }
    }
    else
    {
        sim_method_enum = INVALID_SIM_SOLVER; // discrete
    }

    if (!strcmp(dyn_type, "explicit"))
    {
        if (!strcmp(sim_method, "erk"))
        {
            sim_method_enum = ERK;
        }
        else
        {
            MEX_FIELD_NOT_SUPPORTED_GIVEN(fun_name, "sim_method", sim_method, "explicit dynamics", "erk");
        }
    }
    else if (!strcmp(dyn_type, "implicit"))
    {
        if (!strcmp(sim_method, "irk"))
        {
            sim_method_enum = IRK;
        }
        else if (!strcmp(sim_method, "irk_gnsf"))
        {
            sim_method_enum = GNSF;
        }
        else
        {
            MEX_FIELD_NOT_SUPPORTED_GIVEN(fun_name, "sim_method", sim_method, "implicit dynamics", "irk, irk_gnsf");
        }
    }
    for (int ii=0; ii<N; ii++)
    {
        plan->sim_solver_plan[ii].sim_solver = sim_method_enum;
    }


    // constraints
    char *constr_type;
    if (mxGetField( matlab_model, 0, "constr_type" )!=NULL)
    {
        constr_type = mxArrayToString( mxGetField( matlab_model, 0, "constr_type" ) );
    }
    else
    {
        MEX_MISSING_ARGUMENT(fun_name, "constr_type");
    }
    if (!strcmp(constr_type, "bgh"))
    {
        for (int ii=0; ii<N; ii++)
        {
            plan->nlp_constraints[ii] = BGH;
        }
    }
    else
    {
        MEX_FIELD_VALUE_NOT_SUPPORTED_SUGGEST(fun_name, "constr_type", constr_type, "bgh");
    }

    char *constr_type_e;
    if (mxGetField( matlab_model, 0, "constr_type_e" )!=NULL)
    {
        constr_type_e = mxArrayToString( mxGetField( matlab_model, 0, "constr_type_e" ) );
    }
    else
    {
        MEX_MISSING_ARGUMENT(fun_name, "constr_type_e");
    }
    if (!strcmp(constr_type_e, "bgh"))
    {
        plan->nlp_constraints[N] = BGH;
    }
    else
    {
        MEX_FIELD_VALUE_NOT_SUPPORTED_SUGGEST(fun_name, "constr_type_e", constr_type_e, "bgh");
    }

    // qp solver
    char *qp_solver;
    if (mxGetField( matlab_opts, 0, "qp_solver" )!=NULL)
    {
        qp_solver = mxArrayToString( mxGetField( matlab_opts, 0, "qp_solver" ) );
    }
    else
    {
        MEX_MISSING_ARGUMENT(fun_name, "qp_solver");
    }

    if (!strcmp(qp_solver, "partial_condensing_hpipm"))
    {
        plan->ocp_qp_solver_plan.qp_solver = PARTIAL_CONDENSING_HPIPM;
    }
    else if (!strcmp(qp_solver, "full_condensing_hpipm"))
    {
        plan->ocp_qp_solver_plan.qp_solver = FULL_CONDENSING_HPIPM;
    }
#if defined(ACADOS_WITH_QPOASES)
    else if (!strcmp(qp_solver, "full_condensing_qpoases"))
    {
        plan->ocp_qp_solver_plan.qp_solver = FULL_CONDENSING_QPOASES;
    }
#endif
#if defined(ACADOS_WITH_HPMPC)
    else if (!strcmp(qp_solver, "partial_condensing_hpmpc"))
    {
        plan->ocp_qp_solver_plan.qp_solver = PARTIAL_CONDENSING_HPMPC;
    }
#endif
#if defined(ACADOS_WITH_OSQP)
    else if (!strcmp(qp_solver, "partial_condensing_osqp"))
    {
        plan->ocp_qp_solver_plan.qp_solver = PARTIAL_CONDENSING_OSQP;
    }
#endif
    else
    {
        MEX_FIELD_VALUE_NOT_SUPPORTED_SUGGEST(fun_name, "qp_solver", qp_solver,
             "partial_condensing_hpipm, full_condensing_hpipm, full_condensing_qpoases, partial_condensing_osqp, partial_condensing_hpmpc");
    }


    // regularization
    char *regularize_method;
    if (mxGetField( matlab_opts, 0, "regularize_method" )!=NULL)
    {
        regularize_method = mxArrayToString( mxGetField( matlab_opts, 0, "regularize_method" ) );
    }
    else
    {
        MEX_MISSING_ARGUMENT(fun_name, "regularize_method");
    }

    if (!strcmp(regularize_method, "no_regularize"))
    {
        plan->regularization = NO_REGULARIZE;
    }
    else if (!strcmp(regularize_method, "mirror"))
    {
        plan->regularization = MIRROR;
    }
    else if (!strcmp(regularize_method, "project"))
    {
        plan->regularization = PROJECT;
    }
    else if (!strcmp(regularize_method, "project_reduc_hess"))
    {
        plan->regularization = PROJECT_REDUC_HESS;
    }
    else if (!strcmp(regularize_method, "convexify"))
    {
        plan->regularization = CONVEXIFY;
    }
    else
    {
        MEX_FIELD_VALUE_NOT_SUPPORTED_SUGGEST(fun_name, "regularize_method", regularize_method,
             "no_regularize, mirror, project, project_reduc_hess, convexify");
    }

    ocp_nlp_config *config = ocp_nlp_config_create(*plan);


    /* dims */
    int nx, nu;
    int nz = 0;
    int ny = 0;
    int ny_e = 0;
    int nbx;
    int nbx_0;
    int nbx_e;
    int nbu;
    int ng;
    int ng_e;
    int nh;
    int nh_e;
    int ns = 0;
    int ns_e = 0;
    int nsbu = 0;
    int nsbx = 0;
    int nsbx_e = 0;
    int nsg = 0;
    int nsg_e = 0;
    int nsh = 0;
    int nsh_e = 0;
    int nbxe_0;

    ocp_nlp_dims *dims = ocp_nlp_dims_create(config);

    // nx
    if (mxGetField( matlab_model, 0, "dim_nx" )!=NULL)
    {
        nx = mxGetScalar( mxGetField( prhs[0], 0, "dim_nx" ) );
    }
    else
    {
        MEX_MISSING_ARGUMENT(fun_name, "dim_nx");
    }

    i_ptr = (int *) malloc((N+1)*sizeof(int));
    for (int ii=0; ii<=N; ii++)
        i_ptr[ii] = nx;
    ocp_nlp_dims_set_opt_vars(config, dims, "nx", i_ptr);

    // nu
    if (mxGetField( matlab_model, 0, "dim_nu" )!=NULL)
    {
        nu = mxGetScalar( mxGetField( prhs[0], 0, "dim_nu" ) );
    }
    else
    {
        MEX_MISSING_ARGUMENT(fun_name, "dim_nu");
    }
    for (int ii=0; ii<N; ii++)
        i_ptr[ii] = nu;
    i_ptr[N] = 0;
    ocp_nlp_dims_set_opt_vars(config, dims, "nu", i_ptr);
    
    // nz
    if (mxGetField( prhs[0], 0, "dim_nz" )!=NULL)
    {
        nz = mxGetScalar( mxGetField( prhs[0], 0, "dim_nz" ) );
        for (int ii=0; ii<=N; ii++)
            i_ptr[ii] = nz;
        ocp_nlp_dims_set_opt_vars(config, dims, "nz", i_ptr);
    }
    free(i_ptr);

    // gnsf stuff
    if (!strcmp(sim_method, "irk_gnsf"))
    {
        int gnsf_nx1, gnsf_nz1, gnsf_nuhat, gnsf_ny, gnsf_nout;
        // nx1
        if (mxGetField( matlab_model, 0, "dim_gnsf_nx1" )!=NULL)
        {
            gnsf_nx1 = mxGetScalar( mxGetField( matlab_model, 0, "dim_gnsf_nx1" ) );
        }
        else
        {
            MEX_MISSING_ARGUMENT_MODULE(fun_name, "dim_gnsf_nx1", "irk_gnsf");
        }
        // nz1
        if (mxGetField( matlab_model, 0, "dim_gnsf_nz1" )!=NULL)
        {
            gnsf_nz1 = mxGetScalar( mxGetField( matlab_model, 0, "dim_gnsf_nz1" ) );
        }
        else
        {
            MEX_MISSING_ARGUMENT_MODULE(fun_name, "dim_gnsf_nz1", "irk_gnsf");
        }
        // nuhat
        if (mxGetField( matlab_model, 0, "dim_gnsf_nuhat" )!=NULL)
        {
            gnsf_nuhat = mxGetScalar( mxGetField( matlab_model, 0, "dim_gnsf_nuhat" ) );
        }
        else
        {
            MEX_MISSING_ARGUMENT_MODULE(fun_name, "dim_gnsf_nuhat", "irk_gnsf");
        }
        // ny
        if (mxGetField( matlab_model, 0, "dim_gnsf_ny" )!=NULL)
        {
            gnsf_ny = mxGetScalar( mxGetField( matlab_model, 0, "dim_gnsf_ny" ) );
        }
        else
        {
            MEX_MISSING_ARGUMENT_MODULE(fun_name, "dim_gnsf_ny", "irk_gnsf");
        }
        // nout
        if (mxGetField( matlab_model, 0, "dim_gnsf_nout" )!=NULL)
        {
            gnsf_nout = mxGetScalar( mxGetField( matlab_model, 0, "dim_gnsf_nout" ) );
        }
        else
        {
            MEX_MISSING_ARGUMENT_MODULE(fun_name, "dim_gnsf_nout", "irk_gnsf");
        }
        // assign
        for (int ii=0; ii<N; ii++)
        {
            ocp_nlp_dims_set_dynamics(config, dims, ii, "gnsf_nx1", &gnsf_nx1);
            ocp_nlp_dims_set_dynamics(config, dims, ii, "gnsf_nz1", &gnsf_nz1);
            ocp_nlp_dims_set_dynamics(config, dims, ii, "gnsf_nuhat", &gnsf_nuhat);
            ocp_nlp_dims_set_dynamics(config, dims, ii, "gnsf_ny", &gnsf_ny);
            ocp_nlp_dims_set_dynamics(config, dims, ii, "gnsf_nout", &gnsf_nout);
        }
    }


    // ny
    if (mxGetField( matlab_model, 0, "dim_ny" )!=NULL)
    {
        ny = mxGetScalar( mxGetField( matlab_model, 0, "dim_ny" ) );
        for (int ii=0; ii<N; ii++)
        {
            ocp_nlp_dims_set_cost(config, dims, ii, "ny", &ny);
        }
    }
    // ny_e
    if (mxGetField( matlab_model, 0, "dim_ny_e" )!=NULL)
    {
        ny_e = mxGetScalar( mxGetField( matlab_model, 0, "dim_ny_e" ) );
        ocp_nlp_dims_set_cost(config, dims, N, "ny", &ny_e);
    }

    // constraint dims
    // nbx_0
    if (mxGetField( matlab_model, 0, "dim_nbx_0" )!=NULL)
    {
        nbx_0 = mxGetScalar( mxGetField( matlab_model, 0, "dim_nbx_0" ) );
        ocp_nlp_dims_set_constraints(config, dims, 0, "nbx", &nbx_0);
    }
    // nbx
    if (mxGetField( matlab_model, 0, "dim_nbx" )!=NULL)
    {
        nbx = mxGetScalar( mxGetField( matlab_model, 0, "dim_nbx" ) );
        for (int ii=1; ii<N; ii++)
        {
            ocp_nlp_dims_set_constraints(config, dims, ii, "nbx", &nbx);
        }
    }
    // nbx_e
    if (mxGetField( matlab_model, 0, "dim_nbx_e" )!=NULL)
    {
        nbx_e = mxGetScalar( mxGetField( matlab_model, 0, "dim_nbx_e" ) );
        ocp_nlp_dims_set_constraints(config, dims, N, "nbx", &nbx_e);
    }
    // nbu
    if (mxGetField( matlab_model, 0, "dim_nbu" )!=NULL)
    {
        nbu = mxGetScalar( mxGetField( matlab_model, 0, "dim_nbu" ) );
        for (int ii=0; ii<N; ii++)
        {
            ocp_nlp_dims_set_constraints(config, dims, ii, "nbu", &nbu);
        }
    }
    // ng
    if (mxGetField( matlab_model, 0, "dim_ng" )!=NULL)
    {
        ng = mxGetScalar( mxGetField( matlab_model, 0, "dim_ng" ) );
        for (int ii=0; ii<N; ii++)
        {
            ocp_nlp_dims_set_constraints(config, dims, ii, "ng", &ng);
        }
    }
    // ng_e
    if (mxGetField( matlab_model, 0, "dim_ng_e" )!=NULL)
    {
        ng_e = mxGetScalar( mxGetField( matlab_model, 0, "dim_ng_e" ) );
        ocp_nlp_dims_set_constraints(config, dims, N, "ng", &ng_e);
    }

    // nh
    if (mxGetField( matlab_model, 0, "dim_nh" )!=NULL)
    {
        nh = mxGetScalar( mxGetField( matlab_model, 0, "dim_nh" ) );
        for (int ii=0; ii<N; ii++)
        {
            ocp_nlp_dims_set_constraints(config, dims, ii, "nh", &nh);
        }
    }
    // nh_e
    if (mxGetField( matlab_model, 0, "dim_nh_e" )!=NULL)
    {
        nh_e = mxGetScalar( mxGetField( matlab_model, 0, "dim_nh_e" ) );
        ocp_nlp_dims_set_constraints(config, dims, N, "nh", &nh_e);
    }

    // slack dims
    // ns
    if (mxGetField( matlab_model, 0, "dim_ns" )!=NULL)
    {
        ns = mxGetScalar( mxGetField( matlab_model, 0, "dim_ns" ) );
    }

    // ns_e
    if (mxGetField( matlab_model, 0, "dim_ns_e" )!=NULL)
    {
        ns_e = mxGetScalar( mxGetField( matlab_model, 0, "dim_ns_e" ) );
    }
    // nsbu
    if (mxGetField( matlab_model, 0, "dim_nsbu" )!=NULL)
    {
        nsbu = mxGetScalar( mxGetField( matlab_model, 0, "dim_nsbu" ) );
        for (int ii=0; ii<N; ii++)
        {
            ocp_nlp_dims_set_constraints(config, dims, ii, "nsbu", &nsbu);
        }
    }
    // nsbx
    if (mxGetField( matlab_model, 0, "dim_nsbx" )!=NULL)
    {
        nsbx = mxGetScalar( mxGetField( matlab_model, 0, "dim_nsbx" ) );
        for (int ii=1; ii<N; ii++) // TODO fix stage 0 !!!!!
        {
            ocp_nlp_dims_set_constraints(config, dims, ii, "nsbx", &nsbx);
        }
    }
    if (mxGetField( matlab_model, 0, "dim_nsbx_e" )!=NULL)
    {
        nsbx_e = mxGetScalar( mxGetField( matlab_model, 0, "dim_nsbx_e" ) );
        ocp_nlp_dims_set_constraints(config, dims, N, "nsbx", &nsbx);
    }
    // nsg
    if (mxGetField( matlab_model, 0, "dim_nsg" )!=NULL)
    {
        nsg = mxGetScalar( mxGetField( matlab_model, 0, "dim_nsg" ) );
        for (int ii=0; ii<N; ii++)
        {
            ocp_nlp_dims_set_constraints(config, dims, ii, "nsg", &nsg);
        }
    }
    // nsg_e
    if (mxGetField( matlab_model, 0, "dim_nsg_e" )!=NULL)
    {
        nsg_e = mxGetScalar( mxGetField( matlab_model, 0, "dim_nsg_e" ) );
        ocp_nlp_dims_set_constraints(config, dims, N, "nsg", &nsg_e);
    }
    // nsh
    if (mxGetField( matlab_model, 0, "dim_nsh" )!=NULL)
    {
        nsh = mxGetScalar( mxGetField( matlab_model, 0, "dim_nsh" ) );
        for (int ii=0; ii<N; ii++)
        {
            ocp_nlp_dims_set_constraints(config, dims, ii, "nsh", &nsh);
        }
    }
    // nsh_e
    if (mxGetField( matlab_model, 0, "dim_nsh_e" )!=NULL)
    {
        nsh_e = mxGetScalar( mxGetField( matlab_model, 0, "dim_nsh_e" ) );
        ocp_nlp_dims_set_constraints(config, dims, N, "nsh", &nsh_e);
    }
    // ns
    if (ns!=nsbu+nsbx+nsg+nsh)
    {
        sprintf(buffer,"ocp_create: ns!=nsbu+nsbx+nsg+nsh, got ns=%d, nsbu=%d, nsbx=%d nsg=%d nsh=%d\n",
            ns, nsbu, nsbx, nsg, nsh);
        mexErrMsgTxt(buffer);
    }
    if (ns_e!=nsbx_e+nsg_e+nsh_e)
    {
        sprintf(buffer,"ocp_create: ns_e!=nsbx_e+nsg_e+nsh_e, got ns_e=%d, nsbx_e=%d nsg_e=%d nsh_e=%d\n",
            ns_e, nsbx_e, nsg_e, nsh_e);
        mexErrMsgTxt(buffer);
    }

    i_ptr = (int *) malloc((N+1)*sizeof(int));
    // TODO fix stage 0 !!!!!!!!
    i_ptr[0] = nsbu+nsg+nsh; // XXX not nsbx !!!!!!!!!!
    for (int ii=1; ii<N; ii++)
        i_ptr[ii] = ns;
    i_ptr[N] = ns_e;
    ocp_nlp_dims_set_opt_vars(config, dims, "ns", i_ptr);
    free(i_ptr);

    // nbxe_0
    if (mxGetField( matlab_model, 0, "dim_nbxe_0" )!=NULL)
    {
        nbxe_0 = mxGetScalar( mxGetField( matlab_model, 0, "dim_nbxe_0" ) );
        ocp_nlp_dims_set_constraints(config, dims, 0, "nbxe", &nbxe_0);
    }


    /* opts */
    void *opts = ocp_nlp_solver_opts_create(config, dims);

    // nlp solver exact hessian
    bool nlp_solver_exact_hessian = false;
    c_ptr = mxArrayToString( mxGetField( matlab_opts, 0, "nlp_solver_exact_hessian" ) );
    MEX_STR_TO_BOOL(fun_name, c_ptr, &nlp_solver_exact_hessian, "nlp_solver_exact_hessian");
    // TODO: this if should not be needed! however, calling the setter with false leads to weird behavior. Investigate!
    if (nlp_solver_exact_hessian)
    {
        ocp_nlp_solver_opts_set(config, opts, "exact_hess", &nlp_solver_exact_hessian);
        if (sim_method_enum == GNSF)
            MEX_FIELD_NOT_SUPPORTED_GIVEN(fun_name, "nlp_solver_exact_hessian",
                 "true", "irk_gnsf", "false")

        // int exact_hess_dyn = {{ solver_options.exact_hess_dyn }};
        // ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "exact_hess_dyn", &exact_hess_dyn);
        // nlp solver max iter
        if (mxGetField( matlab_opts, 0, "exact_hess_dyn" )!=NULL)
        {
            int exact_hess_dyn = mxGetScalar( mxGetField( matlab_opts, 0, "exact_hess_dyn" ) );
            ocp_nlp_solver_opts_set(config, opts, "exact_hess_dyn", &exact_hess_dyn);
        }

        if (mxGetField( matlab_opts, 0, "exact_hess_cost" )!=NULL)
        {
            int exact_hess_cost = mxGetScalar( mxGetField( matlab_opts, 0, "exact_hess_cost" ) );
            ocp_nlp_solver_opts_set(config, opts, "exact_hess_cost", &exact_hess_cost);
        }

        if (mxGetField( matlab_opts, 0, "exact_hess_constr" )!=NULL)
        {
            int exact_hess_constr = mxGetScalar( mxGetField( matlab_opts, 0, "exact_hess_constr" ) );
            ocp_nlp_solver_opts_set(config, opts, "exact_hess_constr", &exact_hess_constr);
        }
    }

    if ( plan->nlp_solver != SQP_RTI )
    {
        // nlp solver max iter
        if (mxGetField( matlab_opts, 0, "nlp_solver_max_iter" )!=NULL)
        {
            int nlp_solver_max_iter = mxGetScalar( mxGetField( matlab_opts, 0, "nlp_solver_max_iter" ) );
            ocp_nlp_solver_opts_set(config, opts, "max_iter", &nlp_solver_max_iter);
        }

        // nlp solver exit tolerances
        if (mxGetField( matlab_opts, 0, "nlp_solver_tol_stat" )!=NULL)
        {
            double nlp_solver_tol_stat = mxGetScalar( mxGetField( matlab_opts, 0, "nlp_solver_tol_stat" ) );
            ocp_nlp_solver_opts_set(config, opts, "tol_stat", &nlp_solver_tol_stat);
        }
        if (mxGetField( matlab_opts, 0, "nlp_solver_tol_eq" )!=NULL)
        {
            double nlp_solver_tol_eq = mxGetScalar( mxGetField( matlab_opts, 0, "nlp_solver_tol_eq" ) );
            ocp_nlp_solver_opts_set(config, opts, "tol_eq", &nlp_solver_tol_eq);
        }
        if (mxGetField( matlab_opts, 0, "nlp_solver_tol_ineq" )!=NULL)
        {
            double nlp_solver_tol_ineq = mxGetScalar( mxGetField( matlab_opts, 0, "nlp_solver_tol_ineq" ) );
            ocp_nlp_solver_opts_set(config, opts, "tol_ineq", &nlp_solver_tol_ineq);
        }
        if (mxGetField( matlab_opts, 0, "nlp_solver_tol_comp" )!=NULL)
        {
            double nlp_solver_tol_comp = mxGetScalar( mxGetField( matlab_opts, 0, "nlp_solver_tol_comp" ) );
            ocp_nlp_solver_opts_set(config, opts, "tol_comp", &nlp_solver_tol_comp);
        }
    }

    // nlp solver ext qp res
    if (mxGetField( matlab_opts, 0, "nlp_solver_ext_qp_res" )!=NULL)
    {
        int nlp_solver_ext_qp_res = mxGetScalar( mxGetField( matlab_opts, 0, "nlp_solver_ext_qp_res" ) );
        ocp_nlp_solver_opts_set(config, opts, "ext_qp_res", &nlp_solver_ext_qp_res);
    }
    // nlp solver step length
    if (mxGetField( matlab_opts, 0, "nlp_solver_step_length" )!=NULL)
    {
        double nlp_solver_step_length = mxGetScalar( mxGetField( matlab_opts, 0, "nlp_solver_step_length" ) );
        ocp_nlp_solver_opts_set(config, opts, "step_length", &nlp_solver_step_length);
    }
    // RTI phase
    if (mxGetField( matlab_opts, 0, "rti_phase" )!=NULL)
    {
        int rti_phase = mxGetScalar( mxGetField( matlab_opts, 0, "rti_phase" ) );
        ocp_nlp_solver_opts_set(config, opts, "rti_phase", &rti_phase);
    }
    // nlp solver: warm start first
    if (mxGetField( matlab_opts, 0, "nlp_solver_warm_start_first_qp" )!=NULL)
    {
        int nlp_solver_warm_start_first_qp = mxGetScalar( mxGetField( matlab_opts, 0, "nlp_solver_warm_start_first_qp" ) );
        ocp_nlp_solver_opts_set(config, opts, "warm_start_first_qp", &nlp_solver_warm_start_first_qp);
    }
    // iter_max
    if (mxGetField( matlab_opts, 0, "qp_solver_iter_max" )!=NULL)
    {
        int qp_solver_iter_max = mxGetScalar( mxGetField( matlab_opts, 0, "qp_solver_iter_max" ) );
        ocp_nlp_solver_opts_set(config, opts, "qp_iter_max", &qp_solver_iter_max);
    }

	// qp solver exit tolerances
	if (mxGetField( matlab_opts, 0, "qp_solver_tol_stat" )!=NULL)
	{
		double qp_solver_tol_stat = mxGetScalar( mxGetField( matlab_opts, 0, "qp_solver_tol_stat" ) );
		ocp_nlp_solver_opts_set(config, opts, "qp_tol_stat", &qp_solver_tol_stat);
	}
	if (mxGetField( matlab_opts, 0, "qp_solver_tol_eq" )!=NULL)
	{
		double qp_solver_tol_eq = mxGetScalar( mxGetField( matlab_opts, 0, "qp_solver_tol_eq" ) );
		ocp_nlp_solver_opts_set(config, opts, "qp_tol_eq", &qp_solver_tol_eq);
	}
	if (mxGetField( matlab_opts, 0, "qp_solver_tol_ineq" )!=NULL)
	{
		double qp_solver_tol_ineq = mxGetScalar( mxGetField( matlab_opts, 0, "qp_solver_tol_ineq" ) );
		ocp_nlp_solver_opts_set(config, opts, "qp_tol_ineq", &qp_solver_tol_ineq);
	}
	if (mxGetField( matlab_opts, 0, "qp_solver_tol_comp" )!=NULL)
	{
		double qp_solver_tol_comp = mxGetScalar( mxGetField( matlab_opts, 0, "qp_solver_tol_comp" ) );
		ocp_nlp_solver_opts_set(config, opts, "qp_tol_comp", &qp_solver_tol_comp);
	}

    // N_part_cond
    if (mxGetField( matlab_opts, 0, "qp_solver_cond_N" )!=NULL)
    {
        int qp_solver_cond_N = mxGetScalar( mxGetField( matlab_opts, 0, "qp_solver_cond_N" ) );
        if (qp_solver_cond_N > N)
        {
            sprintf(buffer, "%s setting qp_solver_cond_N = %d not possible, since N = %d. Require qp_solver_cond_N <= N",
                    fun_name, qp_solver_cond_N, N);
            mexErrMsgTxt(buffer);
        }
        else if ( plan->ocp_qp_solver_plan.qp_solver == PARTIAL_CONDENSING_HPIPM )
        {
            ocp_nlp_solver_opts_set(config, opts, "qp_cond_N", &qp_solver_cond_N);
        }
        #if defined( ACADOS_WITH_HPMPC )
        else if ( plan->ocp_qp_solver_plan.qp_solver == PARTIAL_CONDENSING_HPMPC )
        {
            ocp_nlp_solver_opts_set(config, opts, "qp_cond_N", &qp_solver_cond_N);
        }
        #endif
        #if defined( ACADOS_WITH_OOQP )
        else if ( plan->ocp_qp_solver_plan.qp_solver == PARTIAL_CONDENSING_OOQP )
        {
            ocp_nlp_solver_opts_set(config, opts, "qp_cond_N", &qp_solver_cond_N);
        }
        #endif
        #if defined( ACADOS_WITH_OSQP )
        else if ( plan->ocp_qp_solver_plan.qp_solver == PARTIAL_CONDENSING_OSQP )
        {
            ocp_nlp_solver_opts_set(config, opts, "qp_cond_N", &qp_solver_cond_N);
        }
        #endif
        #if defined( ACADOS_WITH_QPDUNES )
        else if ( plan->ocp_qp_solver_plan.qp_solver == PARTIAL_CONDENSING_QPDUNES )
        {
            ocp_nlp_solver_opts_set(config, opts, "qp_cond_N", &qp_solver_cond_N);
        }
        #endif
    }
    else if (plan->ocp_qp_solver_plan.qp_solver == PARTIAL_CONDENSING_HPIPM)
    {
        MEX_MISSING_ARGUMENT_MODULE(fun_name, "qp_solver_cond_N", "partial_condensing_hpipm");
    }
    // cond riccati-like algorithm
    if (mxGetField( matlab_opts, 0, "qp_solver_cond_ric_alg" )!=NULL)
    {
        int qp_solver_cond_ric_alg = mxGetScalar( mxGetField( matlab_opts, 0, "qp_solver_cond_ric_alg" ) );
        ocp_nlp_solver_opts_set(config, opts, "qp_cond_ric_alg", &qp_solver_cond_ric_alg);
    }
    // hpipm: riccati algorithm
    if (mxGetField( matlab_opts, 0, "qp_solver_ric_alg" )!=NULL
        && plan->ocp_qp_solver_plan.qp_solver == PARTIAL_CONDENSING_HPIPM)
    {
        int qp_solver_ric_alg = mxGetScalar( mxGetField( matlab_opts, 0, "qp_solver_ric_alg" ) );
        ocp_nlp_solver_opts_set(config, opts, "qp_ric_alg", &qp_solver_ric_alg);
    }
    // qp solver: warm start
    if (mxGetField( matlab_opts, 0, "qp_solver_warm_start" )!=NULL)
    {
        int qp_solver_warm_start = mxGetScalar( mxGetField( matlab_opts, 0, "qp_solver_warm_start" ) );
        ocp_nlp_solver_opts_set(config, opts, "qp_warm_start", &qp_solver_warm_start);
    }
    // print_level
    if (mxGetField( matlab_opts, 0, "print_level" )!=NULL)
    {
        int print_level = mxGetScalar( mxGetField( matlab_opts, 0, "print_level" ) );
        ocp_nlp_solver_opts_set(config, opts, "print_level", &print_level);
    }

    if (strcmp(dyn_type, "discrete"))
    {
        // sim_method_num_stages
        sprintf(matlab_field_name, "sim_method_num_stages");
        const mxArray *matlab_array;
        matlab_array = mxGetField( matlab_opts, 0, matlab_field_name );
        if (matlab_array!=NULL)
        {
            int matlab_size = (int) mxGetNumberOfElements( matlab_array );
            MEX_DIM_CHECK_VEC_TWO(fun_name, matlab_field_name, matlab_size, 1, N);

            if (matlab_size == 1)
            {
                int sim_method_num_stages = mxGetScalar( matlab_array );
                for (int ii=0; ii<N; ii++)
                {
                    ocp_nlp_solver_opts_set_at_stage(config, opts, ii, "dynamics_num_stages", &sim_method_num_stages);
                }
            }
            else
            {
                double *values = mxGetPr(matlab_array);
                for (int ii=0; ii<N; ii++)
                {
                    int sim_method_num_stages = (int) values[ii];
                    // mexPrintf("\nsim_method_num_stages[%d] = %d", ii, sim_method_num_stages);
                    ocp_nlp_solver_opts_set_at_stage(config, opts, ii, "dynamics_num_stages", &sim_method_num_stages);
                }
            }
        }
        // sim_method_num_steps
        sprintf(matlab_field_name, "sim_method_num_steps");
        matlab_array = mxGetField( matlab_opts, 0, matlab_field_name );
        if (matlab_array!=NULL)
        {
            int matlab_size = (int) mxGetNumberOfElements( matlab_array );
            MEX_DIM_CHECK_VEC_TWO(fun_name, matlab_field_name, matlab_size, 1, N);

            if (matlab_size == 1)
            {
                int sim_method_num_steps = mxGetScalar( matlab_array );
                for (int ii=0; ii<N; ii++)
                {
                    ocp_nlp_solver_opts_set_at_stage(config, opts, ii, "dynamics_num_steps", &sim_method_num_steps);
                }
            }
            else
            {
                double *values = mxGetPr(matlab_array);
                for (int ii=0; ii<N; ii++)
                {
                    int sim_method_num_steps = (int) values[ii];
                    // mexPrintf("\nsim_method_num_steps[%d] = %d", ii, sim_method_num_steps);
                    ocp_nlp_solver_opts_set_at_stage(config, opts, ii, "dynamics_num_steps", &sim_method_num_steps);
                }
            }
        }
        // sim_method_newton_iter
        sprintf(matlab_field_name, "sim_method_newton_iter");
        matlab_array = mxGetField( matlab_opts, 0, matlab_field_name );
        if (matlab_array!=NULL)
        {
            int matlab_size = (int) mxGetNumberOfElements( matlab_array );
            MEX_DIM_CHECK_VEC_TWO(fun_name, matlab_field_name, matlab_size, 1, N);

            if (matlab_size == 1)
            {
                int sim_method_newton_iter = mxGetScalar( matlab_array );
                for (int ii=0; ii<N; ii++)
                {
                    ocp_nlp_solver_opts_set_at_stage(config, opts, ii, "dynamics_newton_iter", &sim_method_newton_iter);
                }
            }
            else
            {
                double *values = mxGetPr(matlab_array);
                for (int ii=0; ii<N; ii++)
                {
                    int sim_method_newton_iter = (int) values[ii];
                    // mexPrintf("\nsim_method_newton_iter[%d] = %d", ii, sim_method_newton_iter);
                    ocp_nlp_solver_opts_set_at_stage(config, opts, ii, "dynamics_newton_iter", &sim_method_newton_iter);
                }
            }
        }
        // sim_method_jac_reuse
        sprintf(matlab_field_name, "sim_method_jac_reuse");
        matlab_array = mxGetField( matlab_opts, 0, matlab_field_name );
        if (matlab_array!=NULL)
        {
            int matlab_size = (int) mxGetNumberOfElements( matlab_array );
            MEX_DIM_CHECK_VEC_TWO(fun_name, matlab_field_name, matlab_size, 1, N);

            if (matlab_size == 1)
            {
                bool sim_method_jac_reuse = mxGetScalar( matlab_array );
                for (int ii=0; ii<N; ii++)
                {
                    ocp_nlp_solver_opts_set_at_stage(config, opts, ii, "dynamics_jac_reuse", &sim_method_jac_reuse);
                }
            }
            else
            {
                double *values = mxGetPr(matlab_array);
                for (int ii=0; ii<N; ii++)
                {
                    bool sim_method_jac_reuse = (int) values[ii];
                    // mexPrintf("\nsim_method_jac_reuse[%d] = %d", ii, sim_method_jac_reuse);
                    ocp_nlp_solver_opts_set_at_stage(config, opts, ii, "dynamics_jac_reuse", &sim_method_jac_reuse);
                }
            }
        }
        // sim_method_exact_z_output
        sprintf(matlab_field_name, "sim_method_exact_z_output");
        matlab_array = mxGetField( matlab_opts, 0, matlab_field_name );
        if (matlab_array!=NULL)
        {
            int matlab_size = (int) mxGetNumberOfElements( matlab_array );
            MEX_DIM_CHECK_VEC_TWO(fun_name, matlab_field_name, matlab_size, 1, N);

            if (matlab_size == 1)
            {
                bool sim_method_exact_z_output = mxGetScalar( matlab_array );
                for (int ii=0; ii<N; ii++)
                {
                    ocp_nlp_solver_opts_set_at_stage(config, opts, ii, "dynamics_exact_z_output", &sim_method_exact_z_output);
                }
            }
            else
            {
                double *values = mxGetPr(matlab_array);
                for (int ii=0; ii<N; ii++)
                {
                    bool sim_method_exact_z_output = (int) values[ii];
                    // mexPrintf("\nsim_method_exact_z_output[%d] = %d", ii, sim_method_exact_z_output);
                    ocp_nlp_solver_opts_set_at_stage(config, opts, ii, "dynamics_exact_z_output", &sim_method_exact_z_output);
                }
            }
        }
    }

    // levenberg_marquardt regularization
    if (mxGetField( matlab_opts, 0, "levenberg_marquardt" )!=NULL)
    {
        double levenberg_marquardt = mxGetScalar( mxGetField( matlab_opts, 0, "levenberg_marquardt" ) );
        ocp_nlp_solver_opts_set(config, opts, "levenberg_marquardt", &levenberg_marquardt);
    }


    /* in */
    ocp_nlp_in *in = ocp_nlp_in_create(config, dims);

    // time T
    double T;
    if (mxGetField( matlab_model, 0, "T" )!=NULL)
    {
        T = mxGetScalar( mxGetField( matlab_model, 0, "T" ) );
    }
    else
    {
        MEX_MISSING_ARGUMENT(fun_name, "T");
    }

    /* discretization grid */
    // time_steps
    double *time_steps;

    if (mxGetField( matlab_opts, 0, "time_steps" )!=NULL)
    {
        int matlab_size = (int) mxGetNumberOfElements( mxGetField( matlab_opts, 0, "time_steps" ));
        int acados_size = N;
        MEX_DIM_CHECK_VEC(fun_name, "time_steps", matlab_size, acados_size);

        time_steps = mxGetPr( mxGetField( matlab_opts, 0, "time_steps" ) );
        for (int ii=0; ii<N; ii++)
        {
            ocp_nlp_in_set(config, dims, in, ii, "Ts", &time_steps[ii]);
            ocp_nlp_cost_model_set(config, dims, in, ii, "scaling", &time_steps[ii]);
        }
    }
    else
    {
        MEX_MISSING_ARGUMENT(fun_name, "time_steps")
    }


    /* out */
    ocp_nlp_out *out = ocp_nlp_out_create(config, dims);

    /* solver */
    ocp_nlp_solver *solver = ocp_nlp_solver_create(config, dims, opts);

    /* sens_out */
    ocp_nlp_out *sens_out = ocp_nlp_out_create(config, dims);

    /* populate output struct */
    // plan
    mxArray *plan_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(plan_mat);
    l_ptr[0] = (long long) plan;
    mxSetField(plhs[0], 0, "plan", plan_mat);

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

    // sens_out
    mxArray *sens_out_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(sens_out_mat);
    l_ptr[0] = (long long) sens_out;
    mxSetField(plhs[0], 0, "sens_out", sens_out_mat);

    /* NUMERICAL DATA */
    /* cost */
    // lagrange term
    if ((cost_type_enum == LINEAR_LS) || (cost_type_enum == NONLINEAR_LS))
    {

        const mxArray *W_matlab = mxGetField( matlab_model, 0, "cost_W" );
        if (W_matlab!=NULL)
        {
            int nrow = (int) mxGetM( W_matlab );
            int ncol = (int) mxGetN( W_matlab );
            MEX_DIM_CHECK_MAT(fun_name, "cost_W", nrow, ncol, ny, ny);
            double *W = mxGetPr( W_matlab );
            for (int ii=0; ii<N; ii++)
            {
                ocp_nlp_cost_model_set(config, dims, in, ii, "W", W);
            }
        }
        else
        {
            MEX_MISSING_ARGUMENT_MODULE(fun_name, "cost_W", "linear_ls or nonlinear_ls");
        }

        if (mxGetField( matlab_model, 0, "cost_y_ref" )!=NULL)
        {
            int matlab_size = (int) mxGetNumberOfElements( mxGetField( matlab_model, 0, "cost_y_ref" ) );
            int acados_size = ny;
            double *yr = mxGetPr( mxGetField( matlab_model, 0, "cost_y_ref" ) );
            for (int ii=0; ii<N; ii++)
            {
                ocp_nlp_cost_model_set(config, dims, in, ii, "y_ref", yr);
            }
        }
        // else: is set to zero by default in C
    }
    if (cost_type_enum == LINEAR_LS)
    {
        const mxArray *Vu_matlab = mxGetField( matlab_model, 0, "cost_Vu" );
        if (Vu_matlab!=NULL)
        {
            int nrow = (int) mxGetM( Vu_matlab );
            int ncol = (int) mxGetN( Vu_matlab );
            MEX_DIM_CHECK_MAT(fun_name, "cost_Vu", nrow, ncol, ny, nu);
            double *Vu = mxGetPr( Vu_matlab );
            if (nrow > 0 && ncol > 0)
            {
                for (int ii=0; ii<N; ii++)
                {
                    ocp_nlp_cost_model_set(config, dims, in, ii, "Vu", Vu);
                }
            }
        }
        // else: set to zero by default in C
        const mxArray *Vx_matlab = mxGetField( matlab_model, 0, "cost_Vx" );
        if (Vx_matlab!=NULL)
        {
            int nrow = (int) mxGetM( Vx_matlab );
            int ncol = (int) mxGetN( Vx_matlab );
            MEX_DIM_CHECK_MAT(fun_name, "cost_Vx", nrow, ncol, ny, nx);
            double *Vx = mxGetPr( Vx_matlab );
            for (int ii=0; ii<N; ii++)
            {
                ocp_nlp_cost_model_set(config, dims, in, ii, "Vx", Vx);
            }
        }
        // else: set to zero by default in C
        const mxArray *Vz_matlab = mxGetField( matlab_model, 0, "cost_Vz" );
        if (Vz_matlab!=NULL)
        {
            int nrow = (int) mxGetM( Vz_matlab );
            int ncol = (int) mxGetN( Vz_matlab );
            MEX_DIM_CHECK_MAT(fun_name, "cost_Vz", nrow, ncol, ny, nz);
            double *Vz = mxGetPr( Vz_matlab );
            if (nrow > 0 && ncol > 0)
            {
                for (int ii=0; ii<N; ii++)
                {
                    ocp_nlp_cost_model_set(config, dims, in, ii, "Vz", Vz);
                }
            }
        }
        // else: set to zero by default in C
        if ((Vz_matlab==NULL) && (Vx_matlab==NULL) && (Vu_matlab==NULL))
        {
            sprintf(buffer, "%s: setting linear_ls cost: one of the following must be provided:\
                             cost_Vx, cost_Vu, cost_Vz", fun_name);
            mexErrMsgTxt(buffer);
        }
    }

    // mayer term
    if ((cost_type_e_enum == LINEAR_LS) || (cost_type_e_enum == NONLINEAR_LS))
    {
        const mxArray *W_e_matlab = mxGetField( matlab_model, 0, "cost_W_e" );
        if (W_e_matlab!=NULL)
        {
            int nrow = (int) mxGetM( W_e_matlab );
            int ncol = (int) mxGetN( W_e_matlab );
            MEX_DIM_CHECK_MAT(fun_name, "cost_W_e", nrow, ncol, ny_e, ny_e);
            double *W_e = mxGetPr( W_e_matlab );
            ocp_nlp_cost_model_set(config, dims, in, N, "W", W_e);
        }
        else
        {
            MEX_MISSING_ARGUMENT_MODULE(fun_name, "cost_W_e", "linear_ls or nonlinear_ls");
        }
        if (mxGetField( matlab_model, 0, "cost_y_ref_e" )!=NULL)
        {
            int matlab_size = (int) mxGetNumberOfElements( mxGetField( matlab_model, 0, "cost_y_ref_e" ) );
            int acados_size = ny_e;
            MEX_DIM_CHECK_VEC(fun_name, "cost_y_ref_e", matlab_size, acados_size);
            double *yr_e = mxGetPr( mxGetField( matlab_model, 0, "cost_y_ref_e" ) );
            ocp_nlp_cost_model_set(config, dims, in, N, "y_ref", yr_e);
        }
        // else set to zero by default in C
    }
    if (cost_type_e_enum == LINEAR_LS)
    {
        const mxArray *Vx_e_matlab = mxGetField( matlab_model, 0, "cost_Vx_e" );
        if (Vx_e_matlab!=NULL)
        {
            int nrow = (int) mxGetM( Vx_e_matlab );
            int ncol = (int) mxGetN( Vx_e_matlab );
            MEX_DIM_CHECK_MAT(fun_name, "cost_Vx_e", nrow, ncol, ny_e, nx);
            double *Vx_e = mxGetPr( Vx_e_matlab );
            ocp_nlp_cost_model_set(config, dims, in, N, "Vx", Vx_e);
        }
        else
        {
            MEX_MISSING_ARGUMENT_MODULE(fun_name, "cost_Vx_e", "cost_type_e: linear_ls");
        }
        
    }

    // slacks
    if (mxGetField( matlab_model, 0, "cost_Z" )!=NULL)
    {
        int matlab_size = (int) mxGetNumberOfElements( mxGetField( matlab_model, 0, "cost_Z" ) );
        int acados_size = ns*ns;
        MEX_DIM_CHECK_VEC(fun_name, "cost_Z", matlab_size, acados_size);
        double *Z = mxGetPr( mxGetField( matlab_model, 0, "cost_Z" ) );
        MEX_CHECK_DIAGONALITY(fun_name, ns, Z, "cost_Z");

        d_ptr = malloc(ns*sizeof(double));
        for (int ii=0; ii<ns; ii++)
        {
            d_ptr[ii] = Z[ii+ns*ii];
        }
        for (int ii=0; ii<N; ii++)
        {
            ocp_nlp_cost_model_set(config, dims, in, ii, "Z", d_ptr);
        }
        free(d_ptr);
    }

    if (mxGetField( matlab_model, 0, "cost_Z_e" )!=NULL)
    {
        int matlab_size = (int) mxGetNumberOfElements( mxGetField( matlab_model, 0, "cost_Z_e" ) );
        int acados_size = ns_e*ns_e;
        MEX_DIM_CHECK_VEC(fun_name, "cost_Z_e", matlab_size, acados_size);
        double *Z_e = mxGetPr( mxGetField( matlab_model, 0, "cost_Z_e" ) );
        MEX_CHECK_DIAGONALITY(fun_name, ns_e, Z_e, "cost_Z_e");
        d_ptr = malloc(ns_e*sizeof(double));
        for (int ii=0; ii<ns_e; ii++)
        {
            d_ptr[ii] = Z_e[ii+ns_e*ii];
        }
        ocp_nlp_cost_model_set(config, dims, in, N, "Z", d_ptr);
        free(d_ptr);
    }

    if (mxGetField( matlab_model, 0, "cost_Zl" )!=NULL)
    {
        int matlab_size = (int) mxGetNumberOfElements( mxGetField( matlab_model, 0, "cost_Zl" ) );
        int acados_size = ns*ns;
        MEX_DIM_CHECK_VEC(fun_name, "cost_Zl", matlab_size, acados_size);
        double *Zl = mxGetPr( mxGetField( matlab_model, 0, "cost_Zl" ) );
        d_ptr = malloc(ns*sizeof(double));
        for (int ii=0; ii<ns; ii++)
        {
            d_ptr[ii] = Zl[ii+ns*ii];
        }
        for (int ii=0; ii<N; ii++)
        {
            ocp_nlp_cost_model_set(config, dims, in, ii, "Zl", d_ptr);
        }
        free(d_ptr);
    }
    if (mxGetField( matlab_model, 0, "cost_Zl_e" )!=NULL)
    {
        int matlab_size = (int) mxGetNumberOfElements( mxGetField( matlab_model, 0, "cost_Zl_e" ) );
        int acados_size = ns_e*ns_e;
        MEX_DIM_CHECK_VEC(fun_name, "cost_Zl_e", matlab_size, acados_size);
        double *Zl_e = mxGetPr( mxGetField( matlab_model, 0, "cost_Zl_e" ) );
        d_ptr = malloc(ns_e*sizeof(double));
        for (int ii=0; ii<ns_e; ii++)
        {
            d_ptr[ii] = Zl_e[ii+ns_e*ii];
        }
        ocp_nlp_cost_model_set(config, dims, in, N, "Zl", d_ptr);
        free(d_ptr);
    }

    if (mxGetField( matlab_model, 0, "cost_Zu" )!=NULL)
    {
        int matlab_size = (int) mxGetNumberOfElements( mxGetField( matlab_model, 0, "cost_Zu" ) );
        int acados_size = ns*ns;
        MEX_DIM_CHECK_VEC(fun_name, "cost_Zu", matlab_size, acados_size);
        double *Zu = mxGetPr( mxGetField( matlab_model, 0, "cost_Zu" ) );
        d_ptr = malloc(ns*sizeof(double));
        for (int ii=0; ii<ns; ii++)
        {
            d_ptr[ii] = Zu[ii+ns*ii];
        }
        for (int ii=0; ii<N; ii++)
        {
            ocp_nlp_cost_model_set(config, dims, in, ii, "Zu", d_ptr);
        }
        free(d_ptr);
    }
    if (mxGetField( matlab_model, 0, "cost_Zu_e" )!=NULL)
    {
        int matlab_size = (int) mxGetNumberOfElements( mxGetField( matlab_model, 0, "cost_Zu_e" ) );
        int acados_size = ns_e*ns_e;
        MEX_DIM_CHECK_VEC(fun_name, "cost_Zu_e", matlab_size, acados_size);
        double *Zu_e = mxGetPr( mxGetField( matlab_model, 0, "cost_Zu_e" ) );
        d_ptr = malloc(ns_e*sizeof(double));
        for (int ii=0; ii<ns_e; ii++)
        {
            d_ptr[ii] = Zu_e[ii+ns_e*ii];
        }
        ocp_nlp_cost_model_set(config, dims, in, N, "Zu", d_ptr);
        free(d_ptr);
    }
    // z
    if (mxGetField( matlab_model, 0, "cost_z" )!=NULL)
    {
        int matlab_size = (int) mxGetNumberOfElements( mxGetField( matlab_model, 0, "cost_z" ) );
        int acados_size = ns;
        MEX_DIM_CHECK_VEC(fun_name, "cost_z", matlab_size, acados_size);
        double *z = mxGetPr( mxGetField( matlab_model, 0, "cost_z" ) );
        for (int ii=0; ii<N; ii++)
        {
            ocp_nlp_cost_model_set(config, dims, in, ii, "z", z);
        }
    }
    // z_e
    if (mxGetField( matlab_model, 0, "cost_z_e" )!=NULL)
    {
        int matlab_size = (int) mxGetNumberOfElements( mxGetField( matlab_model, 0, "cost_z_e" ) );
        int acados_size = ns_e;
        MEX_DIM_CHECK_VEC(fun_name, "cost_z_e", matlab_size, acados_size);
        double *z_e = mxGetPr( mxGetField( matlab_model, 0, "cost_z_e" ) );
        ocp_nlp_cost_model_set(config, dims, in, N, "z", z_e);
    }
    // zl
    if (mxGetField( matlab_model, 0, "cost_zl" )!=NULL)
    {
        int matlab_size = (int) mxGetNumberOfElements( mxGetField( matlab_model, 0, "cost_zl" ) );
        int acados_size = ns;
        MEX_DIM_CHECK_VEC(fun_name, "cost_zl", matlab_size, acados_size);
        double *zl = mxGetPr( mxGetField( matlab_model, 0, "cost_zl" ) );
        for (int ii=0; ii<N; ii++)
        {
            ocp_nlp_cost_model_set(config, dims, in, ii, "zl", zl);
        }
    }
    // zl_e
    if (mxGetField( matlab_model, 0, "cost_zl_e" )!=NULL)
    {
        int matlab_size = (int) mxGetNumberOfElements( mxGetField( matlab_model, 0, "cost_zl_e" ) );
        int acados_size = ns_e;
        MEX_DIM_CHECK_VEC(fun_name, "cost_zl_e", matlab_size, acados_size);
        double *zl_e = mxGetPr( mxGetField( matlab_model, 0, "cost_zl_e" ) );
        ocp_nlp_cost_model_set(config, dims, in, N, "zl", zl_e);
    }
    // zu
    if (mxGetField( matlab_model, 0, "cost_zu" )!=NULL)
    {
        int matlab_size = (int) mxGetNumberOfElements( mxGetField( matlab_model, 0, "cost_zu" ) );
        int acados_size = ns;
        MEX_DIM_CHECK_VEC(fun_name, "cost_zu", matlab_size, acados_size);
        double *zu = mxGetPr( mxGetField( matlab_model, 0, "cost_zu" ) );
        for (int ii=0; ii<N; ii++)
        {
            ocp_nlp_cost_model_set(config, dims, in, ii, "zu", zu);
        }
    }
    // zu_e
    if (mxGetField( matlab_model, 0, "cost_zu_e" )!=NULL)
    {
        int matlab_size = (int) mxGetNumberOfElements( mxGetField( matlab_model, 0, "cost_zu_e" ) );
        int acados_size = ns_e;
        MEX_DIM_CHECK_VEC(fun_name, "cost_zu_e", matlab_size, acados_size);
        double *zu_e = mxGetPr( mxGetField( matlab_model, 0, "cost_zu_e" ) );
        ocp_nlp_cost_model_set(config, dims, in, N, "zu", zu_e);
    }


    /* dynamics */
    // trajectory initialization
    // u_init
    if (mxGetField( matlab_model, 0, "u_init" )!=NULL)
    {
        double *u_init = mxGetPr( mxGetField( matlab_model, 0, "u_init" ) );
        for (int ii=0; ii<N; ii++)
        {
            ocp_nlp_out_set(config, dims, out, ii, "u", u_init+ii*nu);
        }
    }
    else // initialize to zero
    {
        double *u_init = calloc(nu, sizeof(double));
        for (int ii=0; ii<N; ii++)
        {
            ocp_nlp_out_set(config, dims, out, ii, "u", u_init);
        }
        free(u_init);
    }
    // x_init
    if (mxGetField( matlab_model, 0, "x_init" )!=NULL)
    {
        double *x_init = mxGetPr( mxGetField( matlab_model, 0, "x_init" ) );
        for (int ii=0; ii<=N; ii++)
        {
            ocp_nlp_out_set(config, dims, out, ii, "x", x_init+ii*nx);
        }
    }
    else // initialize to zero
    {
        double *x_init = calloc(nx, sizeof(double));
        for (int ii=0; ii<=N; ii++)
        {
            ocp_nlp_out_set(config, dims, out, ii, "x", x_init);
        }
        free(x_init);
    }


    /* constraints */
    // lbx
    double *lbx;    bool set_lbx = false;
    if (mxGetField( matlab_model, 0, "constr_lbx" )!=NULL)
    {
        int matlab_size = (int) mxGetNumberOfElements( mxGetField( matlab_model, 0, "constr_lbx" ) );
        int acados_size = nbx;
        MEX_DIM_CHECK_VEC(fun_name, "constr_lbx", matlab_size, acados_size);
        set_lbx = true;
        lbx = mxGetPr( mxGetField( matlab_model, 0, "constr_lbx" ) );
        for (int ii=1; ii<N; ii++)
        {
            ocp_nlp_constraints_model_set(config, dims, in, ii, "lbx", lbx);
        }
    }
    // ubx
    double *ubx;    bool set_ubx = false;
    if (mxGetField( matlab_model, 0, "constr_ubx" )!=NULL)
    {
        int matlab_size = (int) mxGetNumberOfElements( mxGetField( matlab_model, 0, "constr_ubx" ) );
        int acados_size = nbx;
        MEX_DIM_CHECK_VEC(fun_name, "constr_ubx", matlab_size, acados_size);
        set_ubx = true;
        ubx = mxGetPr( mxGetField( matlab_model, 0, "constr_ubx" ) );
        for (int ii=1; ii<N; ii++)
        {
            ocp_nlp_constraints_model_set(config, dims, in, ii, "ubx", ubx);
        }
    }

    // lbx_e
    double *lbx_e;
    if (mxGetField( matlab_model, 0, "constr_lbx_e" )!=NULL)
    {
        int matlab_size = (int) mxGetNumberOfElements( mxGetField( matlab_model, 0, "constr_lbx_e" ) );
        int acados_size = nbx_e;
        MEX_DIM_CHECK_VEC(fun_name, "constr_lbx_e", matlab_size, acados_size);
        lbx_e = mxGetPr( mxGetField( matlab_model, 0, "constr_lbx_e" ) );
        ocp_nlp_constraints_model_set(config, dims, in, N, "lbx", lbx_e);
    }
    // ubx_e
    double *ubx_e;
    if (mxGetField( matlab_model, 0, "constr_ubx_e" )!=NULL)
    {
        int matlab_size = (int) mxGetNumberOfElements( mxGetField( matlab_model, 0, "constr_ubx_e" ) );
        int acados_size = nbx_e;
        MEX_DIM_CHECK_VEC(fun_name, "constr_ubx_e", matlab_size, acados_size);
        ubx_e = mxGetPr( mxGetField( matlab_model, 0, "constr_ubx_e" ) );
        ocp_nlp_constraints_model_set(config, dims, in, N, "ubx", ubx_e);
    }

    // Jbx_0
    tmp_idx = malloc(nbx_0*sizeof(int));

    double *Jbx_0;
    const mxArray *Jbx_0_matlab = mxGetField( matlab_model, 0, "constr_Jbx_0" );
    if (Jbx_0_matlab!=NULL)
    {
        int nrow = (int) mxGetM( Jbx_0_matlab );
        int ncol = (int) mxGetN( Jbx_0_matlab );
        MEX_DIM_CHECK_MAT(fun_name, "constr_Jbx_0", nrow, ncol, nbx_0, nx);
        Jbx_0 = mxGetPr( Jbx_0_matlab );
        for (int ii=0; ii<nrow; ii++)
        {
            int nnz_row = 0;
            for (int jj=0; jj<ncol; jj++)
            {
                if (Jbx_0[ii+nrow*jj]==1.0)
                {
                    tmp_idx[ii] = jj;
                    nnz_row++;
                }
                else if (Jbx_0[ii+nrow*jj]!=0.0)
                {
                    MEX_NONBINARY_MAT(fun_name, "constr_Jbx_0");
                }
            }
            if (nnz_row > 1)
            {
                MEX_MULTIPLE_ONES_IN_ROW(fun_name, "constr_Jbx_0");
            }
        }
        ocp_nlp_constraints_model_set(config, dims, in, 0, "idxbx", tmp_idx);
    }
    else
    {
        MEX_MISSING_ARGUMENT(fun_name, "constr_Jbx_0");
    }
    free(tmp_idx);


    // Jbx
    tmp_idx = malloc(nbx*sizeof(int));
    double *Jbx;    bool set_Jbx = false;
    const mxArray *Jbx_matlab = mxGetField( matlab_model, 0, "constr_Jbx" );
    if (Jbx_matlab!=NULL)
    {
        int nrow = (int) mxGetM( Jbx_matlab );
        int ncol = (int) mxGetN( Jbx_matlab );
        MEX_DIM_CHECK_MAT(fun_name, "constr_Jbx", nrow, ncol, nbx, nx);
        set_Jbx = true;
        Jbx = mxGetPr( Jbx_matlab );
        for (int ii=0; ii<nrow; ii++)
        {
            int nnz_row = 0;
            for (int jj=0; jj<ncol; jj++)
            {
                if (Jbx[ii+nrow*jj]==1.0)
                {
                    tmp_idx[ii] = jj;
                    nnz_row++;
                }
                else if (Jbx[ii+nrow*jj]!=0.0)
                {
                    MEX_NONBINARY_MAT(fun_name, "constr_Jbx");
                }
            }
            if (nnz_row > 1)
            {
                MEX_MULTIPLE_ONES_IN_ROW(fun_name, "constr_Jbx");
            }
        }
        for (int ii=1; ii<N; ii++)
        {
            ocp_nlp_constraints_model_set(config, dims, in, ii, "idxbx", tmp_idx);
        }
    }
    free(tmp_idx);

    // lbx_0
    if (mxGetField( matlab_model, 0, "constr_lbx_0" )!=NULL)
    {
        int matlab_size = (int) mxGetNumberOfElements( mxGetField( matlab_model, 0, "constr_lbx_0" ) );
        int acados_size = nbx_0;
        MEX_DIM_CHECK_VEC(fun_name, "constr_lbx_0", matlab_size, acados_size);
        double *lbx_0 = mxGetPr( mxGetField( matlab_model, 0, "constr_lbx_0" ) );

        ocp_nlp_constraints_model_set(config, dims, in, 0, "lbx", lbx_0);
    }
    else
    {
        MEX_MISSING_ARGUMENT_NOTE(fun_name, "constr_lbx_0", "can be updated after creation");
    }
    // ubx_0
    if (mxGetField( matlab_model, 0, "constr_ubx_0" )!=NULL)
    {
        int matlab_size = (int) mxGetNumberOfElements( mxGetField( matlab_model, 0, "constr_ubx_0" ) );
        int acados_size = nbx_0;
        MEX_DIM_CHECK_VEC(fun_name, "constr_ubx_0", matlab_size, acados_size);
        double *ubx_0 = mxGetPr( mxGetField( matlab_model, 0, "constr_ubx_0" ) );

        ocp_nlp_constraints_model_set(config, dims, in, 0, "ubx", ubx_0);
    }
    else
    {
        MEX_MISSING_ARGUMENT_NOTE(fun_name, "constr_ubx_0", "can be updated after creation");
    }

    // Jbx_e
    tmp_idx = malloc(nbx_e*sizeof(int));
    const mxArray *Jbx_e_matlab = mxGetField( matlab_model, 0, "constr_Jbx_e" );
    if (Jbx_e_matlab!=NULL)
    {
        int nrow = (int) mxGetM( Jbx_e_matlab );
        int ncol = (int) mxGetN( Jbx_e_matlab );
        MEX_DIM_CHECK_MAT(fun_name, "constr_Jbx_e", nrow, ncol, nbx_e, nx);
        i_ptr = malloc(nrow*sizeof(int));
        double *Jbx_e = mxGetPr( Jbx_e_matlab );
        for (int ii=0; ii<nrow; ii++)
        {
            int nnz_row = 0;
            for (int jj=0; jj<ncol; jj++)
            {
                if (Jbx_e[ii+nrow*jj]==1.0)
                {
                    tmp_idx[ii] = jj;
                    nnz_row++;
                }
                else if (Jbx_e[ii+nrow*jj]!=0.0)
                {
                    MEX_NONBINARY_MAT(fun_name, "constr_Jbx_e");
                }
            }
            if (nnz_row > 1)
            {
                MEX_MULTIPLE_ONES_IN_ROW(fun_name, "constr_Jbx_e");
            }
        }
        ocp_nlp_constraints_model_set(config, dims, in, N, "idxbx", tmp_idx);
    }
    free(tmp_idx);

    // Jbu
    const mxArray *Jbu_matlab = mxGetField( matlab_model, 0, "constr_Jbu" );
    if (Jbu_matlab != NULL)
    {
        int nrow = (int) mxGetM( Jbu_matlab );
        int ncol = (int) mxGetN( Jbu_matlab );
        MEX_DIM_CHECK_MAT(fun_name, "constr_Jbu", nrow, ncol, nbu, nu);

        double *Jbu = mxGetPr( mxGetField( matlab_model, 0, "constr_Jbu" ) );
        i_ptr = malloc(nrow*sizeof(int));
        for (int ii=0; ii<nrow; ii++)
        {
            int nnz_row = 0;
            for (int jj=0; jj<ncol; jj++)
            {
                if (Jbu[ii+nrow*jj]!=0.0)
                {
                    i_ptr[ii] = jj;
                    nnz_row++;
                }
                else if (Jbu[ii+nrow*jj]!=0.0)
                {
                    MEX_NONBINARY_MAT(fun_name, "constr_Jbu");
                }
            }
            if (nnz_row > 1)
            {
                MEX_MULTIPLE_ONES_IN_ROW(fun_name, "constr_Jbu");
            }
        }
        for (int ii=0; ii<N; ii++)
        {
            ocp_nlp_constraints_model_set(config, dims, in, ii, "idxbu", i_ptr);
        }
        free(i_ptr);
    }
    // TODO: else complain?

    // lbu
    if (mxGetField( matlab_model, 0, "constr_lbu" )!=NULL)
    {
        int matlab_size = (int) mxGetNumberOfElements( mxGetField( matlab_model, 0, "constr_lbu" ) );
        int acados_size = nbu;
        MEX_DIM_CHECK_VEC(fun_name, "constr_lbu", matlab_size, acados_size);

        double *lbu = mxGetPr( mxGetField( matlab_model, 0, "constr_lbu" ) );
        for (int ii=0; ii<N; ii++)
        {
            ocp_nlp_constraints_model_set(config, dims, in, ii, "lbu", lbu);
        }
    }
    // ubu
    if (mxGetField( matlab_model, 0, "constr_ubu" )!=NULL)
    {
        int matlab_size = (int) mxGetNumberOfElements( mxGetField( matlab_model, 0, "constr_ubu" ) );
        int acados_size = nbu;
        MEX_DIM_CHECK_VEC(fun_name, "constr_ubu", matlab_size, acados_size);

        double *ubu = mxGetPr( mxGetField( matlab_model, 0, "constr_ubu" ) );
        for (int ii=0; ii<N; ii++)
        {
            ocp_nlp_constraints_model_set(config, dims, in, ii, "ubu", ubu);
        }
    }


    // C
    const mxArray *C_matlab = mxGetField( matlab_model, 0, "constr_C" );
    if (C_matlab != NULL)
    {
        int nrow = (int) mxGetM( C_matlab );
        int ncol = (int) mxGetN( C_matlab );
        MEX_DIM_CHECK_MAT(fun_name, "constr_C", nrow, ncol, ng, nx);

        double *C = mxGetPr( C_matlab );
        for (int ii=0; ii<N; ii++)
        {
            ocp_nlp_constraints_model_set(config, dims, in, ii, "C", C);
        }
    }
    // D
    const mxArray *D_matlab = mxGetField( matlab_model, 0, "constr_D" );
    if (D_matlab!=NULL)
    {
        int nrow = (int) mxGetM( D_matlab );
        int ncol = (int) mxGetN( D_matlab );
        MEX_DIM_CHECK_MAT(fun_name, "constr_D", nrow, ncol, ng, nu);

        double *D = mxGetPr( D_matlab );
        for (int ii=0; ii<N; ii++)
        {
            ocp_nlp_constraints_model_set(config, dims, in, ii, "D", D);
        }
    }
    // lg
    if (mxGetField( matlab_model, 0, "constr_lg" )!=NULL)
    {
        int matlab_size = (int) mxGetNumberOfElements( mxGetField( matlab_model, 0, "constr_lg" ) );
        int acados_size = ng;
        MEX_DIM_CHECK_VEC(fun_name, "constr_lg", matlab_size, acados_size);
        double *lg = mxGetPr( mxGetField( matlab_model, 0, "constr_lg" ) );
        for (int ii=0; ii<N; ii++)
        {
            ocp_nlp_constraints_model_set(config, dims, in, ii, "lg", lg);
        }
    }
    // ug
    if (mxGetField( matlab_model, 0, "constr_ug" )!=NULL)
    {
        int matlab_size = (int) mxGetNumberOfElements( mxGetField( matlab_model, 0, "constr_ug" ) );
        int acados_size = ng;
        MEX_DIM_CHECK_VEC(fun_name, "constr_ug", matlab_size, acados_size);
        double *ug = mxGetPr( mxGetField( matlab_model, 0, "constr_ug" ) );
        for (int ii=0; ii<N; ii++)
        {
            ocp_nlp_constraints_model_set(config, dims, in, ii, "ug", ug);
        }
    }
    // C_e
    if (mxGetField( matlab_model, 0, "constr_C_e" )!=NULL)
    {
        double *C_e = mxGetPr( mxGetField( matlab_model, 0, "constr_C_e" ) );
        ocp_nlp_constraints_model_set(config, dims, in, N, "C", C_e);
    }
    // lg_e
    if (mxGetField( matlab_model, 0, "constr_lg_e" )!=NULL)
        {
        double *lg_e = mxGetPr( mxGetField( matlab_model, 0, "constr_lg_e" ) );
        ocp_nlp_constraints_model_set(config, dims, in, N, "lg", lg_e);
        }
    // ug_e
    if (mxGetField( matlab_model, 0, "constr_ug_e" )!=NULL)
    {
        double *ug_e = mxGetPr( mxGetField( matlab_model, 0, "constr_ug_e" ) );
        ocp_nlp_constraints_model_set(config, dims, in, N, "ug", ug_e);
    }

    // lh
    if (mxGetField( matlab_model, 0, "constr_lh" )!=NULL)
    {
        int matlab_size = (int) mxGetNumberOfElements( mxGetField( matlab_model, 0, "constr_lh" ) );
        int acados_size = nh;
        MEX_DIM_CHECK_VEC(fun_name, "constr_lh", matlab_size, acados_size);
        double *lh = mxGetPr( mxGetField( matlab_model, 0, "constr_lh" ) );
        for (int ii=0; ii<N; ii++)
        {
            ocp_nlp_constraints_model_set(config, dims, in, ii, "lh", lh);
        }
    }
    // uh
    if (mxGetField( matlab_model, 0, "constr_uh" )!=NULL)
    {
        int matlab_size = (int) mxGetNumberOfElements( mxGetField( matlab_model, 0, "constr_uh" ) );
        int acados_size = nh;
        MEX_DIM_CHECK_VEC(fun_name, "constr_uh", matlab_size, acados_size);
        double *uh = mxGetPr( mxGetField( matlab_model, 0, "constr_uh" ) );
        for (int ii=0; ii<N; ii++)
        {
            ocp_nlp_constraints_model_set(config, dims, in, ii, "uh", uh);
        }
    }
    // lh_e
    if (mxGetField( matlab_model, 0, "constr_lh_e" )!=NULL)
    {
        int matlab_size = (int) mxGetNumberOfElements( mxGetField( matlab_model, 0, "constr_lh_e" ) );
        int acados_size = nh_e;
        MEX_DIM_CHECK_VEC(fun_name, "constr_lh_e", matlab_size, acados_size);
        double *lh_e = mxGetPr( mxGetField( matlab_model, 0, "constr_lh_e" ) );
        ocp_nlp_constraints_model_set(config, dims, in, N, "lh", lh_e);
    }
    // uh_e
    if (mxGetField( matlab_model, 0, "constr_uh_e" )!=NULL)
    {
        int matlab_size = (int) mxGetNumberOfElements( mxGetField( matlab_model, 0, "constr_uh_e" ) );
        int acados_size = nh_e;
        MEX_DIM_CHECK_VEC(fun_name, "constr_uh_e", matlab_size, acados_size);
        double *uh_e = mxGetPr( mxGetField( matlab_model, 0, "constr_uh_e" ) );
        ocp_nlp_constraints_model_set(config, dims, in, N, "uh", uh_e);
    }

    // Jsbu
    const mxArray *Jsbu_matlab = mxGetField( matlab_model, 0, "constr_Jsbu" );
    if (Jsbu_matlab!=NULL)
    {
        int nrow = (int) mxGetM( Jsbu_matlab );
        int ncol = (int) mxGetN( Jsbu_matlab );
        MEX_DIM_CHECK_MAT(fun_name, "constr_Jsbu", nrow, ncol, nbu, nsbu);

        double *Jsbu = mxGetPr( mxGetField( matlab_model, 0, "constr_Jsbu" ) );

        i_ptr = malloc(nsbu*sizeof(int));
        for (int ii=0; ii<nsbu; ii++)
        {
            for (int jj=0; jj<nbu; jj++)
            {
                if (Jsbu[jj+nbu*ii]==1.0)
                {
                    i_ptr[ii] = jj;
                }
                else if (Jsbu[jj+nbu*ii]!=0.0)
                {
                    MEX_NONBINARY_MAT(fun_name, "constr_Jsbu");
                }
            }
        }
        for (int ii=0; ii<N; ii++)
        {
            ocp_nlp_constraints_model_set(config, dims, in, ii, "idxsbu", i_ptr);
        }
        free(i_ptr);
    }
    // Jsbx
    const mxArray *Jsbx_matlab = mxGetField( matlab_model, 0, "constr_Jsbx" );
    if (Jsbx_matlab!=NULL)
    {
        int nrow = (int) mxGetM( Jsbx_matlab );
        int ncol = (int) mxGetN( Jsbx_matlab );
        MEX_DIM_CHECK_MAT(fun_name, "constr_Jsbx", nrow, ncol, nbx, nsbx);

        double *Jsbx = mxGetPr( mxGetField( matlab_model, 0, "constr_Jsbx" ) );
        i_ptr = malloc(nsbx*sizeof(int));
        for (int ii=0; ii<nsbx; ii++)
        {
            for (int jj=0; jj<nbx; jj++)
            {
                if (Jsbx[jj+nbx*ii]==1.0)
                {
                    i_ptr[ii] = jj;
                }
                else if (Jsbx[jj+nbx*ii]!=0.0)
                {
                    MEX_NONBINARY_MAT(fun_name, "constr_Jsbx");
                }
            }
        }
        for (int ii=1; ii<N; ii++) // TODO stage 0 !!!!!!!!!!
        {
            ocp_nlp_constraints_model_set(config, dims, in, ii, "idxsbx", i_ptr);
        }
        free(i_ptr);
    }
    // Jsbx_e
    const mxArray *Jsbx_e_matlab = mxGetField( matlab_model, 0, "constr_Jsbx_e" );
    if (Jsbx_e_matlab!=NULL)
    {
        int nrow = (int) mxGetM( Jsbx_matlab );
        int ncol = (int) mxGetN( Jsbx_matlab );
        MEX_DIM_CHECK_MAT(fun_name, "constr_Jsbx_e", nrow, ncol, nbx_e, nsbx_e);

        double *Jsbx_e = mxGetPr( mxGetField( matlab_model, 0, "constr_Jsbx_e" ) );
        i_ptr = malloc(nsbx_e*sizeof(int));
        for (int ii=0; ii<nsbx_e; ii++)
        {
            for (int jj=0; jj<nbx_e; jj++)
            {
                if (Jsbx_e[jj+nbx_e*ii]==1.0)
                {
                    i_ptr[ii] = jj;
                }
                else if (Jsbx_e[jj+nbx_e*ii]!=0.0)
                {
                    MEX_NONBINARY_MAT(fun_name, "constr_Jsbx_e");
                }
            }
        }
        ocp_nlp_constraints_model_set(config, dims, in, N, "idxsbx", i_ptr);
        free(i_ptr);
    }

    // Jsg
    const mxArray *Jsg_matlab = mxGetField( matlab_model, 0, "constr_Jsg" );
    if (Jsg_matlab!=NULL)
    {
        int nrow = (int) mxGetM( Jsg_matlab );
        int ncol = (int) mxGetN( Jsg_matlab );
        MEX_DIM_CHECK_MAT(fun_name, "constr_Jsg", nrow, ncol, ng, nsg);

        double *Jsg = mxGetPr( mxGetField( matlab_model, 0, "constr_Jsg" ) );
        i_ptr = malloc(nsg*sizeof(int));
        for (int ii=0; ii<nsg; ii++)
        {
            for (int jj=0; jj<ng; jj++)
            {
                if (Jsg[jj+ng*ii]==1.0)
                {
                    i_ptr[ii] = jj;
                }
                else if (Jsg[jj+ng*ii]!=0.0)
                {
                    MEX_NONBINARY_MAT(fun_name, "constr_Jsg");
                }
            }
        }
        for (int ii=0; ii<N; ii++)
        {
            ocp_nlp_constraints_model_set(config, dims, in, ii, "idxsg", i_ptr);
        }
        free(i_ptr);
    }
    // Jsg_e
    const mxArray *Jsg_e_matlab = mxGetField( matlab_model, 0, "constr_Jsg_e" );
    if (Jsg_e_matlab!=NULL)
    {
        int nrow = (int) mxGetM( Jsg_e_matlab );
        int ncol = (int) mxGetN( Jsg_e_matlab );
        MEX_DIM_CHECK_MAT(fun_name, "constr_Jsg_e", nrow, ncol, ng_e, nsg_e);
        
        double *Jsg_e = mxGetPr( mxGetField( matlab_model, 0, "constr_Jsg_e" ) );

        i_ptr = malloc(nsg_e*sizeof(int));
        for (int ii=0; ii<nsg_e; ii++)
        {
            for (int jj=0; jj<ng_e; jj++)
            {
                if (Jsg_e[jj+ng_e*ii]==1.0)
                {
                    i_ptr[ii] = jj;
                }
                else if (Jsg_e[jj+ng_e*ii]!=0.0)
                {
                    MEX_NONBINARY_MAT(fun_name, "constr_Jsg_e");
                }
            }
        }
        ocp_nlp_constraints_model_set(config, dims, in, N, "idxsg", i_ptr);
        free(i_ptr);
    }

    // Jsh
    const mxArray *Jsh_matlab = mxGetField( matlab_model, 0, "constr_Jsh" );
    if (Jsh_matlab!=NULL)
    {
        int nrow = (int) mxGetM( Jsh_matlab );
        int ncol = (int) mxGetN( Jsh_matlab );
        MEX_DIM_CHECK_MAT(fun_name, "constr_Jsh", nrow, ncol, nh, nsh);

        double *Jsh = mxGetPr( mxGetField( matlab_model, 0, "constr_Jsh" ) );
        i_ptr = malloc(nsh*sizeof(int));
        for (int ii=0; ii<nsh; ii++)
        {
            for (int jj=0; jj<nh; jj++)
            {
                if (Jsh[jj+nh*ii]==1.0)
                {
                    i_ptr[ii] = jj;
                }
                else if (Jsh[jj+nh*ii]!=0.0)
                {
                    MEX_NONBINARY_MAT(fun_name, "constr_Jsh");
                }
            }
        }
        for (int ii=0; ii<N; ii++)
        {
            ocp_nlp_constraints_model_set(config, dims, in, ii, "idxsh", i_ptr);
        }
        free(i_ptr);
    }

    // Jsh_e
    const mxArray *Jsh_e_matlab = mxGetField( matlab_model, 0, "constr_Jsh_e" );
    if (Jsh_e_matlab!=NULL)
    {
        int nrow = (int) mxGetM( Jsh_e_matlab );
        int ncol = (int) mxGetN( Jsh_e_matlab );
        MEX_DIM_CHECK_MAT(fun_name, "constr_Jsh_e", nrow, ncol, nh_e, nsh_e);
        
        double *Jsh_e = mxGetPr( mxGetField( matlab_model, 0, "constr_Jsh_e" ) );

        i_ptr = malloc(nsh_e*sizeof(int));
        for (int ii=0; ii<nsh_e; ii++)
        {
            for (int jj=0; jj<nh_e; jj++)
            {
                if (Jsh_e[jj+nh_e*ii]==1.0)
                {
                    i_ptr[ii] = jj;
                }
                else if (Jsh_e[jj+nh_e*ii]!=0.0)
                {
                    MEX_NONBINARY_MAT(fun_name, "constr_Jsh_e");
                }
            }
        }
        ocp_nlp_constraints_model_set(config, dims, in, N, "idxsh", i_ptr);
        free(i_ptr);
    }
    // mexPrintf("\nocp_create end\n");

    // idxbxe_0
    tmp_idx = malloc(nbxe_0*sizeof(int));

    double *idxbxe_0;
    const mxArray *idxbxe_0_matlab = mxGetField( matlab_model, 0, "constr_idxbxe_0" );
    if (idxbxe_0_matlab!=NULL)
    {
        int matlab_size = (int) mxGetNumberOfElements( idxbxe_0_matlab );
        MEX_DIM_CHECK_VEC(fun_name, "constr_idxbxe_0", matlab_size, nbxe_0);
        idxbxe_0 = mxGetPr( idxbxe_0_matlab );
        for (int ii=0; ii<nbxe_0; ii++)
        {
			tmp_idx[ii] = (int) idxbxe_0[ii];
        }
        ocp_nlp_constraints_model_set(config, dims, in, 0, "idxbxe", tmp_idx);
    }
    free(tmp_idx);


    return;

}
// TODO: maybe add the following: lower bounds for slack variables
//    double *lsbu;    bool set_lsbu = false;
//    double *usbu;    bool set_usbu = false;
//    double *lsbx;    bool set_lsbx = false;
//    double *usbx;    bool set_usbx = false;
//    double *lsg;    bool set_lsg = false;
//    double *usg;    bool set_usg = false;
//    double *lsg_e;    bool set_lsg_e = false;
//    double *usg_e;    bool set_usg_e = false;
//    double *lsh;    bool set_lsh = false;
//    double *ush;    bool set_ush = false;
//    double *lsh_e;    bool set_lsh_e = false;
//    double *ush_e;    bool set_ush_e = false;

    // lsbu
//    if (mxGetField( matlab_model, 0, "constr_lsbu" )!=NULL)
//        {
//        set_lsbu = true;
//        lsbu = mxGetPr( mxGetField( matlab_model, 0, "constr_lsbu" ) );
//        }
    // usbu
//    if (mxGetField( matlab_model, 0, "constr_usbu" )!=NULL)
//        {
//        set_usbu = true;
//        usbu = mxGetPr( mxGetField( matlab_model, 0, "constr_usbu" ) );
//        }
    // lsbx
//    if (mxGetField( matlab_model, 0, "constr_lsbx" )!=NULL)
//        {
//        set_lsbx = true;
//        lsbx = mxGetPr( mxGetField( matlab_model, 0, "constr_lsbx" ) );
//        }
    // usbx
//    if (mxGetField( matlab_model, 0, "constr_usbx" )!=NULL)
//        {
//        set_usbx = true;
//        usbx = mxGetPr( mxGetField( matlab_model, 0, "constr_usbx" ) );
//        }

    // lsg
//    if (mxGetField( matlab_model, 0, "constr_lsg" )!=NULL)
//        {
//        set_lsg = true;
//        lsg = mxGetPr( mxGetField( matlab_model, 0, "constr_lsg" ) );
//        }
    // usg
//    if (mxGetField( matlab_model, 0, "constr_usg" )!=NULL)
//        {
//        set_usg = true;
//        usg = mxGetPr( mxGetField( matlab_model, 0, "constr_usg" ) );
//        }
    // lsg_e
//    if (mxGetField( matlab_model, 0, "constr_lsg_e" )!=NULL)
//        {
//        set_lsg_e = true;
//        lsg_e = mxGetPr( mxGetField( matlab_model, 0, "constr_lsg_e" ) );
//        }
    // usg_e
//    if (mxGetField( matlab_model, 0, "constr_usg_e" )!=NULL)
//        {
//        set_usg_e = true;
//        usg_e = mxGetPr( mxGetField( matlab_model, 0, "constr_usg_e" ) );
//        }
    // lsh
//    if (mxGetField( matlab_model, 0, "constr_lsh" )!=NULL)
//        {
//        set_lsh = true;
//        lsh = mxGetPr( mxGetField( matlab_model, 0, "constr_lsh" ) );
//        }
    // ush
//    if (mxGetField( matlab_model, 0, "constr_ush" )!=NULL)
//        {
//        set_ush = true;
//        ush = mxGetPr( mxGetField( matlab_model, 0, "constr_ush" ) );
//        }
    // lsh_e
//    if (mxGetField( matlab_model, 0, "constr_lsh_e" )!=NULL)
//        {
//        set_lsh_e = true;
//        lsh_e = mxGetPr( mxGetField( matlab_model, 0, "constr_lsh_e" ) );
//        }
    // ush_e
//    if (mxGetField( matlab_model, 0, "constr_ush_e" )!=NULL)
//        {
//        set_ush_e = true;
//        ush_e = mxGetPr( mxGetField( matlab_model, 0, "constr_ush_e" ) );
//        }
//    if (set_lsbu)
//        {
//        for (int ii=0; ii<N; ii++)
//            {
//            ocp_nlp_constraints_model_set(config, dims, in, ii, "lsbu", lsbu);
//            }
//        }
//    if (set_usbu)
//        {
//        for (int ii=0; ii<N; ii++)
//            {
//            ocp_nlp_constraints_model_set(config, dims, in, ii, "usbu", usbu);
//            }
//        }
//    if (set_lsbx)
//        {
//        for (int ii=0; ii<=N; ii++)
//        for (int ii=1; ii<=N; ii++) // TODO stage 0 !!!!!!!!!!
//            {
//            ocp_nlp_constraints_model_set(config, dims, in, ii, "lsbx", lsbx);
//            }
//        }
//    if (set_usbx)
//        {
//        for (int ii=0; ii<=N; ii++)
//        for (int ii=1; ii<=N; ii++) // TODO stage 0 !!!!!!!!!!
//            {
//            ocp_nlp_constraints_model_set(config, dims, in, ii, "usbx", usbx);
//            }
//        }
//    if (set_lsg)
//        {
//        for (int ii=0; ii<N; ii++)
//            {
//            ocp_nlp_constraints_model_set(config, dims, in, ii, "lsg", lsg);
//            }
//        }
//    if (set_usg)
//        {
//        for (int ii=0; ii<N; ii++)
//            {
//            ocp_nlp_constraints_model_set(config, dims, in, ii, "usg", usg);
//            }
//        }
//    if (set_lsg_e)
//        {
//        ocp_nlp_constraints_model_set(config, dims, in, N, "lsg", lsg_e);
//        }
//    if (set_usg_e)
//        {
//        ocp_nlp_constraints_model_set(config, dims, in, N, "usg", usg_e);
//        }
//    if (set_lsh)
//        {
//        for (int ii=0; ii<N; ii++)
//            {
//            ocp_nlp_constraints_model_set(config, dims, in, ii, "lsh", lsh);
//            }
//        }
//    if (set_ush)
//        {
//        for (int ii=0; ii<N; ii++)
//            {
//            ocp_nlp_constraints_model_set(config, dims, in, ii, "ush", ush);
//            }
//        }
//    if (set_lsh_e)
//        {
//        ocp_nlp_constraints_model_set(config, dims, in, N, "lsh", lsh_e);
//        }
//    if (set_ush_e)
//        {
//        ocp_nlp_constraints_model_set(config, dims, in, N, "ush", ush_e);
//        }
