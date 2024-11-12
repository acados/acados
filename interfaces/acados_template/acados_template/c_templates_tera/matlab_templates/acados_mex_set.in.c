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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// acados
#include "acados/utils/print.h"
#include "acados_c/ocp_nlp_interface.h"
#include "acados_solver_{{ name }}.h"

// mex
#include "mex.h"
#include "mex_macros.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    long long *ptr;
    int acados_size;
    mxArray *mex_field;
    char fun_name[20] = "ocp_set";
    char buffer [500]; // for error messages

    /* RHS */
    int min_nrhs = 3;

    // C ocp
    const mxArray *C_ocp = prhs[0];
    // capsule
    ptr = (long long *) mxGetData( mxGetField( C_ocp, 0, "capsule" ) );
    {{ name }}_solver_capsule *capsule = ({{ name }}_solver_capsule *) ptr[0];
    // plan
    ptr = (long long *) mxGetData( mxGetField( C_ocp, 0, "plan" ) );
    ocp_nlp_plan_t *plan = (ocp_nlp_plan_t *) ptr[0];
    // config
    ptr = (long long *) mxGetData( mxGetField( C_ocp, 0, "config" ) );
    ocp_nlp_config *config = (ocp_nlp_config *) ptr[0];
    // dims
    ptr = (long long *) mxGetData( mxGetField( C_ocp, 0, "dims" ) );
    ocp_nlp_dims *dims = (ocp_nlp_dims *) ptr[0];
    // opts
    ptr = (long long *) mxGetData( mxGetField( C_ocp, 0, "opts" ) );
    void *opts = (void *) ptr[0];
    // in
    ptr = (long long *) mxGetData( mxGetField( C_ocp, 0, "in" ) );
    ocp_nlp_in *in = (ocp_nlp_in *) ptr[0];
    // out
    ptr = (long long *) mxGetData( mxGetField( C_ocp, 0, "out" ) );
    ocp_nlp_out *out = (ocp_nlp_out *) ptr[0];
    // solver
    ptr = (long long *) mxGetData( mxGetField( C_ocp, 0, "solver" ) );
    ocp_nlp_solver *solver = (ocp_nlp_solver *) ptr[0];

    // field
    char *field = mxArrayToString( prhs[1] );
    // value
    double *value = mxGetPr( prhs[2] );

    // for checks
    int matlab_size = (int) mxGetNumberOfElements( prhs[2] );
    int nrow = (int) mxGetM( prhs[2] );
    int ncol = (int) mxGetN( prhs[2] );

    int N = dims->N;
    int tmp_int, offset;

    // stage
    int s0, se;
    if (nrhs == min_nrhs)
    {
        s0 = 0;
        se = N;
    }
    else if (nrhs == min_nrhs+1)
    {
        s0 = mxGetScalar( prhs[3] );
        if (s0 > N)
        {
            sprintf(buffer, "ocp_set: N < specified stage = %d\n", s0);
            mexErrMsgTxt(buffer);
        }
        se = s0 + 1;
    }
    else
    {
        sprintf(buffer, "ocp_set: wrong nrhs: %d\n", nrhs);
        mexErrMsgTxt(buffer);
    }

    /* Set value */
    // constraints
    if (!strcmp(field, "constr_x0"))
    {
        int nbx = ocp_nlp_dims_get_from_attr(config, dims, out, 0, "lbx");
        acados_size = nbx;
        MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
        ocp_nlp_constraints_model_set(config, dims, in, 0, "lbx", value);
        ocp_nlp_constraints_model_set(config, dims, in, 0, "ubx", value);
    }
    else if (!strcmp(field, "constr_C"))
    {
        for (int ii=s0; ii<se; ii++)
        {
            int ng = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "ug");
            int nx = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "x");
            MEX_DIM_CHECK_MAT(fun_name, "constr_C", nrow, ncol, ng, nx);
            if (matlab_size != 0)
                ocp_nlp_constraints_model_set(config, dims, in, ii, "C", value);
        }
    }
    else if (!strcmp(field, "constr_lbx"))
    {
        for (int ii=s0; ii<se; ii++)
        {
            acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "lbx");
            MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
            if (matlab_size != 0)
                ocp_nlp_constraints_model_set(config, dims, in, ii, "lbx", value);
        }
    }
    else if (!strcmp(field, "constr_ubx"))
    {
        for (int ii=s0; ii<se; ii++)
        {
            acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "ubx");
            MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
            if (matlab_size != 0)
                ocp_nlp_constraints_model_set(config, dims, in, ii, "ubx", value);
        }
    }
    else if (!strcmp(field, "constr_lbu"))
    {
        for (int ii=s0; ii<se; ii++)
        {
            acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "lbu");
            MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
            if (matlab_size != 0)
                ocp_nlp_constraints_model_set(config, dims, in, ii, "lbu", value);
        }
    }
    else if (!strcmp(field, "constr_ubu"))
    {
        for (int ii=s0; ii<se; ii++)
        {
            acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "ubu");
            MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
            if (matlab_size != 0)
                ocp_nlp_constraints_model_set(config, dims, in, ii, "ubu", value);
        }
    }
    else if (!strcmp(field, "constr_D"))
    {
        for (int ii=s0; ii<se; ii++)
        {
            int ng = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "ug");
            int nu = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "u");
            MEX_DIM_CHECK_MAT(fun_name, "constr_D", nrow, ncol, ng, nu);
            if (matlab_size != 0)
                ocp_nlp_constraints_model_set(config, dims, in, ii, "D", value);
        }
    }
    else if (!strcmp(field, "constr_lg"))
    {
        for (int ii=s0; ii<se; ii++)
        {
            acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "lg");
            MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
            if (matlab_size != 0)
                ocp_nlp_constraints_model_set(config, dims, in, ii, "lg", value);
        }
    }
    else if (!strcmp(field, "constr_ug"))
    {
        for (int ii=s0; ii<se; ii++)
        {
            acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "ug");
            MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
            if (matlab_size != 0)
                ocp_nlp_constraints_model_set(config, dims, in, ii, "ug", value);
        }
    }
    else if (!strcmp(field, "constr_lh"))
    {
        for (int ii=s0; ii<se; ii++)
        {
            acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "lh");
            MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
            if (matlab_size != 0)
                ocp_nlp_constraints_model_set(config, dims, in, ii, "lh", value);
        }
    }
    else if (!strcmp(field, "constr_uh"))
    {
        for (int ii=s0; ii<se; ii++)
        {
            acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "uh");
            MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
            if (matlab_size != 0)
                ocp_nlp_constraints_model_set(config, dims, in, ii, "uh", value);
        }
    }
    // cost:
    else if (!strcmp(field, "cost_y_ref"))
    {
        for (int ii=s0; ii<se; ii++)
        {
            if ((plan->nlp_cost[ii] == LINEAR_LS) || (plan->nlp_cost[ii] == NONLINEAR_LS))
            {
                acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "y_ref");
                MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
                if (matlab_size != 0)
                    ocp_nlp_cost_model_set(config, dims, in, ii, "y_ref", value);
            }
            else
            {
                MEX_FIELD_NOT_SUPPORTED_FOR_COST_STAGE(fun_name, field, plan->nlp_cost[ii], ii);
            }
        }
    }
    else if (!strcmp(field, "cost_y_ref_e"))
    {
        acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, N, "y_ref");
        MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
        if (matlab_size != 0)
            ocp_nlp_cost_model_set(config, dims, in, N, "y_ref", value);
    }
    else if (!strcmp(field, "cost_Vu"))
    {
        for (int ii=s0; ii<se; ii++)
        {
            if ((plan->nlp_cost[ii] == LINEAR_LS) || (plan->nlp_cost[ii] == NONLINEAR_LS))
            {
                int ny = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "y_ref");
                int nu = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "u");
                acados_size = ny * nu;
                MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
                if (matlab_size != 0)
                    ocp_nlp_cost_model_set(config, dims, in, ii, "Vu", value);
            }
            else
            {
                MEX_FIELD_NOT_SUPPORTED_FOR_COST_STAGE(fun_name, field, plan->nlp_cost[ii], ii);
            }
        }
    }
    else if (!strcmp(field, "cost_Vx"))
    {
        for (int ii=s0; ii<se; ii++)
        {
            if ((plan->nlp_cost[ii] == LINEAR_LS) || (plan->nlp_cost[ii] == NONLINEAR_LS))
            {
                int ny = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "y_ref");
                int nx = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "x");
                acados_size = ny * nx;
                MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
                if (matlab_size != 0)
                    ocp_nlp_cost_model_set(config, dims, in, ii, "Vx", value);
            }
            else
            {
                MEX_FIELD_NOT_SUPPORTED_FOR_COST_STAGE(fun_name, field, plan->nlp_cost[ii], ii);
            }
        }
    }
    else if (!strcmp(field, "cost_W"))
    {
        for (int ii=s0; ii<se; ii++)
        {
            if ((plan->nlp_cost[ii] == LINEAR_LS) || (plan->nlp_cost[ii] == NONLINEAR_LS))
            {
                int ny = ocp_nlp_dims_get_from_attr(config, dims, out, s0, "y_ref");
                acados_size = ny * ny;
                MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
                if (matlab_size != 0)
                    ocp_nlp_cost_model_set(config, dims, in, ii, "W", value);
            }
            else
            {
                MEX_FIELD_NOT_SUPPORTED_FOR_COST_STAGE(fun_name, field, plan->nlp_cost[ii], ii);
            }
        }
    }
    else if (!strcmp(field, "cost_Z"))
    {
        for (int ii=s0; ii<se; ii++)
        {
            acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, s0, "cost_Z");
            MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
            if (matlab_size != 0)
                ocp_nlp_cost_model_set(config, dims, in, ii, "Z", value);
        }
    }
    else if (!strcmp(field, "cost_Zl"))
    {
        for (int ii=s0; ii<se; ii++)
        {
            acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, s0, "Zl");
            MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
            if (matlab_size != 0)
                ocp_nlp_cost_model_set(config, dims, in, ii, "Zl", value);
        }
    }
    else if (!strcmp(field, "cost_Zu"))
    {
        for (int ii=s0; ii<se; ii++)
        {
            acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, s0, "Zu");
            MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
            if (matlab_size != 0)
                ocp_nlp_cost_model_set(config, dims, in, ii, "Zu", value);
        }
    }
    else if (!strcmp(field, "cost_z"))
    {
        for (int ii=s0; ii<se; ii++)
        {
            acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, s0, "cost_z");
            MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
            if (matlab_size != 0)
                ocp_nlp_cost_model_set(config, dims, in, ii, "z", value);
        }
    }
    else if (!strcmp(field, "cost_zl"))
    {
        for (int ii=s0; ii<se; ii++)
        {
            acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, s0, "zl");
            MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
            if (matlab_size != 0)
                ocp_nlp_cost_model_set(config, dims, in, ii, "zl", value);
        }
    }
    else if (!strcmp(field, "cost_zu"))
    {
        for (int ii=s0; ii<se; ii++)
        {
            acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, s0, "zu");
            MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
            if (matlab_size != 0)
                ocp_nlp_cost_model_set(config, dims, in, ii, "zu", value);
        }
    }
    // initializations
    else if (!strcmp(field, "init_x") || !strcmp(field, "x"))
    {
        if (nrhs == min_nrhs)
        {
            acados_size = 0;
            for (int ii=0; ii<=N; ii++)
            {
                tmp_int = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "x");
                acados_size += tmp_int;
            }
            MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
            offset = 0;
            for (int ii=0; ii<=N; ii++)
            {
                ocp_nlp_out_set(config, dims, out, ii, "x", value+offset);
                tmp_int = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "x");
                offset += tmp_int;
            }
        }
        else // (nrhs == min_nrhs + 1)
        {
            acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, s0, "x");
            MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
            ocp_nlp_out_set(config, dims, out, s0, "x", value);
        }
    }
    else if (!strcmp(field, "init_u") || !strcmp(field, "u"))
    {
        if (nrhs == min_nrhs)
        {
            acados_size = 0;
            for (int ii=0; ii<N; ii++)
            {
                tmp_int = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "u");
                acados_size += tmp_int;
            }
            MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
            offset = 0;
            for (int ii=0; ii<N; ii++)
            {
                ocp_nlp_out_set(config, dims, out, ii, "u", value+offset);
                tmp_int = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "u");
                offset += tmp_int;
            }
        }
        else // (nrhs == min_nrhs + 1)
        {
            acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, s0, "u");
            MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
            ocp_nlp_out_set(config, dims, out, s0, "u", value);
        }
    }
    else if (!strcmp(field, "init_z")||!strcmp(field, "z"))
    {
{% if problem_class == "MOCP" %}
        MEX_FIELD_NOT_SUPPORTED(fun_name, field);
{% else %}
        sim_solver_plan_t sim_plan = plan->sim_solver_plan[0];
        sim_solver_t type = sim_plan.sim_solver;
        if (type == IRK)
        {
            int nz = ocp_nlp_dims_get_from_attr(config, dims, out, 0, "z");
            if (nrhs == min_nrhs)
            {
                acados_size = N*nz;
                MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
                for (int ii=0; ii<N; ii++)
                {
                    ocp_nlp_set(solver, ii, "z_guess", value+ii*nz);
                }
            }
            else // (nrhs == min_nrhs+1)
            {
                acados_size = nz;
                MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
                ocp_nlp_set(solver, s0, "z_guess", value);
            }
        }
        else
        {
            MEX_FIELD_ONLY_SUPPORTED_FOR_SOLVER(fun_name, "init_z", "irk")
        }
{% endif %}
    }
    else if (!strcmp(field, "init_xdot")||!strcmp(field, "xdot"))
    {
{% if problem_class == "MOCP" %}
        MEX_FIELD_NOT_SUPPORTED(fun_name, field);
{% else %}
        sim_solver_plan_t sim_plan = plan->sim_solver_plan[0];
        sim_solver_t type = sim_plan.sim_solver;
        if (type == IRK)
        {
            int nx = ocp_nlp_dims_get_from_attr(config, dims, out, 0, "x");
            if (nrhs == min_nrhs)
            {
                acados_size = N*nx;
                MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
                for (int ii=0; ii<N; ii++)
                {
                    ocp_nlp_set(solver, ii, "xdot_guess", value+ii*nx);
                }
            }
            else // nrhs == min_nrhs+1)
            {
                acados_size = nx;
                MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
                ocp_nlp_set(solver, s0, "xdot_guess", value);
            }
        }
        else
        {
            MEX_FIELD_ONLY_SUPPORTED_FOR_SOLVER(fun_name, "init_z", "irk")
        }
{% endif %}

    }
    else if (!strcmp(field, "init_gnsf_phi")||!strcmp(field, "gnsf_phi"))
    {
{% if problem_class == "MOCP" %}
        MEX_FIELD_NOT_SUPPORTED(fun_name, field);
{% else %}
        sim_solver_plan_t sim_plan = plan->sim_solver_plan[0];
        sim_solver_t type = sim_plan.sim_solver;
        if (type == GNSF)
        {
            int nout = ocp_nlp_dims_get_from_attr(config, dims, out, 0, "init_gnsf_phi");

            if (nrhs == min_nrhs)
            {
                acados_size = N*nout;
                MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
                for (int ii=0; ii<N; ii++)
                {
                    ocp_nlp_set(solver, ii, "gnsf_phi_guess", value+ii*nout);
                }
            }
            else // (nrhs == min_nrhs+1)
            {
                acados_size = nout;
                MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
                ocp_nlp_set(solver, s0, "gnsf_phi_guess", value);
            }
        }
        else
        {
            MEX_FIELD_ONLY_SUPPORTED_FOR_SOLVER(fun_name, "init_gnsf_phi", "irk_gnsf")
        }
{% endif %}
    }
    else if (!strcmp(field, "init_pi")||!strcmp(field, "pi"))
    {
        if (nrhs == min_nrhs)
        {
            acados_size = 0;
            for (int ii=0; ii<N; ii++)
            {
                tmp_int = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "pi");
                acados_size += tmp_int;
            }
            MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
            offset = 0;
            for (int ii=0; ii<N; ii++)
            {
                ocp_nlp_out_set(config, dims, out, ii, "pi", value+offset);
                tmp_int = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "pi");
                offset += tmp_int;
            }
        }
        else // (nrhs == min_nrhs + 1)
        {
            acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, s0, "pi");
            MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
            ocp_nlp_out_set(config, dims, out, s0, "pi", value);
        }
    }
    else if (!strcmp(field, "init_lam")||!strcmp(field, "lam"))
    {
        if (nrhs == min_nrhs)
        {
            MEX_SETTER_NO_ALL_STAGES_SUPPORT(fun_name, field)
        }
        else //(nrhs == min_nrhs+1)
        {
            acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, s0, "lam");
            MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
            ocp_nlp_out_set(config, dims, out, s0, "lam", value);
        }
    }
    else if (!strcmp(field, "init_sl")||!strcmp(field, "sl"))
    {
        if (nrhs == min_nrhs)
        {
            MEX_SETTER_NO_ALL_STAGES_SUPPORT(fun_name, field)
        }
        else //(nrhs == min_nrhs+1)
        {
            acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, s0, "sl");
            MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
            ocp_nlp_out_set(config, dims, out, s0, field, value);
        }
    }
    else if (!strcmp(field, "init_su")||!strcmp(field, "su"))
    {
        if (nrhs == min_nrhs)
        {
            MEX_SETTER_NO_ALL_STAGES_SUPPORT(fun_name, field)
        }
        else //(nrhs == min_nrhs+1)
        {
            acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, s0, "su");
            MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
            ocp_nlp_out_set(config, dims, out, s0, field, value);
        }
    }
    else if (!strcmp(field, "p"))
    {
        if (nrhs == min_nrhs) // all stages
        {
            for (int ii=0; ii<=N; ii++)
            {
                acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "p");
                MEX_DIM_CHECK_VEC_STAGE(fun_name, field, ii, matlab_size, acados_size)
                {{ name }}_acados_update_params(capsule, ii, value, matlab_size);
            }
        }
        else if (nrhs == min_nrhs+1) // one stage
        {
            int stage = mxGetScalar( prhs[3] );
            {{ name }}_acados_update_params(capsule, stage, value, matlab_size);
        }
    }
    else if (!strcmp(field, "p_global"))
    {
        if (nrhs == min_nrhs)
        {
            {{ name }}_acados_set_p_global_and_precompute_dependencies(capsule, value, matlab_size);
        }
        else if (nrhs > min_nrhs)
        {
            sprintf(buffer, "ocp_set: p_global cannot be set stage-wise.");
            mexErrMsgTxt(buffer);
        }
    }
    else if (!strcmp(field, "params_sparse"))
    {
{%- if problem_class == "MOCP" %}
    {%- set np_values = [] -%}
    {%- for jj in range(end=n_phases) %}
        {%- set_global np_values = np_values | concat(with=(phases_dims[jj].np)) %}
    {%- endfor %}
    {%- set np_max = np_values | sort | last %}
{% else %}
    {%- set np_max = dims.np %}
{%- endif %}

{% if np_max > 0 %}
        MEX_DIM_CHECK_MAT(fun_name, field, nrow, ncol, nrow, 2);
        // create int index vector
        int idx_tmp[{{ np_max }}];
        for (int ip = 0; ip<nrow; ip++)
            idx_tmp[ip] = (int) value[ip];

        if (nrhs == min_nrhs) // all stages
        {
            for (int ii=0; ii<=N; ii++)
            {
                {{ name }}_acados_update_params_sparse(capsule, ii, idx_tmp, value+nrow, nrow);
            }
        }
        else if (nrhs == min_nrhs+1) // one stage
        {
            int stage = mxGetScalar( prhs[3] );
                {{ name }}_acados_update_params_sparse(capsule, stage, idx_tmp, value+nrow, nrow);
        }
{% endif %}
    }
/* OPTIONS */
    else if (!strcmp(field, "nlp_solver_max_iter"))
    {
        acados_size = 1;
        MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
        int nlp_solver_max_iter = (int) value[0];
        ocp_nlp_solver_opts_set(config, opts, "max_iter", &nlp_solver_max_iter);
    }
    else if (!strcmp(field, "rti_phase"))
    {
        acados_size = 1;
        MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
        int rti_phase = (int) value[0];
        if (plan->nlp_solver == SQP)
        {
            MEX_FIELD_ONLY_SUPPORTED_FOR_SOLVER(fun_name, field, "sqp_rti")
        }
        ocp_nlp_solver_opts_set(config, opts, "rti_phase", &rti_phase);
    }
    else if (!strcmp(field, "qp_warm_start"))
    {
        acados_size = 1;
        MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
        int qp_warm_start = (int) value[0];
        ocp_nlp_solver_opts_set(config, opts, "qp_warm_start", &qp_warm_start);
    }
    else if (!strcmp(field, "qp_mu0") || !strcmp(field, "qp_solver_mu0"))
    {
        if (!(plan->ocp_qp_solver_plan.qp_solver == FULL_CONDENSING_HPIPM || plan->ocp_qp_solver_plan.qp_solver == PARTIAL_CONDENSING_HPIPM))
        {
            MEX_FIELD_ONLY_SUPPORTED_FOR_SOLVER(fun_name, "qp_mu0", "HPIPM")
        }
        acados_size = 1;
        MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
        double qp_mu0 = (double) value[0];
        ocp_nlp_solver_opts_set(config, opts, "qp_mu0", &qp_mu0);
    }
    else if (!strcmp(field, "qp_print_level"))
    {
        if (!(plan->ocp_qp_solver_plan.qp_solver == FULL_CONDENSING_HPIPM || plan->ocp_qp_solver_plan.qp_solver == PARTIAL_CONDENSING_HPIPM))
        {
            MEX_FIELD_ONLY_SUPPORTED_FOR_SOLVER(fun_name, "qp_print_level", "HPIPM")
        }
        acados_size = 1;
        MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
        int qp_print_level = (int) value[0];
        ocp_nlp_solver_opts_set(config, opts, "qp_print_level", &qp_print_level);
    }
    else if (!strcmp(field, "warm_start_first_qp"))
    {
        acados_size = 1;
        MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
        int warm_start_first_qp = (int) value[0];
        ocp_nlp_solver_opts_set(config, opts, "warm_start_first_qp", &warm_start_first_qp);
    }
    else if (!strcmp(field, "print_level"))
    {
        acados_size = 1;
        MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
        int print_level = (int) value[0];
        ocp_nlp_solver_opts_set(config, opts, "print_level", &print_level);
    }
    else
    {
        MEX_FIELD_NOT_SUPPORTED_SUGGEST(fun_name, field, "p, constr_x0,\
 constr_lbx, constr_ubx, constr_C, constr_D, constr_lg, constr_ug, constr_lh, constr_uh,\
 constr_lbu, constr_ubu, cost_y_ref[_e], sl, su, x, xdot, u, pi, lam, z, \
 cost_Vu, cost_Vx, cost_Vz, cost_W, cost_Z, cost_Zl, cost_Zu, cost_z,\
 cost_zl, cost_zu, init_x, init_u, init_z, init_xdot, init_gnsf_phi,\
 init_pi, nlp_solver_max_iter, qp_warm_start, qp_solver_mu0, qp_print_level, warm_start_first_qp, print_level");
    }

    return;
}

