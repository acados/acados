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
#include "acados/utils/strsep.h"
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

    char *ptr_field_name = NULL;
    int field_name_length = 0;
    char field_name[128];

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
    int ii;

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
        for (ii=s0; ii<se; ii++)
        {
            int ng = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "ug");
            int nx = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "x");
            MEX_DIM_CHECK_MAT(fun_name, "constr_C", nrow, ncol, ng, nx);
            if (matlab_size != 0)
                ocp_nlp_constraints_model_set(config, dims, in, ii, "C", value);
        }
    }
    else if (!strcmp(field, "constr_D"))
    {
        for (ii=s0; ii<se; ii++)
        {
            int ng = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "ug");
            int nu = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "u");
            MEX_DIM_CHECK_MAT(fun_name, "constr_D", nrow, ncol, ng, nu);
            if (matlab_size != 0)
                ocp_nlp_constraints_model_set(config, dims, in, ii, "D", value);
        }
    }
    else if (!strcmp(field, "constr_lbx") || !strcmp(field, "constr_ubx") ||
             !strcmp(field, "constr_lh") || !strcmp(field, "constr_uh") ||
             !strcmp(field, "constr_lg") || !strcmp(field, "constr_ug") ||
             !strcmp(field, "constr_lbu") || !strcmp(field, "constr_ubu"))
    {
        extract_field_name(field, field_name, &field_name_length, &ptr_field_name);

        if (nrhs == min_nrhs) // all stages
        {
            acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, 0, field_name);

            if (acados_size == matlab_size) // set the same value for all stages for which the dimension is not 0
            {
                for (ii=0; ii<=N; ii++)
                {
                    acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, ii, field_name);
                    if (matlab_size != 0)
                    {
                        // NOTE: checking only stages with dimension > 0 allows this to work
                        // also for lbu/ubu. for which dimension at terminal stage is 0
                        MEX_DIM_CHECK_VEC_STAGE(fun_name, field, ii, matlab_size, acados_size)
                        ocp_nlp_constraints_model_set(config, dims, in, ii, field_name, value);
                    }
                }
            }
            else
            {
                acados_size = ocp_nlp_dims_get_total_from_attr(config, dims, out, field_name);
                MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
                offset = 0;
                for (ii=0; ii<=N; ii++) // TODO implement set_all
                {
                    ocp_nlp_constraints_model_set(config, dims, in, ii, field_name, value+offset);
                    tmp_int = ocp_nlp_dims_get_from_attr(config, dims, out, ii, field_name);
                    offset += tmp_int;
                }
            }
        }
        else if (nrhs == min_nrhs + 1) // single stage
        {
            acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, s0, field_name);
            MEX_DIM_CHECK_VEC_STAGE(fun_name, field, s0, matlab_size, acados_size)
            if (matlab_size != 0)
                ocp_nlp_constraints_model_set(config, dims, in, s0, field_name, value);
        }
    }
    // cost:
    else if (!strcmp(field, "cost_y_ref"))
    {
        for (ii=s0; ii<se; ii++)
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
        for (ii=s0; ii<se; ii++)
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
        for (ii=s0; ii<se; ii++)
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
        for (ii=s0; ii<se; ii++)
        {
            if ((plan->nlp_cost[ii] == LINEAR_LS) || (plan->nlp_cost[ii] == NONLINEAR_LS))
            {
                int ny = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "y_ref");
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
    else if (!strcmp(field, "cost_z") || !strcmp(field, "cost_Z")) // NOTE this cannot be merged with next case due to dimension getter
    {
        extract_field_name(field, field_name, &field_name_length, &ptr_field_name);

        if (nrhs == min_nrhs) // all stages
        {
            acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, 0, field);

            if (acados_size == matlab_size) // set the same value for all stages
            {
                for (ii=0; ii<=N; ii++)
                {
                    acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, ii, field);
                    MEX_DIM_CHECK_VEC_STAGE(fun_name, field, ii, matlab_size, acados_size)
                    if (matlab_size != 0)
                        ocp_nlp_cost_model_set(config, dims, in, ii, field_name, value);
                }
            }
            else
            {
                acados_size = ocp_nlp_dims_get_total_from_attr(config, dims, out, field);
                MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
                offset = 0;
                for (ii=0; ii<=N; ii++) // TODO implement set_all
                {
                    ocp_nlp_cost_model_set(config, dims, in, ii, field_name, value+offset);
                    tmp_int = ocp_nlp_dims_get_from_attr(config, dims, out, ii, field);
                    offset += tmp_int;
                }
            }
        }
        else if (nrhs == min_nrhs + 1) // single stage
        {
            acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, s0, field);
            MEX_DIM_CHECK_VEC_STAGE(fun_name, field, s0, matlab_size, acados_size)
            if (matlab_size != 0)
                ocp_nlp_cost_model_set(config, dims, in, s0, field_name, value);
        }
    }
    else if (!strcmp(field, "cost_zl") || !strcmp(field, "cost_zu") ||
             !strcmp(field, "cost_Zl") || !strcmp(field, "cost_Zu"))
    {

        extract_field_name(field, field_name, &field_name_length, &ptr_field_name);

        if (nrhs == min_nrhs) // all stages
        {
            acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, 0, field_name);

            if (acados_size == matlab_size) // set the same value for all stages
            {
                for (ii=0; ii<=N; ii++)
                {
                    acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, ii, field_name);
                    MEX_DIM_CHECK_VEC_STAGE(fun_name, field, ii, matlab_size, acados_size)
                    if (matlab_size != 0)
                    {
                        ocp_nlp_cost_model_set(config, dims, in, ii, field_name, value);
                    }
                }
            }
            else
            {
                acados_size = ocp_nlp_dims_get_total_from_attr(config, dims, out, field_name);
                MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
                offset = 0;
                for (ii=0; ii<=N; ii++) // TODO implement set_all
                {
                    ocp_nlp_cost_model_set(config, dims, in, ii, field_name, value+offset);
                    tmp_int = ocp_nlp_dims_get_from_attr(config, dims, out, ii, field_name);
                    offset += tmp_int;
                }
            }
        }
        else if (nrhs == min_nrhs + 1) // single stage
        {
            acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, s0, field_name);
            MEX_DIM_CHECK_VEC_STAGE(fun_name, field, s0, matlab_size, acados_size)
            if (matlab_size != 0)
                ocp_nlp_cost_model_set(config, dims, in, s0, field_name, value);
        }
    }
    // initializations
    else if (!strcmp(field, "init_x") || !strcmp(field, "x"))
    {
        if (nrhs == min_nrhs)
        {
            acados_size = ocp_nlp_dims_get_total_from_attr(config, dims, out, "x");
            MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
            ocp_nlp_set_all(solver, in, out, "x", value);
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
            acados_size = ocp_nlp_dims_get_total_from_attr(config, dims, out, "u");
            MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
            ocp_nlp_set_all(solver, in, out, "u", value);
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
        sim_solver_plan_t sim_plan = plan->sim_solver_plan[0];
        sim_solver_t type = sim_plan.sim_solver;
        if (type == IRK)
        {
            if (nrhs == min_nrhs)
            {
                acados_size = ocp_nlp_dims_get_total_from_attr(config, dims, out, "z");
                MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
                ocp_nlp_set_all(solver, in, out, "z", value);
            }
            else // (nrhs == min_nrhs+1)
            {
                acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, s0, "z");
                MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
                ocp_nlp_set(solver, s0, "z_guess", value);
            }
        }
        else
        {
            MEX_FIELD_ONLY_SUPPORTED_FOR_SOLVER(fun_name, "init_z", "irk")
        }
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
                for (ii=0; ii<N; ii++)
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
                for (ii=0; ii<N; ii++)
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
            acados_size = ocp_nlp_dims_get_total_from_attr(config, dims, out, "pi");
            MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
            ocp_nlp_set_all(solver, in, out, "pi", value);
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
            acados_size = ocp_nlp_dims_get_total_from_attr(config, dims, out, "lam");
            MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
            ocp_nlp_set_all(solver, in, out, "lam", value);
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
            acados_size = ocp_nlp_dims_get_total_from_attr(config, dims, out, "sl");
            MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
            ocp_nlp_set_all(solver, in, out, "sl", value);
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
            acados_size = ocp_nlp_dims_get_total_from_attr(config, dims, out, "su");
            MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
            ocp_nlp_set_all(solver, in, out, "su", value);
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
            acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, 0, "p");

            if (acados_size == matlab_size) // setting the same value for all stages
            {
                for (ii=0; ii<=N; ii++)
                {
                    acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "p");
                    MEX_DIM_CHECK_VEC_STAGE(fun_name, field, ii, matlab_size, acados_size);
                    {{ name }}_acados_update_params(capsule, ii, value, matlab_size);
                }
            }
            else
            {
                acados_size = ocp_nlp_dims_get_total_from_attr(config, dims, out, "p");
                MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size);
                ocp_nlp_set_all(solver, in, out, "p", value);
            }
        }
        else if (nrhs == min_nrhs+1) // one stage
        {
            acados_size = ocp_nlp_dims_get_from_attr(config, dims, out, s0, "p");
            MEX_DIM_CHECK_VEC_STAGE(fun_name, field, s0, matlab_size, acados_size)
            {{ name }}_acados_update_params(capsule, s0, value, matlab_size);
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
            for (ii=0; ii<=N; ii++)
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

