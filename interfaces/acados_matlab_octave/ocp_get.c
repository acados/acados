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

// system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
// acados
#include "acados_c/ocp_nlp_interface.h"
#include "acados/dense_qp/dense_qp_common.h"
#include "blasfeo_d_aux.h"
// mex
#include "mex.h"
#include "mex_macros.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    long long *ptr;
    char fun_name[50] = "ocp_get";
    char buffer [300]; // for error messages

    int ii, jj, length;

    /* RHS */

    // C_ocp
    const mxArray *C_ocp = prhs[0];
    // config
    ptr = (long long *) mxGetData( mxGetField( C_ocp, 0, "config" ) );
    ocp_nlp_config *config = (ocp_nlp_config *) ptr[0];
    // dims
    ptr = (long long *) mxGetData( mxGetField( C_ocp, 0, "dims" ) );
    ocp_nlp_dims *dims = (ocp_nlp_dims *) ptr[0];
    // out
    ptr = (long long *) mxGetData( mxGetField( C_ocp, 0, "out" ) );
    ocp_nlp_out *out = (ocp_nlp_out *) ptr[0];
    // solver
    ptr = (long long *) mxGetData( mxGetField( C_ocp, 0, "solver" ) );
    ocp_nlp_solver *solver = (ocp_nlp_solver *) ptr[0];
    // in
    ptr = (long long *) mxGetData( mxGetField( prhs[0], 0, "in" ) );
    ocp_nlp_in *in = (ocp_nlp_in *) ptr[0];
    // sens_out
    ptr = (long long *) mxGetData( mxGetField( C_ocp, 0, "sens_out" ) );
    ocp_nlp_out *sens_out = (ocp_nlp_out *) ptr[0];
    // plan
    ptr = (long long *) mxGetData( mxGetField( C_ocp, 0, "plan" ) );
    ocp_nlp_plan_t *plan = (ocp_nlp_plan_t *) ptr[0];

    // field
    char *field = mxArrayToString( prhs[1] );
    // mexPrintf("\nin ocp_get: field %s\n", field);
    if (field == NULL)
    {
        mexErrMsgTxt("got NULL pointer for field, maybe you put in a \" instead of \' in MATLAB.\n");
    }

    int N = dims->N;
    int stage;
    int iteration;

    if (nrhs>=3)
    {
        stage = mxGetScalar( prhs[2] );
        if (stage < 0 || stage > N)
        {
            sprintf(buffer, "\nocp_get: invalid stage index, got %d\n", stage);
            mexErrMsgTxt(buffer);
        }
        else if (stage == N && strcmp(field, "x") &&
                               strcmp(field, "lam") &&
                               strcmp(field, "p") &&
                               strcmp(field, "sens_x") &&
                               strcmp(field, "sl") &&
                               strcmp(field, "su") &&
                               strcmp(field, "qp_Q") &&
                               strcmp(field, "qp_q") &&
                               strcmp(field, "qp_C") &&
                               strcmp(field, "qp_lg") &&
                               strcmp(field, "qp_ug") &&
                               strcmp(field, "qp_lbx") &&
                               strcmp(field, "qp_ubx") &&
                               strcmp(field, "qp_zl") &&
                               strcmp(field, "qp_zu") &&
                               strcmp(field, "qp_Zl") &&
                               strcmp(field, "qpscaling_obj") &&
                               strcmp(field, "qpscaling_constr") &&
                               strcmp(field, "qp_Zu"))
        {
            sprintf(buffer, "\nocp_get: invalid stage index, got stage = %d = N, field = %s, field not available at final shooting node\n", stage, field);
            mexErrMsgTxt(buffer);
        }
    }

    if (nrhs == 4)
    {
        iteration = mxGetScalar(prhs[3]);
        int nlp_iter;
        ocp_nlp_get(solver, "nlp_iter", &nlp_iter);
        if (iteration < 0 || iteration > nlp_iter)
        {
            sprintf(buffer, "\nocp_get: invalid iteration index, got stage = %d, should be nonnegative and <= nlp_iter = %d\n", iteration, nlp_iter);
            mexErrMsgTxt(buffer);
        }
    }

    if (!strcmp(field, "x"))
    {
        int nx;
        if (nrhs==2)
        {
            int nx0 = ocp_nlp_dims_get_from_attr(config, dims, out, 0, "x");
            plhs[0] = mxCreateNumericMatrix(nx0, N+1, mxDOUBLE_CLASS, mxREAL);
            double *x = mxGetPr( plhs[0] );
            for (ii=0; ii<=N; ii++)
            {
                nx = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "x");
                if (nx != nx0)
                {
                    sprintf(buffer, "\nocp_get: cannot get multiple %s values at once due to varying dimension", field);
                    mexErrMsgTxt(buffer);
                }
            }
            ocp_nlp_get_all(solver, in, out, "x", x);

        }
        else if (nrhs==3)
        {
            nx = ocp_nlp_dims_get_from_attr(config, dims, out, stage, "x");
            plhs[0] = mxCreateNumericMatrix(nx, 1, mxDOUBLE_CLASS, mxREAL);
            double *x = mxGetPr( plhs[0] );
            ocp_nlp_out_get(config, dims, out, stage, "x", x);
        }
        else if (nrhs==4)
        {
            nx = ocp_nlp_dims_get_from_attr(config, dims, out, stage, "x");
            plhs[0] = mxCreateNumericMatrix(nx, 1, mxDOUBLE_CLASS, mxREAL);
            double *x = mxGetPr(plhs[0]);
            ocp_nlp_get_from_iterate(solver, iteration, stage, "x", x);
        }
        else
        {
            sprintf(buffer, "\nocp_get: wrong nrhs: %d\n", nrhs);
            mexErrMsgTxt(buffer);
        }
    }
    else if (!strcmp(field, "u"))
    {
        int nu;
        if (nrhs==2)
        {
            int nu0 = ocp_nlp_dims_get_from_attr(config, dims, out, 0, "u");
            plhs[0] = mxCreateNumericMatrix(nu0, N, mxDOUBLE_CLASS, mxREAL);
            double *u = mxGetPr( plhs[0] );
            for (ii=0; ii<N; ii++)
            {
                nu = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "u");
                if (nu != nu0)
                {
                    sprintf(buffer, "\nocp_get: cannot get multiple %s values at once due to varying dimension", field);
                    mexErrMsgTxt(buffer);
                }
            }
            ocp_nlp_get_all(solver, in, out, "u", u);
        }
        else if (nrhs==3)
        {
            nu = ocp_nlp_dims_get_from_attr(config, dims, out, stage, "u");
            plhs[0] = mxCreateNumericMatrix(nu, 1, mxDOUBLE_CLASS, mxREAL);
            double *u = mxGetPr( plhs[0] );
            ocp_nlp_out_get(config, dims, out, stage, "u", u);
        }
        else if (nrhs==4)
        {
            nu = ocp_nlp_dims_get_from_attr(config, dims, out, stage, "u");
            plhs[0] = mxCreateNumericMatrix(nu, 1, mxDOUBLE_CLASS, mxREAL);
            double *u = mxGetPr(plhs[0]);
            ocp_nlp_get_from_iterate(solver, iteration, stage, "u", u);
        }
        else
        {
            sprintf(buffer, "\nocp_get: wrong nrhs: %d\n", nrhs);
            mexErrMsgTxt(buffer);
        }
    }
    else if (!strcmp(field, "sl") || !strcmp(field, "su"))
    {
        if (nrhs==2)
        {
            sprintf(buffer, "\nocp_get: values %s can only be accessed stagewise\n", field);
            mexErrMsgTxt(buffer);
        }
        else if (nrhs==3)
        {
            length = ocp_nlp_dims_get_from_attr(config, dims, out, stage, field);
            plhs[0] = mxCreateNumericMatrix(length, 1, mxDOUBLE_CLASS, mxREAL);
            double *value = mxGetPr( plhs[0] );
            ocp_nlp_out_get(config, dims, out, stage, field, value);
        }
        else if (nrhs==4)
        {
            length = ocp_nlp_dims_get_from_attr(config, dims, out, stage, field);
            plhs[0] = mxCreateNumericMatrix(length, 1, mxDOUBLE_CLASS, mxREAL);
            double *value = mxGetPr(plhs[0]);
            ocp_nlp_get_from_iterate(solver, iteration, stage, field, value);
        }
        else
        {
            sprintf(buffer, "\nocp_get: wrong nrhs: %d\n", nrhs);
            mexErrMsgTxt(buffer);
        }
    }
    else if (!strcmp(field, "z"))
    {
        int nz;
        if (nrhs==2)
        {
            int nz0 = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "z");
            plhs[0] = mxCreateNumericMatrix(nz0, N, mxDOUBLE_CLASS, mxREAL);
            double *z = mxGetPr( plhs[0] );
            for (ii=0; ii<N; ii++)
            {
                nz = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "z");
                if (nz != nz0)
                {
                    sprintf(buffer, "\nocp_get: cannot get multiple %s values at once due to varying dimension", field);
                    mexErrMsgTxt(buffer);
                }
            }
            ocp_nlp_get_all(solver, in, out, "z", z);
        }
        else if (nrhs==3)
        {
            nz = ocp_nlp_dims_get_from_attr(config, dims, out, stage, "z");
            plhs[0] = mxCreateNumericMatrix(nz, 1, mxDOUBLE_CLASS, mxREAL);
            double *z = mxGetPr( plhs[0] );
            ocp_nlp_out_get(config, dims, out, stage, "z", z);
        }
        else if (nrhs==4)
        {
            nz = ocp_nlp_dims_get_from_attr(config, dims, out, stage, "z");
            plhs[0] = mxCreateNumericMatrix(nz, 1, mxDOUBLE_CLASS, mxREAL);
            double *z = mxGetPr(plhs[0]);
            ocp_nlp_get_from_iterate(solver, iteration, stage, "z", z);
        }
        else
        {
            sprintf(buffer, "\nocp_get: wrong nrhs: %d\n", nrhs);
            mexErrMsgTxt(buffer);
        }
    }
    else if (!strcmp(field, "pi"))
    {
        int npi;

        if (nrhs==2)
        {
            int npi0 = ocp_nlp_dims_get_from_attr(config, dims, out, 0, "pi");
            plhs[0] = mxCreateNumericMatrix(npi0, N, mxDOUBLE_CLASS, mxREAL);
            double *pi = mxGetPr( plhs[0] );
            for (ii=0; ii<N; ii++)
            {
                npi = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "pi");
                if (npi != npi0)
                {
                    sprintf(buffer, "\nocp_get: cannot get multiple %s values at once due to varying dimension", field);
                    mexErrMsgTxt(buffer);
                }
            }
            ocp_nlp_get_all(solver, in, out, "pi", pi);
        }
        else if (nrhs==3)
        {
            npi = ocp_nlp_dims_get_from_attr(config, dims, out, stage, "pi");
            plhs[0] = mxCreateNumericMatrix(npi, 1, mxDOUBLE_CLASS, mxREAL);
            double *pi = mxGetPr( plhs[0] );
            ocp_nlp_out_get(config, dims, out, stage, "pi", pi);
        }
        else if (nrhs==4)
        {
            npi = ocp_nlp_dims_get_from_attr(config, dims, out, stage, "pi");
            plhs[0] = mxCreateNumericMatrix(npi, 1, mxDOUBLE_CLASS, mxREAL);
            double *pi = mxGetPr(plhs[0]);
            ocp_nlp_get_from_iterate(solver, iteration, stage, "pi", pi);
        }
        else
        {
            sprintf(buffer, "\nocp_get: wrong nrhs: %d\n", nrhs);
            mexErrMsgTxt(buffer);
        }
    }
    else if (!strcmp(field, "p"))
    {
        if (nrhs==2)
        {
            sprintf(buffer, "\nocp_get: field p: only supported for a single shooting node.\n");
            mexErrMsgTxt(buffer);
        }
        else if (nrhs==3)
        {
            int np = ocp_nlp_dims_get_from_attr(config, dims, out, stage, field);
            plhs[0] = mxCreateNumericMatrix(np, 1, mxDOUBLE_CLASS, mxREAL);
            double *p = mxGetPr( plhs[0] );
            ocp_nlp_in_get(config, dims, in, stage, "p", p);
        }
        else
        {
            sprintf(buffer, "\nocp_get: wrong nrhs: %d\n", nrhs);
            mexErrMsgTxt(buffer);
        }
    }
    else if (!strcmp(field, "lam"))
    {
        if (nrhs==2)
        {
            sprintf(buffer, "\nocp_get: field lam: only supported for a single shooting node.\n");
            mexErrMsgTxt(buffer);
        }
        else if (nrhs==3)
        {
            int nlam = ocp_nlp_dims_get_from_attr(config, dims, out, stage, "lam");
            plhs[0] = mxCreateNumericMatrix(nlam, 1, mxDOUBLE_CLASS, mxREAL);
            double *lam = mxGetPr( plhs[0] );
            ocp_nlp_out_get(config, dims, out, stage, "lam", lam);
        }
        else if (nrhs==4)
        {
            int nlam = ocp_nlp_dims_get_from_attr(config, dims, out, stage, "lam");
            plhs[0] = mxCreateNumericMatrix(nlam, 1, mxDOUBLE_CLASS, mxREAL);
            double *lam = mxGetPr(plhs[0]);
            ocp_nlp_get_from_iterate(solver, iteration, stage, "lam", lam);
        }
        else
        {
            sprintf(buffer, "\nocp_get: wrong nrhs: %d\n", nrhs);
            mexErrMsgTxt(buffer);
        }
    }
    else if (!strcmp(field, "sens_x"))
    {
        int nx;
        if (nrhs==2)
        {
            int nx0 = ocp_nlp_dims_get_from_attr(config, dims, out, 0, "x");
            for (ii=0; ii<=N; ii++)
            {
                nx = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "x");
                if (nx != nx0)
                {
                    sprintf(buffer, "\nocp_get: cannot get multiple %s values at once due to varying dimension", field);
                    mexErrMsgTxt(buffer);
                }
            }

            plhs[0] = mxCreateNumericMatrix(nx, N+1, mxDOUBLE_CLASS, mxREAL);
            double *x = mxGetPr( plhs[0] );
            for (ii=0; ii<=N; ii++)
            {
                ocp_nlp_out_get(config, dims, sens_out, ii, "x", x+ii*nx);
            }
        }
        else if (nrhs==3)
        {
            nx = ocp_nlp_dims_get_from_attr(config, dims, out, stage, "x");
            plhs[0] = mxCreateNumericMatrix(nx, 1, mxDOUBLE_CLASS, mxREAL);
            double *x = mxGetPr( plhs[0] );
            ocp_nlp_out_get(config, dims, sens_out, stage, "x", x);
        }
        else
        {
            sprintf(buffer, "\nocp_get: wrong nrhs: %d\n", nrhs);
            mexErrMsgTxt(buffer);
        }
    }
    else if (!strcmp(field, "sens_u"))
    {
        // TODO: fix for MOCP
        int nu;
        if (nrhs==2)
        {
            int nu0 = ocp_nlp_dims_get_from_attr(config, dims, out, 0, "u");
            for (ii=0; ii<N; ii++)
            {
                nu = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "u");
                if (nu0 != nu)
                {
                    sprintf(buffer, "\nocp_get: cannot get multiple %s values at once due to varying dimension", field);
                    mexErrMsgTxt(buffer);
                }
            }
            plhs[0] = mxCreateNumericMatrix(nu, N, mxDOUBLE_CLASS, mxREAL);
            double *u = mxGetPr( plhs[0] );
            for (ii=0; ii<N; ii++)
            {
                ocp_nlp_out_get(config, dims, sens_out, ii, "u", u+ii*nu);
            }
        }
        else if (nrhs==3)
        {
            nu = ocp_nlp_dims_get_from_attr(config, dims, out, stage, "u");
            plhs[0] = mxCreateNumericMatrix(nu, 1, mxDOUBLE_CLASS, mxREAL);
            double *u = mxGetPr( plhs[0] );
            ocp_nlp_out_get(config, dims, sens_out, stage, "u", u);
        }
        else
        {
            sprintf(buffer, "\nocp_get: wrong nrhs: %d\n", nrhs);
            mexErrMsgTxt(buffer);
        }
    }
    else if (!strcmp(field, "sens_pi"))
    {
        int npi;
        if (nrhs==2)
        {
            int npi0 = ocp_nlp_dims_get_from_attr(config, dims, out, 0, "pi");
            plhs[0] = mxCreateNumericMatrix(npi0, N, mxDOUBLE_CLASS, mxREAL);

            double *pi = mxGetPr( plhs[0] );
            for (ii=0; ii<N; ii++)
            {
                if (npi0 != npi)
                {
                    sprintf(buffer, "\nocp_get: cannot get multiple %s values at once due to varying dimension", field);
                    mexErrMsgTxt(buffer);
                }
                npi = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "pi");
                ocp_nlp_out_get(config, dims, sens_out, ii, "pi", pi+ii*npi);
            }
        }
        else if (nrhs==3)
        {
            npi = ocp_nlp_dims_get_from_attr(config, dims, out, stage, "pi");
            plhs[0] = mxCreateNumericMatrix(npi, 1, mxDOUBLE_CLASS, mxREAL);
            double *pi = mxGetPr( plhs[0] );
            ocp_nlp_out_get(config, dims, sens_out, stage, "pi", pi);
        }
        else
        {
            mexPrintf("\nocp_get: wrong nrhs: %d\n", nrhs);
            return;
        }
    }
    else if (!strcmp(field, "status"))
    {
        plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
        double *mat_ptr = mxGetPr( plhs[0] );
        int status;
        ocp_nlp_get(solver, "status", &status);
        *mat_ptr = (double) status;
    }
    else if (!strcmp(field, "cost_value"))
    {
        plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
        double *out_data = mxGetPr( plhs[0] );
        ocp_nlp_eval_cost(solver, in, out);
        ocp_nlp_get(solver, "cost_value", out_data);
    }
    else if (!strcmp(field, "constraint_violation"))
    {
        int out_dims[2];
        // evaluate
        ocp_nlp_eval_constraints(solver, in, out);
        // create output
        mxArray *cell_array = mxCreateCellMatrix(N+1, 1);
        plhs[0] = cell_array;
        mxArray *tmp_mat;
        for (ii=0; ii<N+1; ii++)
        {
            ocp_nlp_constraint_dims_get_from_attr(config, dims, out, ii, "ineq_fun", out_dims);
            tmp_mat = mxCreateNumericMatrix(out_dims[0], 1, mxDOUBLE_CLASS, mxREAL);
            double *mat_ptr = mxGetPr( tmp_mat );
            ocp_nlp_get_at_stage(solver, ii, "ineq_fun", mat_ptr);
            mxSetCell(cell_array, ii, tmp_mat);
        }
    }
    else if (!strcmp(field, "res_stat_all"))
    {
        int nv;
        mxArray *cell_array = mxCreateCellMatrix(N+1, 1);
        plhs[0] = cell_array;
        mxArray *tmp_mat;
        for (ii=0; ii<N+1; ii++)
        {
            nv = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "nv");
            tmp_mat = mxCreateNumericMatrix(nv, 1, mxDOUBLE_CLASS, mxREAL);
            double *mat_ptr = mxGetPr( tmp_mat );
            ocp_nlp_get_at_stage(solver, ii, "res_stat", mat_ptr);
            mxSetCell(cell_array, ii, tmp_mat);
        }
    }
    else if (!strcmp(field, "res_eq_all"))
    {
        int npi;
        mxArray *cell_array = mxCreateCellMatrix(N, 1);
        plhs[0] = cell_array;
        mxArray *tmp_mat;
        for (ii=0; ii<N; ii++)
        {
            npi = ocp_nlp_dims_get_from_attr(config, dims, out, ii, "pi");
            tmp_mat = mxCreateNumericMatrix(npi, 1, mxDOUBLE_CLASS, mxREAL);
            double *mat_ptr = mxGetPr( tmp_mat );
            ocp_nlp_get_at_stage(solver, ii, "res_eq", mat_ptr);
            mxSetCell(cell_array, ii, tmp_mat);
        }
    }
    else if (!strcmp(field, "sqp_iter") || !strcmp(field, "nlp_iter") || !strcmp(field, "qpscaling_status"))  // int fields
    {
        plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
        double *mat_ptr = mxGetPr( plhs[0] );
        int value;
        ocp_nlp_get(solver, field, &value);
        *mat_ptr = (double) value;
    }
    else if (!strcmp(field, "time_tot") || !strcmp(field, "time_lin") || !strcmp(field, "time_glob") || !strcmp(field, "time_reg") || !strcmp(field, "time_qp_sol") || !strcmp(field, "time_qp_solver_call") || !strcmp(field, "time_qp_solver") || !strcmp(field, "time_qp_xcond") || !strcmp(field, "time_sim") || !strcmp(field, "time_sim_la") || !strcmp(field, "time_sim_ad") || !strcmp(field, "time_qp_scaling"))
    {
        plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
        double *mat_ptr = mxGetPr( plhs[0] );
        ocp_nlp_get(solver, field, mat_ptr);
    }
    else if (!strcmp(field, "qp_iter"))
    {
        plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
        double *mat_ptr = mxGetPr( plhs[0] );
        int qp_iter;
        ocp_nlp_get(solver, "qp_iter", &qp_iter);
        *mat_ptr = (double) qp_iter;
    }
    else if (!strcmp(field, "stat"))
    {
        int nlp_iter;
        int stat_m, stat_n;
        double *stat;
        ocp_nlp_get(solver, "nlp_iter", &nlp_iter);
        ocp_nlp_get(solver, "stat_m", &stat_m);
        ocp_nlp_get(solver, "stat_n", &stat_n);
        ocp_nlp_get(solver, "stat", &stat);
        int min_size = stat_m<nlp_iter+1 ? stat_m : nlp_iter+1;
        plhs[0] = mxCreateNumericMatrix(min_size, stat_n+1, mxDOUBLE_CLASS, mxREAL);
        double *mat_ptr = mxGetPr( plhs[0] );
        for (ii=0; ii<min_size; ii++)
        {
            mat_ptr[ii+0] = ii;
            for (jj=0; jj<stat_n; jj++)
                mat_ptr[ii+(jj+1)*min_size] = stat[jj+ii*stat_n];
        }
    }
    else if (!strcmp(field, "primal_step_norm") || !strcmp(field, "dual_step_norm"))
    {
        int nlp_iter;
        ocp_nlp_get(solver, "nlp_iter", &nlp_iter);
        plhs[0] = mxCreateNumericMatrix(nlp_iter, 1, mxDOUBLE_CLASS, mxREAL);
        double *mat_ptr = mxGetPr( plhs[0] );
        ocp_nlp_get(solver, field, mat_ptr);
    }
    else if (!strcmp(field, "qpscaling_constr"))
    {
        int size = ocp_nlp_dims_get_from_attr(config, dims, out, stage, "qpscaling_constr");
        plhs[0] = mxCreateNumericMatrix(size, 1, mxDOUBLE_CLASS, mxREAL);
        double *mat_ptr = mxGetPr( plhs[0] );
        ocp_nlp_get_at_stage(solver, stage, field, mat_ptr);
    }
    else if (!strcmp(field, "qpscaling_obj"))
    {
        int size = 1;
        plhs[0] = mxCreateNumericMatrix(size, 1, mxDOUBLE_CLASS, mxREAL);
        double *mat_ptr = mxGetPr( plhs[0] );
        ocp_nlp_get_at_stage(solver, stage, field, mat_ptr);
    }
    else if (!strcmp(field, "residuals"))
    {
        if (plan->nlp_solver == SQP_RTI)
            ocp_nlp_eval_residuals(solver, in, out);
        plhs[0] = mxCreateNumericMatrix(4, 1, mxDOUBLE_CLASS, mxREAL);
        double *mat_ptr = mxGetPr( plhs[0] );
        ocp_nlp_get(solver, "res_stat", &mat_ptr[0]);
        ocp_nlp_get(solver, "res_eq", &mat_ptr[1]);
        ocp_nlp_get(solver, "res_ineq", &mat_ptr[2]);
        ocp_nlp_get(solver, "res_comp", &mat_ptr[3]);
    }
    else if (!strcmp(field, "qp_solver_cond_H"))
    {
        void *qp_in_;
        ocp_nlp_get(solver, "qp_xcond_in", &qp_in_);
        int solver_type = 0;
        if (plan->ocp_qp_solver_plan.qp_solver==PARTIAL_CONDENSING_HPIPM)
            solver_type=1;
        if (plan->ocp_qp_solver_plan.qp_solver==FULL_CONDENSING_HPIPM)
            solver_type=2;
#if defined(ACADOS_WITH_QPOASES)
        if (plan->ocp_qp_solver_plan.qp_solver==FULL_CONDENSING_QPOASES)
            solver_type=2;
#endif
#if defined(ACADOS_WITH_DAQP)
        if (plan->ocp_qp_solver_plan.qp_solver==FULL_CONDENSING_DAQP)
            solver_type=2;
#endif
        // ocp solver (not dense)
        if (solver_type==1)
        {
            ocp_qp_in *qp_in = qp_in_;
            int *nu = qp_in->dim->nu;
            int *nx = qp_in->dim->nx;
            int cond_N = qp_in->dim->N;

            mxArray *cell_array = mxCreateCellMatrix(cond_N+1, 1);
            plhs[0] = cell_array;

            mxArray *tmp_mat;

            for (ii=0; ii<=cond_N; ii++)
            {
                tmp_mat = mxCreateNumericMatrix(nu[ii]+nx[ii], nu[ii]+nx[ii], mxDOUBLE_CLASS, mxREAL);
                double *mat_ptr = mxGetPr( tmp_mat );
                blasfeo_unpack_dmat(nu[ii]+nx[ii], nu[ii]+nx[ii], qp_in->RSQrq+ii, 0, 0, mat_ptr, nu[ii]+nx[ii]);

                mxSetCell(cell_array, ii, tmp_mat);
            }
        }
        // dense solver
        else if (solver_type==2)
        {
            dense_qp_in *qp_in = qp_in_;
            int nv = qp_in->dim->nv;
            plhs[0] = mxCreateNumericMatrix(nv, nv, mxDOUBLE_CLASS, mxREAL);
            double *mat_ptr = mxGetPr( plhs[0] );
            blasfeo_unpack_dmat(nv, nv, qp_in->Hv, 0, 0, mat_ptr, nv);
        }
        else
        {
            mexPrintf("\nerror: ocp_get: qp_solver_cond_H: unsupported solver\n");
            exit(1);
        }
    }
    else if (!strcmp(field, "qp_A") || !strcmp(field, "qp_B") || !strcmp(field, "qp_Q") ||
             !strcmp(field, "qp_R") || !strcmp(field, "qp_S") || !strcmp(field, "qp_b") ||
             !strcmp(field, "qp_q") || !strcmp(field, "qp_r") || !strcmp(field, "qp_C") ||
             !strcmp(field, "qp_D") || !strcmp(field, "qp_lg") || !strcmp(field, "qp_ug") ||
             !strcmp(field, "qp_lbx") || !strcmp(field, "qp_ubx") || !strcmp(field, "qp_lbu") ||
             !strcmp(field, "qp_ubu") || !strcmp(field, "qp_zl") || !strcmp(field, "qp_zu") ||
             !strcmp(field, "qp_Zl") || !strcmp(field, "qp_Zu"))
    {
        int out_dims[2];
        if (nrhs==2)
        {
            // fields that dont exist at last node
            int cell_size;
            if (!strcmp(field, "qp_A") || !strcmp(field, "qp_B") || !strcmp(field, "qp_R") || !strcmp(field, "qp_S") ||
                !strcmp(field, "qp_r") || !strcmp(field, "qp_lbu") || !strcmp(field, "qp_ubu"))
            {
                cell_size = N;
            }
            else
            {
                cell_size = N+1;
            }
            mxArray *cell_array = mxCreateCellMatrix(cell_size, 1);
            plhs[0] = cell_array;
            mxArray *tmp_mat;
            for (ii=0; ii<cell_size; ii++)
            {
                ocp_nlp_qp_dims_get_from_attr(config, dims, out, ii, &field[3], out_dims);
                tmp_mat = mxCreateNumericMatrix(out_dims[0], out_dims[1], mxDOUBLE_CLASS, mxREAL);
                double *mat_ptr = mxGetPr( tmp_mat );
                ocp_nlp_get_at_stage(solver, ii, &field[3], mat_ptr);
                mxSetCell(cell_array, ii, tmp_mat);
            }
        }
        else if (nrhs==3)
        {
            ocp_nlp_qp_dims_get_from_attr(config, dims, out, stage, &field[3], out_dims);
            plhs[0] = mxCreateNumericMatrix(out_dims[0], out_dims[1], mxDOUBLE_CLASS, mxREAL);
            double *mat_ptr = mxGetPr( plhs[0] );
            ocp_nlp_get_at_stage(solver, stage, &field[3], mat_ptr);
        }
    }
    else
    {
        MEX_FIELD_NOT_SUPPORTED_SUGGEST(fun_name, field,
             "x, u, z, pi, lam, sl, su, t, sens_x, sens_u, sens_pi, status, sqp_iter, nlp_iter, time_tot, time_lin, time_reg, time_qp_sol, stat, qp_solver_cond_H, qp_A, qp_B, qp_Q, qp_R, qp_S, qp_b, qp_q, qp_r");
    }

    return;
}
