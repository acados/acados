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
#include "blasfeo/include/blasfeo_d_aux.h"
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
    // mexPrintf("\nin ocp_get: field%s\n", field);

    int N = dims->N;
    int nu = dims->nu[0];
    int nx = dims->nx[0];
    int nz = dims->nz[0];

    int stage;

    if (nrhs==3)
    {
        stage = mxGetScalar( prhs[2] );
        if (stage < 0 || stage > N)
        {
            sprintf(buffer, "\nocp_get: invalid stage index, got %d\n", stage);
            mexErrMsgTxt(buffer);
        }
        else if (stage == N && strcmp(field, "x") && strcmp(field, "lam") && strcmp(field, "t") && strcmp(field, "sens_x") && strcmp(field, "sl") && strcmp(field, "su") && strcmp(field, "qp_Q") && strcmp(field, "qp_R") && strcmp(field, "qp_S") && strcmp(field, "qp_q") )
        {
            sprintf(buffer, "\nocp_get: invalid stage index, got stage = %d = N, field = %s, only x, lam, t, slacks available at this stage\n", stage, field);
            mexErrMsgTxt(buffer);
        }
    }

    if (!strcmp(field, "x"))
    {
        if (nrhs==2)
        {
            plhs[0] = mxCreateNumericMatrix(nx, N+1, mxDOUBLE_CLASS, mxREAL);
            double *x = mxGetPr( plhs[0] );
            for (ii=0; ii<=N; ii++)
            {
                ocp_nlp_out_get(config, dims, out, ii, "x", x+ii*nx);
            }
        }
        else if (nrhs==3)
        {
            plhs[0] = mxCreateNumericMatrix(nx, 1, mxDOUBLE_CLASS, mxREAL);
            double *x = mxGetPr( plhs[0] );
            ocp_nlp_out_get(config, dims, out, stage, "x", x);
        }
        else
        {
            sprintf(buffer, "\nocp_get: wrong nrhs: %d\n", nrhs);
            mexErrMsgTxt(buffer);
        }
    }
    else if (!strcmp(field, "u"))
    {
        if (nrhs==2)
        {
            plhs[0] = mxCreateNumericMatrix(nu, N, mxDOUBLE_CLASS, mxREAL);
            double *u = mxGetPr( plhs[0] );
            for (ii=0; ii<N; ii++)
            {
                ocp_nlp_out_get(config, dims, out, ii, "u", u+ii*nu);
            }
        }
        else if (nrhs==3)
        {
            plhs[0] = mxCreateNumericMatrix(nu, 1, mxDOUBLE_CLASS, mxREAL);
            double *u = mxGetPr( plhs[0] );
            ocp_nlp_out_get(config, dims, out, stage, "u", u);
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
        else
        {
            sprintf(buffer, "\nocp_get: wrong nrhs: %d\n", nrhs);
            mexErrMsgTxt(buffer);
        }
    }
    else if (!strcmp(field, "t"))
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
        else
        {
            sprintf(buffer, "\nocp_get: wrong nrhs: %d\n", nrhs);
            mexErrMsgTxt(buffer);
        }
    }
    else if (!strcmp(field, "z"))
    {
        if (nrhs==2)
        {
            plhs[0] = mxCreateNumericMatrix(nz, N, mxDOUBLE_CLASS, mxREAL);
            double *z = mxGetPr( plhs[0] );
            for (ii=0; ii<N; ii++)
            {
                ocp_nlp_out_get(config, dims, out, ii, "z", z+ii*nz);
            }
        }
        else if (nrhs==3)
        {
            plhs[0] = mxCreateNumericMatrix(nz, 1, mxDOUBLE_CLASS, mxREAL);
            double *z = mxGetPr( plhs[0] );
            ocp_nlp_out_get(config, dims, out, stage, "z", z);
        }
        else
        {
            sprintf(buffer, "\nocp_get: wrong nrhs: %d\n", nrhs);
            mexErrMsgTxt(buffer);
        }
    }
    else if (!strcmp(field, "pi"))
    {
        if (nrhs==2)
        {
            plhs[0] = mxCreateNumericMatrix(nx, N, mxDOUBLE_CLASS, mxREAL);
            double *pi = mxGetPr( plhs[0] );
            for (ii=0; ii<N; ii++)
            {
                ocp_nlp_out_get(config, dims, out, ii, "pi", pi+ii*nx);
            }
        }
        else if (nrhs==3)
        {
            plhs[0] = mxCreateNumericMatrix(nx, 1, mxDOUBLE_CLASS, mxREAL);
            double *pi = mxGetPr( plhs[0] );
            ocp_nlp_out_get(config, dims, out, stage, "pi", pi);
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
        else
        {
            sprintf(buffer, "\nocp_get: wrong nrhs: %d\n", nrhs);
            mexErrMsgTxt(buffer);
        }
    }
    else if (!strcmp(field, "sens_x"))
    {
        if (nrhs==2)
        {
            plhs[0] = mxCreateNumericMatrix(nx, N+1, mxDOUBLE_CLASS, mxREAL);
            double *x = mxGetPr( plhs[0] );
            for (ii=0; ii<=N; ii++)
            {
                ocp_nlp_out_get(config, dims, sens_out, ii, "x", x+ii*nx);
            }
        }
        else if (nrhs==3)
        {
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
        if (nrhs==2)
        {
            plhs[0] = mxCreateNumericMatrix(nu, N, mxDOUBLE_CLASS, mxREAL);
            double *u = mxGetPr( plhs[0] );
            for (ii=0; ii<N; ii++)
            {
                ocp_nlp_out_get(config, dims, sens_out, ii, "u", u+ii*nu);
            }
        }
        else if (nrhs==3)
        {
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
        if (nrhs==2)
        {
            plhs[0] = mxCreateNumericMatrix(nx, N, mxDOUBLE_CLASS, mxREAL);
            double *pi = mxGetPr( plhs[0] );
            for (ii=0; ii<N; ii++)
            {
                ocp_nlp_out_get(config, dims, sens_out, ii, "pi", pi+ii*nx);
            }
        }
        else if (nrhs==3)
        {
            plhs[0] = mxCreateNumericMatrix(nx, 1, mxDOUBLE_CLASS, mxREAL);
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
        ocp_nlp_get(config, solver, "status", &status);
        *mat_ptr = (double) status;
    }
    else if (!strcmp(field, "sqp_iter"))
    {
        plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
        double *mat_ptr = mxGetPr( plhs[0] );
        int sqp_iter;
        ocp_nlp_get(config, solver, "sqp_iter", &sqp_iter);
        *mat_ptr = (double) sqp_iter;
    }
    else if (!strcmp(field, "time_tot") || !strcmp(field, "time_lin") || !strcmp(field, "time_glob") || !strcmp(field, "time_reg") || !strcmp(field, "time_qp_sol") || !strcmp(field, "time_qp_solver_call") || !strcmp(field, "time_qp_solver") || !strcmp(field, "time_qp_xcond") || !strcmp(field, "time_sim") || !strcmp(field, "time_sim_la") || !strcmp(field, "time_sim_ad"))
    {
        plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
        double *mat_ptr = mxGetPr( plhs[0] );
        ocp_nlp_get(config, solver, field, mat_ptr);
    }
    else if (!strcmp(field, "qp_iter"))
    {
        plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
        double *mat_ptr = mxGetPr( plhs[0] );
        int qp_iter;
        ocp_nlp_get(config, solver, "qp_iter", &qp_iter);
        *mat_ptr = (double) qp_iter;
    }
    else if (!strcmp(field, "stat"))
    {
        int sqp_iter;
        int stat_m, stat_n;
        double *stat;
        ocp_nlp_get(config, solver, "sqp_iter", &sqp_iter);
        ocp_nlp_get(config, solver, "stat_m", &stat_m);
        ocp_nlp_get(config, solver, "stat_n", &stat_n);
        ocp_nlp_get(config, solver, "stat", &stat);
        int min_size = stat_m<sqp_iter+1 ? stat_m : sqp_iter+1;
        plhs[0] = mxCreateNumericMatrix(min_size, stat_n+1, mxDOUBLE_CLASS, mxREAL);
        double *mat_ptr = mxGetPr( plhs[0] );
        for (ii=0; ii<min_size; ii++)
        {
            mat_ptr[ii+0] = ii;
            for (jj=0; jj<stat_n; jj++)
                mat_ptr[ii+(jj+1)*min_size] = stat[jj+ii*stat_n];
        }
    }
    else if (!strcmp(field, "residuals"))
    {
        if (plan->nlp_solver == SQP_RTI)
            ocp_nlp_eval_residuals(solver, in, out);
        plhs[0] = mxCreateNumericMatrix(4, 1, mxDOUBLE_CLASS, mxREAL);
        double *mat_ptr = mxGetPr( plhs[0] );
        ocp_nlp_get(config, solver, "res_stat", &mat_ptr[0]);
        ocp_nlp_get(config, solver, "res_eq", &mat_ptr[1]);
        ocp_nlp_get(config, solver, "res_ineq", &mat_ptr[2]);
        ocp_nlp_get(config, solver, "res_comp", &mat_ptr[3]);
    }
    else if (!strcmp(field, "qp_solver_cond_H"))
    {
        void *qp_in_;
        ocp_nlp_get(config, solver, "qp_xcond_in", &qp_in_);
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
        if(solver_type==1)
        {
            ocp_qp_in *qp_in = qp_in_;
            int *nu = qp_in->dim->nu;
            int *nx = qp_in->dim->nx;

            mxArray *cell_array = mxCreateCellMatrix(N+1, 1);
            plhs[0] = cell_array;

            mxArray *tmp_mat;

            for (ii=0; ii<=N; ii++)
            {
                tmp_mat = mxCreateNumericMatrix(nu[ii]+nx[ii], nu[ii]+nx[ii], mxDOUBLE_CLASS, mxREAL);
                double *mat_ptr = mxGetPr( tmp_mat );
                blasfeo_unpack_dmat(nu[ii]+nx[ii], nu[ii]+nx[ii], qp_in->RSQrq+ii, 0, 0, mat_ptr, nu[ii]+nx[ii]);

                mxSetCell(cell_array, ii, tmp_mat);
            }
        }
        // dense solver
        else if(solver_type==2)
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
             !strcmp(field, "qp_ubu"))
    {
        int out_dims[2];
        if (nrhs==2)
        {
            mxArray *cell_array = mxCreateCellMatrix(N, 1);
            plhs[0] = cell_array;
            mxArray *tmp_mat;

            for (ii=0; ii<N; ii++)
            {
                ocp_nlp_qp_dims_get_from_attr(config, dims, out, ii, &field[3], out_dims);
                tmp_mat = mxCreateNumericMatrix(out_dims[0], out_dims[1], mxDOUBLE_CLASS, mxREAL);
                double *mat_ptr = mxGetPr( tmp_mat );
                ocp_nlp_get_at_stage(config, dims, solver, ii, &field[3], mat_ptr);
                mxSetCell(cell_array, ii, tmp_mat);
            }
        }
        else if (nrhs==3)
        {
            ocp_nlp_qp_dims_get_from_attr(config, dims, out, stage, &field[3], out_dims);
            plhs[0] = mxCreateNumericMatrix(out_dims[0], out_dims[1], mxDOUBLE_CLASS, mxREAL);
            double *mat_ptr = mxGetPr( plhs[0] );
            ocp_nlp_get_at_stage(config, dims, solver, stage, &field[3], mat_ptr);
        }
    }
    else
    {
        MEX_FIELD_NOT_SUPPORTED_SUGGEST(fun_name, field,
             "x, u, z, pi, lam, sl, su, t, sens_x, sens_u, sens_pi, status, sqp_iter, time_tot, time_lin, time_reg, time_qp_sol, stat, qp_solver_cond_H, qp_A, qp_B, qp_Q, qp_R, qp_S, qp_b, qp_q, qp_r");
    }

    return;

}

