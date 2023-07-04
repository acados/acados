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


#include <stdio.h>

#include <stdlib.h>
#include <stddef.h>

// acados
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/dense_qp/dense_qp_hpipm.h"
#include "acados/utils/print.h"
#include "acados/utils/mem.h"
#include "acados_c/dense_qp_interface.h"
// hpipm
#include "hpipm/include/hpipm_d_dense_qp_kkt.h"

#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"

// static dense_qp_res *dense_qp_res_create(dense_qp_dims *dims)
// {
//     int size = dense_qp_res_calculate_size(dims);
//     void *ptr = acados_malloc(size, 1);
//     dense_qp_res *qp_res = dense_qp_res_assign(dims, ptr);
//     return qp_res;
// }

// static dense_qp_res_ws *dense_qp_res_workspace_create(dense_qp_dims *dims)
// {
//     int size = dense_qp_res_workspace_calculate_size(dims);
//     void *ptr = acados_malloc(size, 1);
//     dense_qp_res_ws *res_ws = dense_qp_res_workspace_assign(dims, ptr);
//     return res_ws;
// }

int main() {

    double H[] = {3., 1., 1., 4.};
    double g[] = {2., 1.};

    int idxb[] = {0, 1};
    double d_lb[] = {-1., -1.};
    double d_ub[] = {1., 1.};

    double C[] = {1.2, 1.4, 1.3, 1.5};
    double d_lg[] = {-2., -2.};
    double d_ug[] = {2., 2.};

    int idxs[] = {1, 2, 3};
    double d_ls[] = {0.0, 0.0, 0.0};
    double d_us[] = {0.0, 0.0, 0.0};

    double Zl[] = {1.0, 1.0, 1.0};
    double zl[] = {1.0, 1.0, 1.0};
    double Zu[] = {1.0, 1.0, 1.0};
    double zu[] = {1.0, 1.0, 1.0};

    dense_qp_solver_plan plan;
//    plan.qp_solver = DENSE_QP_QPOASES;
//    plan.qp_solver = DENSE_QP_DAQP;
    plan.qp_solver = DENSE_QP_HPIPM;

    qp_solver_config *config = dense_qp_config_create(&plan);

    dense_qp_dims dims;
    dims.nv = 2;
    dims.ne = 0;  // TODO(dimitris): is this even supported in acados?
    dims.nb = 2;
    dims.ng = 2;
    dims.ns = 3;
    dims.nsb = 1;
    dims.nsg = 2;

    dense_qp_in *qp_in = dense_qp_in_create(config, &dims);

    d_dense_qp_set_all(H, g, NULL, NULL, idxb, d_lb, d_ub, C, d_lg, d_ug, Zl, Zu, zl, zu, idxs, d_ls, d_us, qp_in);

    print_dense_qp_in(qp_in);

    void *opts = dense_qp_opts_create(config, &dims);

    // overwrite default opts
	// dense_qp_hpipm_opts *hpipm_opts = opts;
	// hpipm_opts->hpipm_opts->mu0 = 1e1;

    dense_qp_out *qp_out = dense_qp_out_create(config, &dims);

    dense_qp_solver *qp_solver = dense_qp_create(config, &dims, opts);

    int acados_return = dense_qp_solve(qp_solver, qp_in, qp_out);

    printf("STATUS: %d\n", acados_return);

    // if (acados_return != ACADOS_SUCCESS)
    //     return -1;

    /************************************************
     * compute inf norm of residuals
     ************************************************/

    double res[4];

#if 0
	blasfeo_print_exp_tran_dvec(dims.nv+2*dims.ns, qp_res->res_g, 0);
	blasfeo_print_exp_tran_dvec(dims.ne, qp_res->res_b, 0);
	blasfeo_print_exp_tran_dvec(2*dims.nb+2*dims.ng+2*dims.ns, qp_res->res_d, 0);
	blasfeo_print_exp_tran_dvec(2*dims.nb+2*dims.ng+2*dims.ns, qp_res->res_m, 0);
#endif

    printf("v=\n");
    blasfeo_print_exp_tran_dvec(dims.nv+2*dims.ns, qp_out->v, 0);

    dense_qp_inf_norm_residuals(&dims, qp_in, qp_out, res);
    printf("\ninf norm res: %e, %e, %e, %e\n\n", res[0], res[1], res[2], res[3]);

    free(config);
    free(opts);
    free(qp_in);
    free(qp_out);
    free(qp_solver);
}
