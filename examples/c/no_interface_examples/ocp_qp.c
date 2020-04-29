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

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include "acados/utils/print.h"
#include "acados_c/ocp_qp_interface.h"

#define QP_HORIZON 5

int main() {

    double A[] = {1, 0, 1, 1};
    double B[] = {0, 1};
    double b[] = {0, 0};

    double Q[] = {1, 0, 0, 1};
    double S[] = {0, 0};
    double R[] = {1};
    double q[] = {1, 1};
    double r[] = {0};

    double x0[] = {1, 1};
    int idxb0[] = {0, 1};


    ocp_qp_solver_plan plan;
    plan.qp_solver = FULL_CONDENSING_HPIPM; // FULL_CONDENSING_QP_OASES

    ocp_qp_xcond_solver_config *config = ocp_qp_xcond_solver_config_create(plan);

    int N = QP_HORIZON;
    ocp_qp_xcond_solver_dims *dims = ocp_qp_xcond_solver_dims_create(config, N);

    int nx = 2;
    int nu = 1;
    int nu_e = 0;
    // here: no general linear constraints (ng), soft constraints (ns, nsbx, nsbu, nsg)

    for (int i = 0; i < N+1; i++)
    {
        ocp_qp_xcond_solver_dims_set(config, dims, i, "nx", &nx);
        ocp_qp_xcond_solver_dims_set(config, dims, i, "nu", &nu);
    }
    // last stage has no controls
    ocp_qp_xcond_solver_dims_set(config, dims, N+1, "nu", &nu_e);

    ocp_qp_xcond_solver_dims_set(config, dims, 0, "nbx", &nx);

    ocp_qp_in *qp_in = ocp_qp_in_create_from_xcond_dims(dims);

    for (int i = 0; i < N; i++)
    {
        ocp_qp_in_set(config, qp_in, i, "A", A);
        ocp_qp_in_set(config, qp_in, i, "B", B);
        ocp_qp_in_set(config, qp_in, i, "b", b);
        ocp_qp_in_set(config, qp_in, i, "Q", Q);
        ocp_qp_in_set(config, qp_in, i, "S", S);
        ocp_qp_in_set(config, qp_in, i, "R", R);
        ocp_qp_in_set(config, qp_in, i, "q", q);
        ocp_qp_in_set(config, qp_in, i, "r", r);
    }
    ocp_qp_in_set(config, qp_in, 0, "idxbx", idxb0);
    ocp_qp_in_set(config, qp_in, 0, "lbx", x0);
    ocp_qp_in_set(config, qp_in, 0, "ubx", x0);

    ocp_qp_in_set(config, qp_in, N, "Q", Q);
    ocp_qp_in_set(config, qp_in, N, "S", S);
    ocp_qp_in_set(config, qp_in, N, "q", q);
    ocp_qp_in_set(config, qp_in, N, "r", r);

    print_ocp_qp_in(qp_in);

    printf("\nqp_in dimensions:\n");
    print_ocp_qp_dims(qp_in->dim);


    void *opts = ocp_qp_xcond_solver_opts_create(config, dims);

    ocp_qp_out *qp_out = ocp_qp_out_create_from_xcond_dims(dims);

    // TODO(dimitris): only have N2 in one place!!
    // printf("N2 in config = %d\n", config->N2);
    // printf("N2 in opts = %d\n", ((ocp_qp_partial_condensing_opts *)(((ocp_qp_partial_condensing_solver_opts *)opts)->pcond_opts))->N2);

    ocp_qp_solver *qp_solver = ocp_qp_create(config, dims, opts);

    int acados_return = ocp_qp_solve(qp_solver, qp_in, qp_out);

    if (acados_return != ACADOS_SUCCESS)
        return -1;

    print_ocp_qp_out(qp_out);

    /************************************************
     * compute inf norm of residuals
     ************************************************/

    double res[4];
    ocp_qp_inf_norm_residuals(dims->orig_dims, qp_in, qp_out, res);
    printf("\ninf norm res: %e, %e, %e, %e\n\n", res[0], res[1], res[2], res[3]);

    // ocp_qp_info *info = (ocp_qp_info *)qp_out->misc;
    // print_ocp_qp_info(info);

    ocp_qp_xcond_solver_dims_free(dims);
    ocp_qp_xcond_solver_config_free(config);
    ocp_qp_xcond_solver_opts_free(opts);
    ocp_qp_in_free(qp_in);
    ocp_qp_out_free(qp_out);
    ocp_qp_solver_destroy(qp_solver);
}
