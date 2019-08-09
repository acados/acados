/*
 * Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren, Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor, Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan, Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

// #include "acados/ocp_qp/ocp_qp_common.h"
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
    plan.qp_solver = FULL_CONDENSING_HPIPM;

    ocp_qp_xcond_solver_config *config = ocp_qp_config_create(plan);

    ocp_qp_dims dims;

    int nx[] = {2, 2, 2, 2, 2, 2};
    int nu[] = {1, 1, 1, 1, 1, 0};
    int nb[] = {2, 0, 0, 0, 0, 0};
    int ng[] = {0, 0, 0, 0, 0, 0};
    int ns[] = {0, 0, 0, 0, 0, 0};
    int nsbx[] = {0, 0, 0, 0, 0, 0};
    int nsbu[] = {0, 0, 0, 0, 0, 0};
    int nsg[] = {0, 0, 0, 0, 0, 0};
    int nbx[] = {2, 0, 0, 0, 0, 0};
    int nbu[] = {0, 0, 0, 0, 0, 0};

    dims.N = QP_HORIZON;
    dims.nx = nx;
    dims.nu = nu;
    dims.nb = nb;
    dims.ng = ng;
    dims.ns = ns;
    dims.nsbx = nsbx;
    dims.nsbu = nsbu;
    dims.nsg = nsg;
    dims.nbx = nbx;
    dims.nbu = nbu;

    ocp_qp_in *qp_in = ocp_qp_in_create(config, &dims);

    double *hA[] = {A, A, A, A, A};
    double *hB[] = {B, B, B, B, B};
    double *hb[] = {b, b, b, b, b};
    double *hQ[] = {Q, Q, Q, Q, Q, Q};
    double *hS[] = {S, S, S, S, S, S};
    double *hR[] = {R, R, R, R, R};
    double *hq[] = {q, q, q, q, q, q};
    double *hr[] = {r, r, r, r, r, r};
    int *hidxbx[] = {idxb0};
    double *hlbx[] = {x0};
    double *hubx[] = {x0};

    d_ocp_qp_set_all(hA, hB, hb, hQ, hS, hR, hq, hr, hidxbx, hlbx, hubx,
                     NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                     NULL, NULL, NULL, NULL, NULL, qp_in);

    print_ocp_qp_in(qp_in);

    void *opts = ocp_qp_opts_create(config, &dims);

    ocp_qp_out *qp_out = ocp_qp_out_create(config, &dims);

    // TODO(dimitris): only have N2 in one place!!
    // printf("N2 in config = %d\n", config->N2);
    // printf("N2 in opts = %d\n", ((ocp_qp_partial_condensing_opts *)(((ocp_qp_partial_condensing_solver_opts *)opts)->pcond_opts))->N2);

    ocp_qp_solver *qp_solver = ocp_qp_create(config, &dims, opts);

    int acados_return = ocp_qp_solve(qp_solver, qp_in, qp_out);

    if (acados_return != ACADOS_SUCCESS)
        return -1;

    print_ocp_qp_out(qp_out);

    /************************************************
     * compute inf norm of residuals
     ************************************************/

    double res[4];
    ocp_qp_inf_norm_residuals(&dims, qp_in, qp_out, res);
    printf("\ninf norm res: %e, %e, %e, %e\n\n", res[0], res[1], res[2], res[3]);

    // ocp_qp_info *info = (ocp_qp_info *)qp_out->misc;
    // print_ocp_qp_info(info);

    ocp_qp_config_free(config);
    ocp_qp_opts_free(opts);
    ocp_qp_in_free(qp_in);
    ocp_qp_out_free(qp_out);
    // TODO(oj): use free from C interface
    free(qp_solver);
}
