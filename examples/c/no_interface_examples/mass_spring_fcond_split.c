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


// external
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

// acados_c
#include "acados_c/dense_qp_interface.h"
#include "acados_c/ocp_qp_interface.h"
#include "acados_c/condensing_interface.h"

// acados
#include "acados/ocp_qp/ocp_qp_common_frontend.h"
#include "acados/utils/timing.h"

// blasfeo
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"

#define NREP 100
#define ELIMINATE_X0

// mass spring
ocp_qp_dims *create_ocp_qp_dims_mass_spring(int N, int nx_, int nu_, int nb_, int ng_, int ngN);
ocp_qp_in *create_ocp_qp_in_mass_spring(void *config, ocp_qp_dims *dims);

int main()
{
    printf("\n");
    printf("\n");
    printf("\n");
    printf(" acados + full condensing + dense solver + expansion\n");
    printf("\n");
    printf("\n");
    printf("\n");

    /************************************************
     * set up dimensions
     ************************************************/

    int nx_ = 8;   // number of states (it has to be even for the mass-spring system test problem)

    int nu_ = 3;   // number of inputs (controllers) (it has to be at least 1 and
                   // at most nx_/2 for the mass-spring system test problem)

    int N = 15;    // horizon length
    int nb_ = 11;  // number of box constrained inputs and states
    int ng_ = 0;   // 4;  // number of general constraints

    #ifdef GENERAL_CONSTRAINT_AT_TERMINAL_STAGE
    int num_of_stages_equal_to_zero = 4;  // number of states to be enforced to zero at last stage
    int ngN = num_of_stages_equal_to_zero;
    #else
    int ngN = 0;
    #endif

	ocp_qp_dims *qp_dims = create_ocp_qp_dims_mass_spring(N, nx_, nu_, nb_, ng_, ngN);

    /************************************************
     * ocp qp in/out
     ************************************************/

    ocp_qp_in *qp_in = create_ocp_qp_in_mass_spring(NULL, qp_dims);
    ocp_qp_out *qp_out = ocp_qp_out_create(NULL, qp_dims);

    /************************************************
    * dense qp in/out
    ************************************************/

    dense_qp_dims ddims;
    compute_dense_qp_dims(qp_in->dim, &ddims);

    dense_qp_in *qpd_in = dense_qp_in_create(NULL, &ddims);
    dense_qp_out *qpd_out = dense_qp_out_create(NULL, &ddims);

    /************************************************
    * condensing
    ************************************************/

    condensing_plan cond_plan;
    cond_plan.condensing_type = FULL_CONDENSING;

    ocp_qp_xcond_config *cond_config = ocp_qp_condensing_config_create(&cond_plan);

    ocp_qp_full_condensing_opts *cond_opts = ocp_qp_condensing_opts_create(cond_config, qp_in->dim);

    condensing_module *cond_module = ocp_qp_condensing_create(cond_config, qp_in->dim, cond_opts);

    /************************************************
    * dense ipm
    ************************************************/

    dense_qp_solver_plan plan;
    plan.qp_solver = DENSE_QP_HPIPM;

    qp_solver_config *config = dense_qp_config_create(&plan);

    void *dopts = dense_qp_opts_create(config, &ddims);

    dense_qp_solver *qp_solver = dense_qp_create(config, &ddims, dopts);

	int acados_return;  // 0 normal; 1 max iter

    acados_timer timer;
    acados_tic(&timer);

	for(int rep = 0; rep < NREP; rep++)
    {
        ocp_qp_condense(cond_module, qp_in, qpd_in);

        acados_return = dense_qp_solve(qp_solver, qpd_in, qpd_out);

        if (acados_return != 0)
            printf("error with dense qp solution\n");

        ocp_qp_expand(cond_module, qpd_out, qp_out);
    }

    real_t time = acados_toc(&timer)/NREP;

    /************************************************
    * extract solution
    ************************************************/

    ocp_qp_dims *dims = qp_in->dim;

    colmaj_ocp_qp_out *sol;
    void *memsol = malloc(colmaj_ocp_qp_out_calculate_size(dims));
    assign_colmaj_ocp_qp_out(dims, &sol, memsol);
    convert_ocp_qp_out_to_colmaj(qp_out, sol);

    /************************************************
     * compute infinity norm of residuals
     ************************************************/

    double res[4];
    ocp_qp_inf_norm_residuals(qp_dims, qp_in, qp_out, res);

    double max_res = 0.0;
    for (int ii = 0; ii < 4; ii++)
        max_res = (res[ii] > max_res) ? res[ii] : max_res;

    assert(max_res <= 1e6*ACADOS_EPS && "The largest KKT residual greater than 1e6*ACADOS_EPS");

    /************************************************
    * print solution and stats
    ************************************************/

    int *nx = qp_dims->nx;
    int *nu = qp_dims->nu;
    int *nb = qp_dims->nb;
    int *ng = qp_dims->ng;

    printf("\nu = \n");
    for (int ii = 0; ii < N; ii++) d_print_mat(1, nu[ii], sol->u[ii], 1);

    printf("\nx = \n");
    for (int ii = 0; ii <= N; ii++) d_print_mat(1, nx[ii], sol->x[ii], 1);

    printf("\npi = \n");
    for (int ii = 0; ii < N; ii++) d_print_mat(1, nx[ii+1], sol->pi[ii], 1);

    printf("\nlam = \n");
    for (int ii = 0; ii <= N; ii++) d_print_mat(1, 2*nb[ii]+2*ng[ii], sol->lam[ii], 1);

    printf("\ninf norm res: %e, %e, %e, %e\n", res[0], res[1], res[2], res[3]);

    dense_qp_info *info = (dense_qp_info *)qpd_out->misc;

    printf("\nSolution time for %d IPM iterations, averaged over %d runs: %5.2e seconds\n\n\n",
        info->num_iter, NREP, time);

    /************************************************
    * free memory
    ************************************************/

    free(qp_dims);
    free(qp_in);
    free(qpd_in);
    free(qp_out);
    free(qpd_out);
    free(sol);
    free(dopts);
    free(config);
    free(qp_solver);

    free(cond_opts);
    free(cond_config);
    free(cond_module);

	return 0;
}
