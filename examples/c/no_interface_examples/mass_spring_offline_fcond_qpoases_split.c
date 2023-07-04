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
#include "acados_c/dense_qp_interface.h"
#include "acados_c/ocp_qp_interface.h"
#include "acados_c/condensing_interface.h"

// acados
#include "acados/ocp_qp/ocp_qp_full_condensing.h"
#include "acados/ocp_qp/ocp_qp_common_frontend.h"
#include "acados/utils/timing.h"

// blasfeo
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_blas.h"

// NOTE(nielsvd): required to cast memory etc. should go.
#include "acados/ocp_qp/ocp_qp_full_condensing_solver.h"
#include "acados/dense_qp/dense_qp_hpipm.h"
#include "acados/dense_qp/dense_qp_qpoases.h"

#define NREP 1  // TODO(dimitris): are timings valid for NREP > 1?
#define ELIMINATE_X0

// TODO(dimitris): when switching this flag, first run gives sometimes a segfault
#define OFFLINE_CONDENSING 1
#define BLASFEO_CHOLESKY 1

// mass spring
ocp_qp_dims *create_ocp_qp_dims_mass_spring(int N, int nx_, int nu_, int nb_, int ng_, int ngN);
ocp_qp_in *create_ocp_qp_in_mass_spring(void *config, ocp_qp_dims *dims);

int main() {
    printf("\n");
    printf("\n");
    printf("\n");
	if(OFFLINE_CONDENSING == 1)
		printf(" acados + offline condensing + qpoases + expansion\n");
	else
		printf(" acados + online condensing + qpoases + expansion\n");
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
    * dense qpoases
    ************************************************/

    dense_qp_solver_plan plan;
    plan.qp_solver = DENSE_QP_QPOASES;

    qp_solver_config *config = dense_qp_config_create(&plan);

    void *dopts = dense_qp_opts_create(config, &ddims);

	dense_qp_qpoases_opts *args = (dense_qp_qpoases_opts *)dopts;

	if (BLASFEO_CHOLESKY == 1)
		args->use_precomputed_cholesky = 1;

	if (OFFLINE_CONDENSING == 1)
		args->hotstart = 1;

	dense_qp_solver *qp_solver = dense_qp_create(config, &ddims, dopts);

	int acados_return;

	int nvd = qpd_in->dim->nv;

	ocp_qp_condense(cond_module, qp_in, qpd_in);

	struct blasfeo_dmat sR;
	blasfeo_allocate_dmat(nvd, nvd, &sR);

	if(OFFLINE_CONDENSING == 1)
    {
        cond_opts->cond_hess = 0;
		cond_opts->expand_dual_sol = 0;

        ocp_qp_full_condensing_opts_update(qp_in->dim, cond_opts);
        ocp_qp_full_condensing_opts_update(qp_in->dim, cond_opts);

		// cholesky factorization of H
		dense_qp_qpoases_memory *qpoases_solver_mem = (dense_qp_qpoases_memory *)qp_solver->mem;
		blasfeo_dpotrf_l(nvd, qpd_in->Hv, 0, 0, &sR, 0, 0);

		// fill in upper triangular of R
		blasfeo_dtrtr_l(nvd, &sR, 0, 0, &sR, 0, 0);

		// extract R
		blasfeo_unpack_dmat(nvd, nvd, &sR, 0, 0, qpoases_solver_mem->R, nvd);
	}

    ocp_qp_condense(cond_module, qp_in, qpd_in);
    // TODO(dimitris): what is this (and the one above) doing here?
    ocp_qp_expand(cond_module, qpd_out, qp_out);

	acados_timer timer;
    acados_tic(&timer);

	for (int rep = 0; rep < NREP; rep++)
    {
		if (OFFLINE_CONDENSING == 0 && BLASFEO_CHOLESKY == 1)
        {
			// cholesky factorization of H
			dense_qp_qpoases_memory *qpoases_solver_mem = (dense_qp_qpoases_memory *)qp_solver->mem;
			blasfeo_dpotrf_l(nvd, qpd_in->Hv, 0, 0, &sR, 0, 0);

			// fill in upper triangular of R
			blasfeo_dtrtr_l(nvd, &sR, 0, 0, &sR, 0, 0);

			// extract R
			blasfeo_unpack_dmat(nvd, nvd, &sR, 0, 0, qpoases_solver_mem->R, nvd);
		}
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

    assert(max_res <= ACADOS_EPS && "The largest KKT residual greater than ACADOS_EPS");

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

    // NOTE(nielsvd): how can we improve/generalize this?
    dense_qp_qpoases_memory *mem = (dense_qp_qpoases_memory *)(qp_solver->mem);

    printf("\ninf norm res: %e, %e, %e, %e\n", res[0], res[1], res[2], res[3]);

    printf("\nSolution time for %d iterations, averaged over %d runs: %5.2e seconds\n\n\n",
        mem->nwsr, NREP, time);

    /************************************************
    * free memory
    ************************************************/

    free(qp_in);
    free(qpd_in);
    free(qp_out);
    free(qpd_out);
    free(sol);
    free(qp_solver);
    free(dopts);

    free(cond_opts);
    free(cond_config);
    free(cond_module);

	return 0;
}
