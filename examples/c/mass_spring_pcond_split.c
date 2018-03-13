/*
 *    This file is part of acados.
 *
 *    acados is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    acados is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with acados; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

// external
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

// acados_c
#include "acados_c/ocp_qp_interface.h"
#include "acados_c/legacy_create.h"
#include "acados_c/options.h"

// acados
#include "acados/ocp_qp/ocp_qp_partial_condensing.h"
#include "acados/ocp_qp/ocp_qp_common_frontend.h"

// blasfeo
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"

#define NREP 100
#define ELIMINATE_X0

// mass spring
ocp_qp_dims *create_ocp_qp_dims_mass_spring(int N, int nx_, int nu_, int nb_, int ng_, int ngN);
ocp_qp_in *create_ocp_qp_in_mass_spring(void *config, int N, int nx_, int nu_, int nb_, int ng_, int ngN);

int main()
{
    printf("\n");
    printf("\n");
    printf("\n");
    printf(" acados + partial condensing + dense solver + expansion\n");
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

    ocp_qp_in *qp_in = create_ocp_qp_in_mass_spring(NULL, N, nx_, nu_, nb_, ng_, ngN);
    ocp_qp_out *qp_out = ocp_qp_out_create(NULL, qp_dims);


    /************************************************
    * partial condensing opts & memory
    ************************************************/

    // TODO(dimitris): rename

    ocp_qp_partial_condensing_opts *pcond_opts =
        ocp_qp_partial_condensing_create_arguments(qp_in->dim);

    pcond_opts->N2 = 4;

    ocp_qp_partial_condensing_memory *pcond_memory =
        ocp_qp_partial_condensing_create_memory(qp_in->dim, pcond_opts);

    /************************************************
    * partially condensed qp in/out
    ************************************************/

    // TODO(dimitris): can't I do the same for fcond?
    ocp_qp_in *qpp_in = ocp_qp_in_create(NULL, pcond_opts->pcond_dims);
    ocp_qp_out *qpp_out = ocp_qp_out_create(NULL, pcond_opts->pcond_dims);

    /************************************************
    * sparse ipm
    ************************************************/

    ocp_qp_solver_plan plan2;
    plan2.qp_solver = PARTIAL_CONDENSING_HPIPM;  // UUUUUPSSSSS

    ocp_qp_xcond_solver_config *config = ocp_qp_config_create(plan2);

    void *popts = ocp_qp_opts_create(config,  pcond_opts->pcond_dims);

    // NOTE(dimitris): NOT TO DO A SECOND PCOND!
    set_option_int(popts, "hpipm.N2", pcond_opts->pcond_dims->N);

    ocp_qp_solver *qp_solver = ocp_qp_create(config, pcond_opts->pcond_dims, popts);

    int acados_return;

    acados_timer timer;
    acados_tic(&timer);

	for(int rep = 0; rep < NREP; rep++)
    {
        ocp_qp_partial_condensing(qp_in, qpp_in, pcond_opts, pcond_memory, NULL);

        acados_return = ocp_qp_solve(qp_solver, qpp_in, qpp_out);

        if (acados_return != 0)
            printf("error with ocp qp solution\n");

        ocp_qp_partial_expansion(qpp_out, qp_out, pcond_opts, pcond_memory, NULL);
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

    ocp_qp_info *info = (ocp_qp_info *)qpp_out->misc;

    printf("\nSolution time for %d IPM iterations, averaged over %d runs: %5.2e seconds\n\n\n",
        info->num_iter, NREP, time);

    /************************************************
    * free memory
    ************************************************/

    free(qp_dims);
    free(qp_in);
    free(qpp_in);
    free(qpp_out);
    free(sol);
    free(popts);
    free(config);
    free(qp_solver);
    free(pcond_memory);
    free(pcond_opts);

	return 0;
}
