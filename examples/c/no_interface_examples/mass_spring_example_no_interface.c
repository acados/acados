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
#include <stdlib.h>

// acados
#include <acados/utils/print.h>
#include <acados/ocp_qp/ocp_qp_partial_condensing_solver.h>
#include <acados/ocp_qp/ocp_qp_full_condensing_solver.h>
#include <acados/ocp_qp/ocp_qp_hpipm.h>
#ifdef ACADOS_WITH_HPMPC
#include <acados/ocp_qp/ocp_qp_hpmpc.h>
#endif
#ifdef ACADOS_WITH_QPDUNES
#include <acados/ocp_qp/ocp_qp_qpdunes.h>
#endif
#include <acados/dense_qp/dense_qp_hpipm.h>
#ifdef ACADOS_WITH_QPOASES
#include <acados/dense_qp/dense_qp_qpoases.h>
#endif
#ifdef ACADOS_WITH_QORE
#include <acados/dense_qp/dense_qp_qore.h>
#endif

// mass spring
ocp_qp_dims *create_ocp_qp_dims_mass_spring(int N, int nx_, int nu_, int nb_, int ng_, int ngN);
ocp_qp_in *create_ocp_qp_in_mass_spring(void *config, ocp_qp_dims *dims);



#ifndef ACADOS_WITH_QPDUNES
#define ELIMINATE_X0
#endif
#define GENERAL_CONSTRAINT_AT_TERMINAL_STAGE

#define NREP 100

int main() {
    printf("\n");
    printf("\n");
    printf("\n");
    printf(" mass spring example: acados ocp_qp solvers\n");
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

    int num_of_stages_equal_to_zero = 4;  // number of states to be enforced to zero at last stage

    #ifdef GENERAL_CONSTRAINT_AT_TERMINAL_STAGE
    int ngN = num_of_stages_equal_to_zero;  // number of general constraints at the last stage
    #else
    int ngN = 0;
    #endif

	ocp_qp_dims *qp_dims = create_ocp_qp_dims_mass_spring(N, nx_, nu_, nb_, ng_, ngN);

    /************************************************
     * ocp qp solvers
     ************************************************/

    // choose values for N2 in partial condensing solvers
    int num_N2_values = 3;
    int N2_values[3] = {15, 5, 3};

	// number of tested solver configurations
	int ii_max = 7;

    for (int ii = 0; ii < ii_max; ii++)
    {

        for (int jj = 0; jj < num_N2_values; jj++)
        {
            int N2 = N2_values[jj];

			int config_size = ocp_qp_xcond_solver_config_calculate_size();
			void *config_mem = malloc(config_size);
			ocp_qp_xcond_solver_config *solver_config = ocp_qp_xcond_solver_config_assign(config_mem);

			int solver_opts_size;
			void *solver_opts_mem;
			void *solver_opts = NULL;
			ocp_qp_xcond_solver_opts *xcond_solver_opts;
			ocp_qp_hpipm_opts *hpipm_opts;
			ocp_qp_partial_condensing_opts *part_cond_opts;
			int solver_mem_size;
			void *solver_mem_mem;
			void *solver_mem;
			int solver_work_size;
			void *solver_work;
			dense_qp_hpipm_opts *dense_hpipm_opts;

			switch(ii)
			{
				case 0: // HPIPM
                    printf("\nHPIPM\n\n");

					// config
					ocp_qp_partial_condensing_solver_config_initialize_default(solver_config);
					ocp_qp_partial_condensing_config_initialize_default(solver_config->xcond);
					ocp_qp_hpipm_config_initialize_default(solver_config->qp_solver);

					// opts
					solver_opts_size = solver_config->opts_calculate_size(solver_config, qp_dims);
// 				printf("\nopts size = %d\n", solver_opts_size);
					solver_opts_mem = malloc(solver_opts_size);
					solver_opts = solver_config->opts_assign(solver_config, qp_dims, solver_opts_mem);
					solver_config->opts_initialize_default(solver_config, qp_dims, solver_opts);
					xcond_solver_opts = solver_opts;
					hpipm_opts = xcond_solver_opts->qp_solver_opts;
					hpipm_opts->hpipm_opts->iter_max = 30;
					part_cond_opts = xcond_solver_opts->xcond_opts;
					part_cond_opts->N2 = N;
					break;

				case 1: // PARTIAL_CONDENSING_HPIPM
                    printf("\nPARTIAL_CONDENSING_HPIPM, N2 = %d\n\n", N2);

					// config
					ocp_qp_partial_condensing_solver_config_initialize_default(solver_config);
					ocp_qp_partial_condensing_config_initialize_default(solver_config->xcond);
					ocp_qp_hpipm_config_initialize_default(solver_config->qp_solver);
					// solver_config->N2 = N2;

					// opts
					solver_opts_size = solver_config->opts_calculate_size(solver_config, qp_dims);
// 				printf("\nopts size = %d\n", solver_opts_size);
					solver_opts_mem = malloc(solver_opts_size);
					solver_opts = solver_config->opts_assign(solver_config, qp_dims, solver_opts_mem);
					solver_config->opts_initialize_default(solver_config, qp_dims, solver_opts);
					xcond_solver_opts = solver_opts;
					hpipm_opts = xcond_solver_opts->qp_solver_opts;
					hpipm_opts->hpipm_opts->iter_max = 30;
					part_cond_opts = xcond_solver_opts->xcond_opts;
					part_cond_opts->N2 = N;

					break;

				case 2: // FULL_CONDENSING_HPIPM
                    printf("\nFULL_CONDENSING_HPIPM\n\n");

					// config
					ocp_qp_full_condensing_solver_config_initialize_default(solver_config);
					ocp_qp_full_condensing_config_initialize_default(solver_config->xcond);
					dense_qp_hpipm_config_initialize_default(solver_config->qp_solver);

					// opts
					solver_opts_size = solver_config->opts_calculate_size(solver_config, qp_dims);
// 				printf("\nopts size = %d\n", solver_opts_size);
					solver_opts_mem = malloc(solver_opts_size);
					solver_opts = solver_config->opts_assign(solver_config, qp_dims, solver_opts_mem);
					solver_config->opts_initialize_default(solver_config, qp_dims, solver_opts);
					xcond_solver_opts = solver_opts;
					dense_hpipm_opts = xcond_solver_opts->qp_solver_opts;
					dense_hpipm_opts->hpipm_opts->iter_max = 30;

					break;

				case 3: // FULL_CONDENSING_QPOASES
#ifdef ACADOS_WITH_QPOASES
                    printf("\nFULL_CONDENSING_QPOASES\n\n");

					// config
					ocp_qp_full_condensing_solver_config_initialize_default(solver_config);
					ocp_qp_full_condensing_config_initialize_default(solver_config->xcond);
					dense_qp_qpoases_config_initialize_default(solver_config->qp_solver);

					// opts
					solver_opts_size = solver_config->opts_calculate_size(solver_config, qp_dims);
// 				printf("\nopts size = %d\n", solver_opts_size);
					solver_opts_mem = malloc(solver_opts_size);
					solver_opts = solver_config->opts_assign(solver_config, qp_dims, solver_opts_mem);
					solver_config->opts_initialize_default(solver_config, qp_dims, solver_opts);
// 				xcond_solver_opts = solver_opts;
// 				dense_hpipm_opts = xcond_solver_opts->qp_solver_opts;
// 				dense_hpipm_opts->hpipm_opts->iter_max = 30;
#endif

					break;

				case 4: // FULL_CONDENSING_QORE
#ifdef ACADOS_WITH_QORE
                    printf("\nFULL_CONDENSING_QORE\n\n");

					// config
					ocp_qp_full_condensing_solver_config_initialize_default(solver_config);
					ocp_qp_full_condensing_config_initialize_default(solver_config->xcond);
					dense_qp_qore_config_initialize_default(solver_config->qp_solver);

					// opts
					solver_opts_size = solver_config->opts_calculate_size(solver_config, qp_dims);
// 				printf("\nopts size = %d\n", solver_opts_size);
					solver_opts_mem = malloc(solver_opts_size);
					solver_opts = solver_config->opts_assign(solver_config, qp_dims, solver_opts_mem);
					solver_config->opts_initialize_default(solver_config, qp_dims, solver_opts);
// 				xcond_solver_opts = solver_opts;
// 				dense_hpipm_opts = xcond_solver_opts->qp_solver_opts;
// 				dense_hpipm_opts->hpipm_opts->iter_max = 30;
#endif

					break;

				case 5: // PARTIAL_CONDENSING_HPMPC
#ifdef ACADOS_WITH_HPMPC
                    printf("\nPARTIAL_CONDENSING_HPMPC, N2 = %d\n\n", N2);

					// config
					ocp_qp_partial_condensing_solver_config_initialize_default(solver_config);
					ocp_qp_partial_condensing_config_initialize_default(solver_config->xcond);
					ocp_qp_hpmpc_config_initialize_default(solver_config->qp_solver);

					// opts
					solver_opts_size = solver_config->opts_calculate_size(solver_config, qp_dims);
// 				printf("\nopts size = %d\n", solver_opts_size);
					solver_opts_mem = malloc(solver_opts_size);
					solver_opts = solver_config->opts_assign(solver_config, qp_dims, solver_opts_mem);
					solver_config->opts_initialize_default(solver_config, qp_dims, solver_opts);
					xcond_solver_opts = solver_opts;
					part_cond_opts = xcond_solver_opts->xcond_opts;
					part_cond_opts->N2 = N;
// 				hpipm_opts = xcond_solver_opts->qp_solver_opts;
// 				hpipm_opts->hpipm_opts->iter_max = 30;
#endif

					break;

				case 6: // PARTIAL_CONDENSING_QPDUNES
#ifdef ACADOS_WITH_QPDUNES
                    printf("\nPARTIAL_CONDENSING_QPDUNES, N2 = %d\n\n", N2);

					// config
					ocp_qp_partial_condensing_solver_config_initialize_default(solver_config);
					ocp_qp_partial_condensing_config_initialize_default(solver_config->xcond);
					ocp_qp_qpdunes_config_initialize_default(solver_config->qp_solver);

					// opts
					solver_opts_size = solver_config->opts_calculate_size(solver_config, qp_dims);
// 				printf("\nopts size = %d\n", solver_opts_size);
					solver_opts_mem = malloc(solver_opts_size);
					solver_opts = solver_config->opts_assign(solver_config, qp_dims, solver_opts_mem);
					solver_config->opts_initialize_default(solver_config, qp_dims, solver_opts);
					xcond_solver_opts = solver_opts;
					part_cond_opts = xcond_solver_opts->xcond_opts;
					part_cond_opts->N2 = N;
// 				hpipm_opts = xcond_solver_opts->qp_solver_opts;
// 				hpipm_opts->hpipm_opts->iter_max = 30;
#endif

					break;

				default:
					printf("\nWrong solver name!!!\n\n");
			}

			// break if there is no solver
#ifndef ACADOS_WITH_QPOASES
			if (ii==3)
				break;
#endif
#ifndef ACADOS_WITH_QORE
			if (ii==4)
				break;
#endif
#ifndef ACADOS_WITH_HPMPC
			if (ii==5)
				break;
#endif
#ifndef ACADOS_WITH_QPDUNES
			if (ii==6)
				break;
#endif

			// update after user-defined opts
			solver_config->opts_update(solver_config, qp_dims, solver_opts);

			/************************************************
			 * ocp qp in
			 ************************************************/

			ocp_qp_in *qp_in = create_ocp_qp_in_mass_spring(solver_config, qp_dims);

			/************************************************
			 * ocp qp out
			 ************************************************/

			int qp_out_size = ocp_qp_out_calculate_size(solver_config, qp_dims);
			void *qp_out_mem = malloc(qp_out_size);
			ocp_qp_out *qp_out = ocp_qp_out_assign(solver_config, qp_dims, qp_out_mem);


            /************************************************
            * memory
            ************************************************/

			solver_mem_size = solver_config->memory_calculate_size(solver_config, qp_dims, solver_opts);
// 		printf("\nmem size = %d\n", solver_mem_size);
			solver_mem_mem = malloc(solver_mem_size);
			solver_mem = solver_config->memory_assign(solver_config, qp_dims, solver_opts, solver_mem_mem);

            /************************************************
            * workspace
            ************************************************/

			solver_work_size = solver_config->workspace_calculate_size(solver_config, qp_dims, solver_opts);
// 		printf("\nwork size = %d\n", solver_work_size);
			solver_work = malloc(solver_work_size);

            /************************************************
            * solve
            ************************************************/

            for (int rep = 0; rep < NREP; rep++)
            {
				solver_config->evaluate(solver_config, qp_in, qp_out, solver_opts, solver_mem, solver_work);
			}

            /************************************************
            * info
            ************************************************/

            ocp_qp_info *info = (ocp_qp_info *)qp_out->misc;
            ocp_qp_info min_info;

			min_info.num_iter = info->num_iter;
			min_info.total_time = info->total_time;
			min_info.condensing_time = info->condensing_time;
			min_info.solve_QP_time = info->solve_QP_time;
			min_info.interface_time = info->interface_time;

            /************************************************
            * print solution
            ************************************************/

// 		 print_ocp_qp_out(qp_out);

            /************************************************
            * compute residuals
            ************************************************/

			int res_size = ocp_qp_res_calculate_size(qp_dims);
// 		printf("\nres size = %d\n", res_size);
			void *res_mem = malloc(res_size);
			ocp_qp_res *qp_res = ocp_qp_res_assign(qp_dims, res_mem);

			int res_work_size = ocp_qp_res_workspace_calculate_size(qp_dims);
// 		printf("\nres work size = %d\n", res_work_size);
			void *res_work_mem = malloc(res_work_size);
			ocp_qp_res_ws *res_ws = ocp_qp_res_workspace_assign(qp_dims, res_work_mem);

            ocp_qp_res_compute(qp_in, qp_out, qp_res, res_ws);

            /************************************************
             * print residuals
             ************************************************/

// 		 print_ocp_qp_res(qp_res);

            /************************************************
            * compute infinity norm of residuals
            ************************************************/

            double res[4];
            ocp_qp_res_compute_nrm_inf(qp_res, res);
            double max_res = 0.0;
            for (int ii = 0; ii < 4; ii++)
                max_res = (res[ii] > max_res) ? res[ii] : max_res;

            /************************************************
            * print stats
            ************************************************/

            printf("\ninf norm res: %e, %e, %e, %e\n", res[0], res[1], res[2], res[3]);

            print_ocp_qp_info(&min_info);

            /************************************************
            * free memory
            ************************************************/

			// TODO
			free(config_mem);



			// just one run if there is not partial condensing
			if ((ii==0) | (ii==2) | (ii==3) | (ii==4))
				break;

		}

	}

/************************************************
* return
************************************************/

	printf("\nsuccess!\n\n");

    return 0;
}

