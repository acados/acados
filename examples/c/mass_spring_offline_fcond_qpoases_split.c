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
#include <acados_c/dense_qp.h>
#include <acados_c/legacy_create.h>
// blasfeo
#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_blas.h"
// acados
#include <acados/ocp_qp/ocp_qp_full_condensing.h>
// NOTE(nielsvd): required to cast memory etc. should go.
#include <acados/ocp_qp/ocp_qp_full_condensing_solver.h>
#include <acados/dense_qp/dense_qp_hpipm.h>
#include <acados/dense_qp/dense_qp_qpoases.h>

#define NREP 1
#define ELIMINATE_X0 
#define OFFLINE_CONDENSING 1
#define BLASFEO_CHOLESKY 0

#include "./mass_spring.c"

int main() {
    printf("\n");
    printf("\n");
    printf("\n");
	if(OFFLINE_CONDENSING == 1) {
		printf(" acados + offline condensing + qpoases + expansion\n");
	} else {
		printf(" acados + online condensing + qpoases + expansion\n");
	}
    printf("\n");
    printf("\n");
    printf("\n");

    /************************************************
    * ocp qp
    ************************************************/

    ocp_qp_in *qp_in = create_ocp_qp_in_mass_spring();

    ocp_qp_dims *qp_dims = qp_in->dim;

    int N = qp_dims->N;
    int *nx = qp_dims->nx;
    int *nu = qp_dims->nu;
    int *nb = qp_dims->nb;
    int *ng = qp_dims->ng;

    /************************************************
    * dense qp
    ************************************************/

    dense_qp_dims ddims;
    compute_dense_qp_dims(qp_in->dim, &ddims);

    dense_qp_in *qpd_in = create_dense_qp_in(&ddims);

    ocp_qp_full_condensing_args *cond_args = ocp_qp_full_condensing_create_arguments(qp_in->dim, NULL);
	ocp_qp_full_condensing_initialize_default_args(cond_args);
	ocp_qp_full_condensing_memory *cond_memory = ocp_qp_full_condensing_create_memory(qp_in->dim, cond_args);

    /************************************************
    * ocp qp solution
    ************************************************/

    ocp_qp_out *qp_out = create_ocp_qp_out(qp_in->dim);

    /************************************************
    * dense sol
    ************************************************/

    dense_qp_out *qpd_out = create_dense_qp_out(&ddims);

    /************************************************
    * dense qpoases
    ************************************************/

    dense_qp_solver_config config;
	config.qp_solver = DENSE_QP_QPOASES;

    dense_qp_solver_fcn_ptrs *fcn_ptrs = create_dense_qp_solver_fcn_ptrs(&config, &ddims);

    void *argd = dense_qp_create_args(fcn_ptrs, &ddims);

	dense_qp_qpoases_args *args = (dense_qp_qpoases_args *)argd; 
	
	if (BLASFEO_CHOLESKY == 1) {
		args->use_precomputed_cholesky = 1;
	}

	if (OFFLINE_CONDENSING == 1) {
		args->hotstart = 1; 
	}
	
	dense_qp_solver *qp_solver = dense_qp_create(fcn_ptrs, &ddims, argd);

    free(fcn_ptrs);
	
	int acados_return;  // 0 normal; 1 max iter

	int nvd = qpd_in->dim->nv;
	
	ocp_qp_full_condensing(qp_in, qpd_in, cond_args, cond_memory, NULL);
	dense_qp_qpoases_memory *qpoases_mem = (dense_qp_qpoases_memory *)qp_solver->mem;
	
	struct blasfeo_dmat sR;
	blasfeo_allocate_dmat(nvd, nvd, &sR);
	
	if(OFFLINE_CONDENSING == 1) {
		cond_args->condense_rhs_only = 1; 
		cond_args->expand_primal_sol_only = 1; 
		
		// cholesky factorization of H
		dense_qp_qpoases_memory *qpoases_solver_mem = (dense_qp_qpoases_memory *)qp_solver->mem;
		blasfeo_dpotrf_l(nvd, qpd_in->Hv, 0, 0, &sR, 0, 0);

		// fill in upper triangular of R 
		blasfeo_dtrtr_l(nvd, &sR, 0, 0, &sR, 0, 0); 

		// extract R 
		blasfeo_unpack_dmat(nvd, nvd, &sR, 0, 0, qpoases_solver_mem->R, nvd); 
	} 

	ocp_qp_full_condensing(qp_in, qpd_in, cond_args, cond_memory, NULL);  
	
	ocp_qp_full_expansion(qpd_out, qp_out, cond_args, cond_memory, NULL);
	
	acados_timer timer;
    acados_tic(&timer);

	for(int rep = 0; rep < NREP; rep++) {
		if(OFFLINE_CONDENSING == 0 && BLASFEO_CHOLESKY == 1) {
			// cholesky factorization of H
			dense_qp_qpoases_memory *qpoases_solver_mem = (dense_qp_qpoases_memory *)qp_solver->mem;
			blasfeo_dpotrf_l(nvd, qpd_in->Hv, 0, 0, &sR, 0, 0);

			// fill in upper triangular of R
			blasfeo_dtrtr_l(nvd, &sR, 0, 0, &sR, 0, 0); 

			// extract R
			blasfeo_unpack_dmat(nvd, nvd, &sR, 0, 0, qpoases_solver_mem->R, nvd); 
		}

        ocp_qp_full_condensing(qp_in, qpd_in, cond_args, cond_memory, NULL); 

        acados_return = dense_qp_solve(qp_solver, qpd_in, qpd_out);

        ocp_qp_full_expansion(qpd_out, qp_out, cond_args, cond_memory, NULL);
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
    * compute residuals
    ************************************************/

    ocp_qp_res *qp_res = create_ocp_qp_res(dims);
    ocp_qp_res_ws *res_ws = create_ocp_qp_res_ws(dims);
    compute_ocp_qp_res(qp_in, qp_out, qp_res, res_ws);

    /************************************************
    * compute infinity norm of residuals
    ************************************************/

    double res[4];
    compute_ocp_qp_res_nrm_inf(qp_res, res);
    double max_res = 0.0;
    for (int ii = 0; ii < 4; ii++) max_res = (res[ii] > max_res) ? res[ii] : max_res;
    // assertion switched off when using primal-only expansion (no multipliers computed)
	if (cond_args->expand_primal_sol_only == 0) {
		assert(max_res <= 1e6*ACADOS_EPS && "The largest KKT residual greater than 1e6*ACADOS_EPS");  
	}

    /************************************************
    * print solution and stats
    ************************************************/

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
    free(qp_res);
    free(res_ws);
    free(argd);
    free(cond_memory);
    free(cond_args);

	return 0;
}
