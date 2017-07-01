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

#include "acados/ocp_qp/ocp_qp_hpipm.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_blas.h"
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_v_aux_ext_dep.h"

#include "hpipm/include/hpipm_d_ocp_qp.h"
#include "hpipm/include/hpipm_d_ocp_qp_sol.h"
#include "hpipm/include/hpipm_d_ocp_qp_ipm_hard.h"

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/types.h"

//void d_print_e_mat(int m, int n, double *a, int lda);



int ocp_qp_hpipm_calculate_workspace_size(ocp_qp_in *qp_in, ocp_qp_hpipm_args *args) {
//
	return 0;
}



int ocp_qp_hpipm_calculate_memory_size(ocp_qp_in *qp_in, ocp_qp_hpipm_args *args) {
//
	int ii;

	int N = qp_in->N;
	int *nx = (int *) qp_in->nx;
	int *nu = (int *) qp_in->nu;
	int *nb = (int *) qp_in->nb;
	int *ng = (int *) qp_in->nc;

	struct d_ocp_qp qp;
	qp.N = N;
	qp.nx = nx;
	qp.nu = nu;
	qp.nb = nb;
	qp.ng = ng;

	struct d_ipm_hard_ocp_qp_arg ipm_arg;
	ipm_arg.iter_max = args->iter_max;

	int size = 0;

	size += d_memsize_ocp_qp(N, nx, nu, nb, ng);
	size += d_memsize_ocp_qp_sol(N, nx, nu, nb, ng);
	size += d_memsize_ipm_hard_ocp_qp(&qp, &ipm_arg);
	size += 4*(N+1)*sizeof(double *); // lam_lb lam_ub lam_lg lam_ug

	size = (size+63)/64*64; // make multipl of typical cache line size
	size += 1*64; // align once to typical cache line size

	return size;

}



void ocp_qp_hpipm_create_memory(ocp_qp_in *qp_in, ocp_qp_hpipm_args *args, ocp_qp_hpipm_memory *hpipm_memory, void *memory) {
//

	// loop indexed
	int ii;


    // extract problem size
    int N = qp_in->N;
    int *nx = (int *) qp_in->nx;
    int *nu = (int *) qp_in->nu;
    int *nb = (int *) qp_in->nb;
    int *ng = (int *) qp_in->nc;


	// char pointer
	char *c_ptr = (char *) memory;


	//
	hpipm_memory->qp = (struct d_ocp_qp *) c_ptr;
	c_ptr += 1*sizeof(struct d_ocp_qp);
	//
	hpipm_memory->qp_sol = (struct d_ocp_qp_sol *) c_ptr;
	c_ptr += 1*sizeof(struct d_ocp_qp_sol);
	//
	hpipm_memory->ipm_arg = (struct d_ipm_hard_ocp_qp_arg *) c_ptr;
	c_ptr += 1*sizeof(struct d_ipm_hard_ocp_qp_arg);
	//
	hpipm_memory->ipm_workspace = (struct d_ipm_hard_ocp_qp_workspace *) c_ptr;
	c_ptr += 1*sizeof(struct d_ipm_hard_ocp_qp_workspace);
	//
	hpipm_memory->hlam_lb = (double **) c_ptr;
	c_ptr += (N+1)*sizeof(double *);
	//
	hpipm_memory->hlam_ub = (double **) c_ptr;
	c_ptr += (N+1)*sizeof(double *);
	//
	hpipm_memory->hlam_lg = (double **) c_ptr;
	c_ptr += (N+1)*sizeof(double *);
	//
	hpipm_memory->hlam_ug = (double **) c_ptr;
	c_ptr += (N+1)*sizeof(double *);


	//
	struct d_ocp_qp *qp = hpipm_memory->qp;
	//
	struct d_ocp_qp_sol *qp_sol = hpipm_memory->qp_sol;
	//
	struct d_ipm_hard_ocp_qp_arg *ipm_arg = hpipm_memory->ipm_arg;
	//
	struct d_ipm_hard_ocp_qp_workspace *ipm_workspace = hpipm_memory->ipm_workspace;


	// align memory to typical cache line size
	size_t s_ptr = (size_t) c_ptr;
	s_ptr = (s_ptr+63)/64*64;
	c_ptr = (char *) s_ptr;


	// ocp qp structure
	d_create_ocp_qp(N, nx, nu, nb, ng, qp, c_ptr);
	c_ptr += qp->memsize;
	// ocp qp sol structure
	d_create_ocp_qp_sol(N, nx, nu, nb, ng, qp_sol, c_ptr);
	c_ptr += qp_sol->memsize;
	// ipm arg structure
	ipm_arg->iter_max = args->iter_max;
	ipm_arg->alpha_min = args->alpha_min;
	ipm_arg->mu_max = args->mu_max;
	ipm_arg->mu0 = args->mu0;
	// ipm workspace structure
	d_create_ipm_hard_ocp_qp(qp, ipm_arg, ipm_workspace, c_ptr);
	c_ptr += ipm_workspace->memsize;


	return;

}



int ocp_qp_hpipm(ocp_qp_in *qp_in, ocp_qp_out *qp_out,
        ocp_qp_hpipm_args *args, ocp_qp_hpipm_memory *memory, void *workspace_) {
//
    // initialize return code
    int acados_status = ACADOS_SUCCESS;

    // loop index
    int ii, jj;

	// extract memory
	double **hlam_lb = memory->hlam_lb;
	double **hlam_ub = memory->hlam_ub;
	double **hlam_lg = memory->hlam_lg;
	double **hlam_ug = memory->hlam_ug;
	struct d_ocp_qp *qp = memory->qp;
	struct d_ocp_qp_sol *qp_sol = memory->qp_sol;
//	struct d_ipm_hard_ocp_qp_arg *ipm_arg = memory->ipm_arg;
	struct d_ipm_hard_ocp_qp_workspace *ipm_workspace = memory->ipm_workspace;

    // extract problem size
    int N = qp_in->N;
    int *nx = (int *) qp_in->nx;
    int *nu = (int *) qp_in->nu;
    int *nb = (int *) qp_in->nb;
    int *ng = (int *) qp_in->nc;

	// extract input data
    double **hA = (double **) qp_in->A;
    double **hB = (double **) qp_in->B;
    double **hb = (double **) qp_in->b;
    double **hQ = (double **) qp_in->Q;
    double **hS = (double **) qp_in->S;
    double **hR = (double **) qp_in->R;
    double **hq = (double **) qp_in->q;
    double **hr = (double **) qp_in->r;
    double **hd_lb = (double **) qp_in->lb;
    double **hd_ub = (double **) qp_in->ub;
    double **hC = (double **) qp_in->Cx;
    double **hD = (double **) qp_in->Cu;
    double **hd_lg = (double **) qp_in->lc;
    double **hd_ug = (double **) qp_in->uc;
    int **hidxb = (int **) qp_in->idxb;

	//  swap x and u in bounds (by updating their indeces)
	int itmp;
	for (ii = 0; ii <= N; ii++) {
		itmp = 0;
		for(jj=0; jj<ii; jj++)
			if(hidxb[ii]==hidxb[jj])
				itmp = 1;
		if(itmp==0) // new physical array
			{
			jj = 0;
			for (; jj < nb[ii]; jj++) {
				if (hidxb[ii][jj] < nx[ii]) {  // state
					hidxb[ii][jj] = hidxb[ii][jj]+nu[ii];
				} else {  // input
					hidxb[ii][jj] = hidxb[ii][jj]-nx[ii];
				}
			}
		}
	}

    // extract output struct members
    double **hx = qp_out->x;
    double **hu = qp_out->u;
    double **hpi = qp_out->pi;
    double **hlam = qp_out->lam;

	//
	for(ii=0; ii<=N; ii++)
		{
		hlam_lb[ii] = hlam[ii];
		hlam_ub[ii] = hlam[ii]+nb[ii];
		hlam_lg[ii] = hlam[ii]+2*nb[ii];
		hlam_ug[ii] = hlam[ii]+2*nb[ii]+ng[ii];
		}

	// ocp qp structure
	d_cvt_colmaj_to_ocp_qp(hA, hB, hb, hQ, hS, hR, hq, hr, hidxb, hd_lb, hd_ub, hC, hD, hd_lg, hd_ug, qp);

	// ocp qp sol structure

	// arg structure // XXX fixed in memory !!!
//	ipm_arg->alpha_min = args->alpha_min;
//	ipm_arg->mu_max = args->mu_max;
//	ipm_arg->iter_max = args->iter_max;
//	ipm_arg->mu0 = args->mu0;

	// ipm structure


	// solve ipm
	d_solve_ipm2_hard_ocp_qp(qp, qp_sol, ipm_workspace);


	// extract solution
	d_cvt_ocp_qp_sol_to_colmaj(qp, qp_sol, hu, hx, hpi, hlam_lb, hlam_ub, hlam_lg, hlam_ug);


	// extract iteration number
	memory->iter = ipm_workspace->iter;


	// compute infinity norm of residuals
	double *inf_norm_res = memory->inf_norm_res;
	double res_tmp;
	//
	double *res_g;
	inf_norm_res[0] = 0;
	for(ii=0; ii<=N; ii++)
		{
		res_g = (ipm_workspace->res_g+ii)->pa;
		for(jj=0; jj<nu[ii]+nx[ii]; jj++)
			{
			res_tmp = fabs(res_g[jj]);
			if(res_tmp>inf_norm_res[0])
				{
				inf_norm_res[0] = res_tmp;
				}
			}
		}
	double *res_b;
	inf_norm_res[1] = 0;
	for(ii=0; ii<N; ii++)
		{
		res_b = (ipm_workspace->res_b+ii)->pa;
		for(jj=0; jj<nx[ii+1]; jj++)
			{
			res_tmp = fabs(res_b[jj]);
			if(res_tmp>inf_norm_res[1])
				{
				inf_norm_res[1] = res_tmp;
				}
			}
		}
	double *res_d;
	inf_norm_res[2] = 0;
	for(ii=0; ii<=N; ii++)
		{
		res_d = (ipm_workspace->res_d_lb+ii)->pa;
		for(jj=0; jj<nb[ii]; jj++)
			{
			res_tmp = fabs(res_d[jj]);
			if(res_tmp>inf_norm_res[2])
				{
				inf_norm_res[2] = res_tmp;
				}
			}
		res_d = (ipm_workspace->res_d_ub+ii)->pa;
		for(jj=0; jj<nb[ii]; jj++)
			{
			res_tmp = fabs(res_d[jj]);
			if(res_tmp>inf_norm_res[2])
				{
				inf_norm_res[2] = res_tmp;
				}
			}
		res_d = (ipm_workspace->res_d_lg+ii)->pa;
		for(jj=0; jj<ng[ii]; jj++)
			{
			res_tmp = fabs(res_d[jj]);
			if(res_tmp>inf_norm_res[2])
				{
				inf_norm_res[2] = res_tmp;
				}
			}
		res_d = (ipm_workspace->res_d_ug+ii)->pa;
		for(jj=0; jj<ng[ii]; jj++)
			{
			res_tmp = fabs(res_d[jj]);
			if(res_tmp>inf_norm_res[2])
				{
				inf_norm_res[2] = res_tmp;
				}
			}
		}
	double *res_m;
	inf_norm_res[3] = 0;
	for(ii=0; ii<=N; ii++)
		{
		res_m = (ipm_workspace->res_m_lb+ii)->pa;
		for(jj=0; jj<nb[ii]; jj++)
			{
			res_tmp = fabs(res_m[jj]);
			if(res_tmp>inf_norm_res[3])
				{
				inf_norm_res[3] = res_tmp;
				}
			}
		res_m = (ipm_workspace->res_m_ub+ii)->pa;
		for(jj=0; jj<nb[ii]; jj++)
			{
			res_tmp = fabs(res_m[jj]);
			if(res_tmp>inf_norm_res[3])
				{
				inf_norm_res[3] = res_tmp;
				}
			}
		res_m = (ipm_workspace->res_m_lg+ii)->pa;
		for(jj=0; jj<ng[ii]; jj++)
			{
			res_tmp = fabs(res_m[jj]);
			if(res_tmp>inf_norm_res[3])
				{
				inf_norm_res[3] = res_tmp;
				}
			}
		res_m = (ipm_workspace->res_m_ug+ii)->pa;
		for(jj=0; jj<ng[ii]; jj++)
			{
			res_tmp = fabs(res_m[jj]);
			if(res_tmp>inf_norm_res[3])
				{
				inf_norm_res[3] = res_tmp;
				}
			}
		}
	inf_norm_res[4] = ipm_workspace->res_mu;


	// max number of iterations
    if (ipm_workspace->iter==args->iter_max)
		acados_status = ACADOS_MAXITER;
	// minimum step length
    if (ipm_workspace->stat[3+(ipm_workspace->iter-1)*5]<args->alpha_min)
		acados_status = ACADOS_MINSTEP;
	

	//  swap (back) x and u in bounds (by updating their indeces)
	for (ii = 0; ii <= N; ii++) {
		itmp = 0;
		for(jj=0; jj<ii; jj++)
			if(hidxb[ii]==hidxb[jj])
				itmp = 1;
		if(itmp==0) // new physical array
			{
			jj = 0;
			for (; jj < nb[ii]; jj++) {
				if (hidxb[ii][jj] >= nu[ii]) {  // state
					hidxb[ii][jj] = hidxb[ii][jj]-nu[ii];
				} else {  // input
					hidxb[ii][jj] = hidxb[ii][jj]+nx[ii];
				}
			}
		}
	}


    // return
    return acados_status;
//
}

