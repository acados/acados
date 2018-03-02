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

#include "acados/ocp_qp/ocp_qp_hpmpc.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_blas.h"
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_v_aux_ext_dep.h"
#include "hpmpc/include/c_interface.h"
#include "hpmpc/include/mpc_solvers.h"
#include "hpmpc/include/lqcp_solvers.h"
#include "hpmpc/include/mpc_aux.h"

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"
#include "acados/utils/mem.h"


int ocp_qp_hpmpc_calculate_args_size(ocp_qp_dims *dims)
{
    int N = dims->N;
    int size = sizeof(ocp_qp_hpmpc_args);
    size += (N + 1) * sizeof(*(((ocp_qp_hpmpc_args *)0)->ux0));
    size += (N + 1) * sizeof(*(((ocp_qp_hpmpc_args *)0)->pi0));
    size += (N + 1) * sizeof(*(((ocp_qp_hpmpc_args *)0)->lam0));
    size += (N + 1) * sizeof(*(((ocp_qp_hpmpc_args *)0)->t0));
    for (int i = 0; i <= N; i++)
    {
        size += (dims->nu[i] + dims->nx[i]) * sizeof(**(((ocp_qp_hpmpc_args *)0)->ux0));
        size += (2 * dims->nb[i] + 2 * dims->ng[i]) * sizeof(**(((ocp_qp_hpmpc_args *)0)->lam0));
        size += (2 * dims->nb[i] + 2 * dims->ng[i]) * sizeof(**(((ocp_qp_hpmpc_args *)0)->t0));
        if (i > 0)  // TODO(dimitris) p0 is not loaded later
            size += (dims->nx[i]) * sizeof(**(((ocp_qp_hpmpc_args *)0)->pi0));
    }
    size += 5 * sizeof(*(((ocp_qp_hpmpc_args *)0)->inf_norm_res));
    size = (size + 63) / 64 * 64;   // make multiple of typical cache line size
    size += 1 * 64;                 // align once to typical cache line size
    return size;
}



void *ocp_qp_hpmpc_assign_args(ocp_qp_dims *dims, void *raw_memory)
{
    ocp_qp_hpmpc_args *args;
    char *c_ptr = (char *) raw_memory;
    args = (ocp_qp_hpmpc_args *) c_ptr;
    // ocp_qp_hpmpc_args **args = (ocp_qp_hpmpc_args **) args_;
    int N = dims->N;

    args->N = dims->N;
    args->M = dims->N;
    args->N2 = dims->N;


    c_ptr += sizeof(ocp_qp_hpmpc_args);

    args->ux0 = (real_t **) c_ptr;
    c_ptr += (N + 1) * sizeof(*(((ocp_qp_hpmpc_args *)0)->ux0));

    args->pi0 = (real_t **) c_ptr;
    c_ptr += (N + 1) * sizeof(*(((ocp_qp_hpmpc_args *)0)->pi0));

    args->lam0 = (real_t **) c_ptr;
    c_ptr += (N + 1) * sizeof(*(((ocp_qp_hpmpc_args *)0)->lam0));

    args->t0 = (real_t **) c_ptr;
    c_ptr += (N + 1) * sizeof(*(((ocp_qp_hpmpc_args *)0)->t0));

    // align memory to typical cache line size
    size_t s_ptr = (size_t) c_ptr;
    s_ptr = (s_ptr + 63) / 64 * 64;
    c_ptr = (char *) s_ptr;

    for (int i = 0; i <= N; i++) {
        args->ux0[i] = (real_t *) c_ptr;
        c_ptr += (dims->nu[i] + dims->nx[i]) * sizeof(**(((ocp_qp_hpmpc_args *)0)->ux0));
    }
    for (int i = 1; i <= N; i++) {
        args->pi0[i] = (real_t *) c_ptr;
        c_ptr += dims->nx[i] * sizeof(**(((ocp_qp_hpmpc_args *)0)->pi0));
    }
    for (int i = 0; i <= N; i++) {
        args->lam0[i] = (real_t *) c_ptr;
        c_ptr += (2 * dims->nb[i] + 2 * dims->ng[i]) * sizeof(**(((ocp_qp_hpmpc_args *)0)->lam0));
    }
    for (int i = 0; i <= N; i++) {
        args->t0[i] = (real_t *) c_ptr;
        c_ptr += (2 * dims->nb[i] + 2 * dims->ng[i]) * sizeof(**(((ocp_qp_hpmpc_args *)0)->t0));
    }

    // TODO(dimitris): is this even written somehwere?
    args->inf_norm_res = (real_t *) c_ptr;
    c_ptr += 5 * sizeof(*(((ocp_qp_hpmpc_args *)0)->inf_norm_res));
    return (void *)args;
}



void ocp_qp_hpmpc_initialize_default_args(void *args_)
{
    ocp_qp_hpmpc_args *args = (ocp_qp_hpmpc_args *)args_;
    args->tol = 1e-8;
    args->max_iter = 100;
    args->mu0 = 1e3;
    args->warm_start = 0;
    args->alpha_min = 1e-8;
}



int ocp_qp_hpmpc_calculate_memory_size(ocp_qp_dims *dims, void *args_)
{
    ocp_qp_hpmpc_args *args = (ocp_qp_hpmpc_args*) args_;

	int N = dims->N;
    int N2 = args->N2;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *nb = dims->nb;
    int *ng = dims->ng;

    int ii;
    int ws_size = sizeof(ocp_qp_hpmpc_memory);

    ws_size += 6*args->max_iter*sizeof(double);  // stats

    ws_size += (N+1)*sizeof(struct blasfeo_dvec);  //hpi
	for (ii = 0; ii <= N; ii++)
        ws_size += blasfeo_memsize_dvec(nx[ii]);  //hpi

    ws_size += hpmpc_d_ip_ocp_hard_tv_work_space_size_bytes_noidxb(N, nx, nu, nb, dims->nbx, dims->nbu, ng, N2);

    // TODO(dimitris): only calculate sizes below when partial tightenging is used
    #if 0
    int ii;
    int M = args->M;

	for ( ii = 0; ii < N; ii++ ) {
		ws_size += sizeof(double)*(nu[ii]+nx[ii]+1)*(nu[ii]+nx[ii]);  // L
		ws_size += sizeof(double)*(nu[ii]+nx[ii]);  // dux
        ws_size += 3*sizeof(double)*(2*nb[ii]+2*ng[ii]);  // dlam, dt, lamt
        ws_size += sizeof(double)*(2*nb[ii]+2*ng[ii]);  // tinv
        ws_size += 2*sizeof(double)*(nb[ii]+ng[ii]);  // Qx, qx
		ws_size += sizeof(double)*(nx[ii+1]);  // Pb
	}
    // TODO(dimitris): put in loop
    ii = N;
	ws_size += sizeof(double)*(nu[ii]+nx[ii]+1)*(nu[ii]+nx[ii]);  // L
    ws_size += sizeof(double)*(nu[ii]+nx[ii]);  // dux
    ws_size += 3*sizeof(double)*(2*nb[ii]+2*ng[ii]);  // dlam, dt, lamt
    ws_size += sizeof(double)*(2*nb[ii]+2*ng[ii]);  // tinv
    ws_size += 2*sizeof(double)*(nb[ii]+ng[ii]);  // Qx, qx
    ws_size += sizeof(double)*(nx[ii+1]);  // Pb

	ws_size += d_back_ric_rec_work_space_size_bytes_libstr(N, nx, nu, nb, ng);
    ws_size += 2*sizeof(double)*(nx[M]+1)*nx[M];  // LxM, PpM
    #endif

    ws_size += 2*64;

    return ws_size;
}



void *ocp_qp_hpmpc_assign_memory(ocp_qp_dims *dims, void *args_, void *raw_memory)
{
    ocp_qp_hpmpc_args *args = (ocp_qp_hpmpc_args*) args_;

    // char pointer
    char *c_ptr = (char *)raw_memory;

    ocp_qp_hpmpc_memory *mem = (ocp_qp_hpmpc_memory *) c_ptr;

    c_ptr += sizeof(ocp_qp_hpmpc_memory);

    mem->hpi = (struct blasfeo_dvec *) c_ptr;
    c_ptr += (dims->N+1)*sizeof(struct blasfeo_dvec);

    mem->stats = (double *)c_ptr;
    c_ptr += 6*args->max_iter*sizeof(double);

    align_char_to(64, &c_ptr);

    // TODO(dimitris): PUT ANY PARTIAL_TIGHTENING DMATS HERE (WITH IF M < N)

	for (int ii = 0; ii <= dims->N; ii++)
        assign_blasfeo_dvec_mem(dims->nx[ii], &mem->hpi[ii], &c_ptr);

    align_char_to(64, &c_ptr);

    mem->hpmpc_work = (void *) c_ptr;

    // TODO(dimitris): add assert, move hpmpc mem to workspace?
    return raw_memory;
}



int ocp_qp_hpmpc_calculate_workspace_size(ocp_qp_dims *dims, void *args_)
{
    return 0;
}



int ocp_qp_hpmpc(ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *args_, void *mem_, void *work_)
{
    ocp_qp_hpmpc_args *hpmpc_args = (ocp_qp_hpmpc_args*) args_;
    ocp_qp_hpmpc_memory *mem = (ocp_qp_hpmpc_memory*) mem_;

    ocp_qp_info *info = (ocp_qp_info *) qp_out->misc;
    acados_timer tot_timer, qp_timer, interface_timer;
    acados_tic(&tot_timer);

    int ii;

    // extract input struct members
    int N = qp_in->dim->N;
    int *nx = qp_in->dim->nx;
    int *nu = qp_in->dim->nu;
    int *nb = qp_in->dim->nb;
    int *ng = qp_in->dim->ng;

    int hpmpc_status = -1;

    int M = hpmpc_args->M;

    // extract args from struct
    double mu_tol = hpmpc_args->tol;
    int k_max = hpmpc_args->max_iter;
    double mu0 = hpmpc_args->mu0;
    int warm_start = hpmpc_args->warm_start;

    //  other solver arguments
	int compute_mult = 1;
    int kk = -1;  // actual number of iterations

    // IPM constants
    // int kk_avg;
    double alpha_min = 1e-8;  // TODO(dimitris): why not in opts?

    acados_tic(&interface_timer);

    #if 0
	struct blasfeo_dmat *hsmatdummy = NULL;
    struct blasfeo_dvec *hsvecdummy = NULL;
    struct blasfeo_dvec hsQx[N+1];
    struct blasfeo_dvec hsqx[N+1];
    struct blasfeo_dvec hstinv[N+1];
    struct blasfeo_dvec hsdux[N+1];

	struct blasfeo_dvec hsdlam[N+1];  // to be checked
	struct blasfeo_dvec hsdt[N+1];
	struct blasfeo_dvec hslamt[N+1];  // to be checked

	// partial tightening-specific
    struct blasfeo_dvec hsPb[N+1];
    struct blasfeo_dmat hsL[N+1];
    //    struct blasfeo_dmat hsLxt[N+1];
    struct blasfeo_dmat hsric_work_mat[2];

	for (ii = 0; ii < N; ii++)
    {
		// partial tightening-specific
		blasfeo_create_dmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], &hsL[ii], ptr_memory);
		ptr_memory += (&hsL[ii])->memsize;
	}
    ii = N;
	// partial tightening-specific
	blasfeo_create_dmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], &hsL[ii], ptr_memory);
	ptr_memory += (&hsL[ii])->memsize;
    #endif

    for (int ii = 0; ii <= N; ++ii)
    {
		// temporarily invert sign of upper bounds
		blasfeo_dvecsc(nb[ii], -1.0, &qp_in->d[ii], nb[ii] + ng[ii]);
		blasfeo_dvecsc(ng[ii], -1.0, &qp_in->d[ii], 2*nb[ii] + ng[ii]);
    }

	// dvec loop
	for ( ii = 0; ii < N; ii++ )
    {
        // TODO(dimitris): why do we _always_ take init. from args? what about warmstart?
        // copy initialization, multipliers and slacks from hpmpc_args
        blasfeo_pack_dvec(nu[ii]+nx[ii], hpmpc_args->ux0[ii], &qp_out->ux[ii], 0);
        blasfeo_pack_dvec(2*nb[ii]+2*ng[ii], hpmpc_args->lam0[ii], &qp_out->lam[ii], 0);
        blasfeo_pack_dvec(2*nb[ii]+2*ng[ii], hpmpc_args->t0[ii], &qp_out->t[ii], 0);
        // TODO(dimitris): pi0 missing

        #if 0
        // initialize hsdux to primal input later usx will be subtracted
        blasfeo_create_dvec(nu[ii]+nx[ii], &hsdux[ii], ptr_memory);
        blasfeo_pack_dvec(nu[ii]+nx[ii], hpmpc_args->ux0[ii], &hsdux[ii], 0);
        ptr_memory += (&hsdux[ii])->memsize;

        blasfeo_create_dvec(2*nb[ii]+2*ng[ii], &hstinv[ii], ptr_memory);
        ptr_memory += (&hstinv[ii])->memsize;
        blasfeo_create_dvec(nb[ii]+ng[ii], &hsQx[ii], ptr_memory);
        ptr_memory += (&hsQx[ii])->memsize;
        blasfeo_create_dvec(nb[ii]+ng[ii], &hsqx[ii], ptr_memory);
        ptr_memory += (&hsqx[ii])->memsize;

        blasfeo_create_dvec(2*nb[ii]+2*ng[ii], &hsdlam[ii], ptr_memory);
        ptr_memory += (&hsdlam[ii])->memsize;
        blasfeo_create_dvec(2*nb[ii]+2*ng[ii], &hsdt[ii], ptr_memory);
        ptr_memory += (&hsdt[ii])->memsize;
        blasfeo_create_dvec(2*nb[ii]+2*ng[ii], &hslamt[ii], ptr_memory);
        ptr_memory += (&hslamt[ii])->memsize;

		// partial tightening specific
		blasfeo_create_dvec(nx[ii+1], &hsPb[ii+1], ptr_memory);
		ptr_memory += (&hsPb[ii+1])->memsize;
        #endif
    }

    ii = N;
	// dvec loop
	// initialize hsdux to primal input later usx will be subtracted
    #if 0
    blasfeo_create_dvec(nu[ii]+nx[ii], &hsdux[ii], ptr_memory);
    blasfeo_pack_dvec(nu[ii]+nx[ii], hpmpc_args->ux0[ii], &hsdux[ii], 0);
    ptr_memory += (&hsdux[ii])->memsize;
    #endif

    // TODO(dimitris): also use nu[ii] for consistency
    blasfeo_pack_dvec(nx[ii], hpmpc_args->ux0[ii], &qp_out->ux[ii], 0);
    blasfeo_pack_dvec(2*nb[ii]+2*ng[ii], hpmpc_args->lam0[ii], &qp_out->lam[ii], 0);
    blasfeo_pack_dvec(2*nb[ii]+2*ng[ii], hpmpc_args->t0[ii], &qp_out->t[ii], 0);

    #if 0
    blasfeo_create_dvec(2*nb[ii]+2*ng[ii], &hslamt[ii], ptr_memory);
    ptr_memory += (&hslamt[ii])->memsize;

    blasfeo_create_dvec(2*nb[ii]+2*ng[ii], &hstinv[ii], ptr_memory);
    ptr_memory += (&hstinv[ii])->memsize;
    blasfeo_create_dvec(nb[ii]+ng[ii], &hsQx[ii], ptr_memory);
    ptr_memory += (&hsQx[ii])->memsize;
    blasfeo_create_dvec(nb[ii]+ng[ii], &hsqx[ii], ptr_memory);
    ptr_memory += (&hsqx[ii])->memsize;
    blasfeo_create_dvec(2*nb[ii]+2*ng[ii], &hsdlam[ii], ptr_memory);
    ptr_memory += (&hsdlam[ii])->memsize;
    blasfeo_create_dvec(2*nb[ii]+2*ng[ii], &hsdt[ii], ptr_memory);
    ptr_memory += (&hsdt[ii])->memsize;

	real_t sigma_mu = hpmpc_args->sigma_mu;

    int nuM;
    int nbM;

    struct blasfeo_dmat sLxM;
    struct blasfeo_dmat sPpM;

	// partial tightening specific

    blasfeo_create_dmat(nx[M]+1, nx[M], &sLxM, ptr_memory);
    ptr_memory += (&sLxM)->memsize;
    blasfeo_create_dmat(nx[M]+1, nx[M], &sPpM, ptr_memory);
    ptr_memory += (&sPpM)->memsize;

    struct blasfeo_dmat hstmpmat0;

    // riccati work space
    void *work_ric;

    work_ric = ptr_memory;
    ptr_memory+=d_back_ric_rec_work_space_size_bytes_libstr(N, nx, nu, nb, ng);
    #endif

    info->interface_time = acados_toc(&interface_timer);
    acados_tic(&qp_timer);

    if (M < N)
    {
        #if 0
        // update cost function matrices and vectors (box constraints)
        d_update_hessian_gradient_mpc_hard_libstr(N-M, &nx[M], &nu[M], &nb[M], &ng[M], \
          &hsd[M], sigma_mu, &hst[M], &hstinv[M], &hslam[M], &hslamt[M], &hsdlam[M], \
          &hsQx[M], &hsqx[M]);

        // backward riccati factorization and solution at the end
        d_back_ric_rec_sv_back_libstr(N-M, &nx[M], &nu[M], &nb[M], &hsidxb[M], &ng[M], \
          0, &hsBAbt[M], hsvecdummy, 1, &hsRSQrq[M], &hsrq[M], &hsDCt[M], &hsQx[M], \
          &hsqx[M], &hsux[M], 1, &hspi[M],  1, &hsPb[M], &hsL[M], work_ric);

        // extract chol factor of [P p; p' *]
        // TODO(Andrea): have m and n !!!!!
        blasfeo_dtrcp_l(nx[M], &hsL[M], nu[M], nu[M], &sLxM, 0, 0);
        blasfeo_dgecp(1, nx[M], &hsL[M], nu[M]+nx[M], nu[M], &sLxM, nx[M], 0);

        // recover [P p; p' *]
        blasfeo_dsyrk_ln_mn(nx[M]+1, nx[M], nx[M], 1.0, &sLxM, 0, 0, &sLxM, 0, 0, 0.0,
          &sPpM, 0, 0, &sPpM, 0, 0);

        // backup stage M
        nuM = nu[M];
        nbM = nb[M];
        hstmpmat0 = hsRSQrq[M];

        // update new terminal cost
        nu[M] = 0;
        nb[M] = 0;
        hsRSQrq[M] = sPpM;
        hsux[M].pa += nuM;

        // IPM at the beginning
        hpmpc_status = d_ip2_res_mpc_hard_libstr(&kk, k_max, mu0, mu_tol, alpha_min,
          warm_start, stat, M, nx, nu, nb, hsidxb, ng, hsBAbt, hsRSQrq, hsDCt,
          hsd, hsux, compute_mult, hspi, hslam, hst, ptr_memory);  // recover original stage M

        nu[M] = nuM;
        nb[M] = nbM;
        hsRSQrq[M] = hstmpmat0;
        hsux[M].pa -= nuM;

        // forward riccati solution at the end
        d_back_ric_rec_sv_forw_libstr(N-M, &nx[M], &nu[M], &nb[M], &hsidxb[M], &ng[M],
          0, &hsBAbt[M], hsvecdummy, 1, &hsRSQrq[M], &hsrq[M], hsmatdummy,
          &hsQx[M], &hsqx[M], &hsux[M], 1, &hspi[M], 1, &hsPb[M], &hsL[M],
          hsric_work_mat);

        // compute alpha, dlam and dt
        real_t alpha = 1.0;
        // compute primal step hsdux for stages M to N
        real_t *temp_p1, *temp_p2;
        for (int i = M; i <= N; i++) {
          // hsdux is initialized to be equal to hpmpc_args.ux0
          temp_p1 = hsdux[i].pa;
          temp_p2 = hsux[i].pa;  // hsux[i].pa;
          for (int j = 0; j < nx[i]+nu[i]; j++) temp_p1[j]= - temp_p1[j] + temp_p2[j];
        }

        d_compute_alpha_mpc_hard_libstr(N-M, &nx[M], &nu[M], &nb[M], &hsidxb[M],
          &ng[M], &alpha, &hst[M], &hsdt[M], &hslam[M], &hsdlam[M], &hslamt[M],
          &hsdux[M], &hsDCt[M], &hsd[M]);

        // update stages M to N
        double mu_scal = 0.0;
        d_update_var_mpc_hard_libstr(N-M, &nx[M], &nu[M], &nb[M], &ng[M],
          &sigma_mu, mu_scal, alpha, &hsux[M], &hsdux[M], &hst[M], &hsdt[M], &hslam[M],
          &hsdlam[M], &hspi[M], &hspi[M]);

        // !!!! TODO(Andrea): equality multipliers are not being updated! Need to
        // define and compute hsdpi (see function prototype).
    #endif
    } else {
        // IPM at the beginning
        hpmpc_status = d_ip2_res_mpc_hard_libstr(&kk, k_max, mu0, mu_tol, alpha_min,
          warm_start, mem->stats, N, nx, nu, nb, qp_in->idxb, ng, qp_in->BAbt, qp_in->RSQrq, qp_in->DCt,
          qp_in->d, qp_out->ux, compute_mult, mem->hpi, qp_out->lam, qp_out->t, mem->hpmpc_work);  // recover original stage M
    }

    info->solve_QP_time = acados_toc(&qp_timer);
    acados_tic(&interface_timer);

	hpmpc_args->out_iter = kk;  // TODO(dimitris): obsolete

    // copy result to qp_out
    for ( ii = 0; ii < N; ii++ )
        blasfeo_dveccp(nx[ii+1], &mem->hpi[ii+1], 0, &qp_out->pi[ii], 0);

	// restore sign of upper bounds
	for(int jj = 0; jj <=N; jj++)
    {
		blasfeo_dvecsc(nb[jj], -1.0, &qp_in->d[jj], nb[jj] + ng[jj]);
		blasfeo_dvecsc(ng[jj], -1.0, &qp_in->d[jj], 2*nb[jj] + ng[jj]);
	}

    info->interface_time += acados_toc(&interface_timer);
    info->total_time = acados_toc(&tot_timer);
    info->num_iter = kk;

    int acados_status = hpmpc_status;
    if (hpmpc_status == 0) acados_status = ACADOS_SUCCESS;
    if (hpmpc_status == 1) acados_status = ACADOS_MAXITER;
    if (hpmpc_status == 2) acados_status = ACADOS_MINSTEP;
    return acados_status;
}



void ocp_qp_hpmpc_config_initialize_default(void *config_)
{

	ocp_qp_solver_config *config = config_;

	config->opts_calculate_size = &ocp_qp_hpmpc_calculate_args_size;
	config->opts_assign = &ocp_qp_hpmpc_assign_args;
	config->opts_initialize_default = &ocp_qp_hpmpc_initialize_default_args;
	config->memory_calculate_size = &ocp_qp_hpmpc_calculate_memory_size;
	config->memory_assign = &ocp_qp_hpmpc_assign_memory;
	config->workspace_calculate_size = &ocp_qp_hpmpc_calculate_workspace_size;
	config->fun = &ocp_qp_hpmpc;

	return;
}
