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
#include <string.h>

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
#include "acados/utils/mem.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"
#include "acados/utils/mem.h"


int ocp_qp_hpmpc_calculate_args_size(ocp_qp_dims *dims, void *submodules_)
{
    int_t N = dims->N;
    int_t size = sizeof(ocp_qp_hpmpc_args);
    size += (N + 1) * sizeof(*(((ocp_qp_hpmpc_args *)0)->ux0));
    size += (N + 1) * sizeof(*(((ocp_qp_hpmpc_args *)0)->pi0));
    size += (N + 1) * sizeof(*(((ocp_qp_hpmpc_args *)0)->lam0));
    size += (N + 1) * sizeof(*(((ocp_qp_hpmpc_args *)0)->t0));
    for (int_t i = 0; i <= N; i++) {
        size += (dims->nu[i] + dims->nx[i]) * sizeof(**(((ocp_qp_hpmpc_args *)0)->ux0));
        if (i > 0) size += (dims->nx[i]) * sizeof(**(((ocp_qp_hpmpc_args *)0)->pi0));
        size += (2 * dims->nb[i] + 2 * dims->ng[i]) * sizeof(**(((ocp_qp_hpmpc_args *)0)->lam0));
        size += (2 * dims->nb[i] + 2 * dims->ng[i]) * sizeof(**(((ocp_qp_hpmpc_args *)0)->t0));
    }
    size += 5 * sizeof(*(((ocp_qp_hpmpc_args *)0)->inf_norm_res));
    size = (size + 63) / 64 * 64;   // make multiple of typical cache line size
    size += 1 * 64;                 // align once to typical cache line size
    return size;
}



void *ocp_qp_hpmpc_assign_args(ocp_qp_dims *dims, void **submodules_, void *raw_memory)
{
    ocp_qp_hpmpc_args *args;
    char *c_ptr = (char *) raw_memory;
    args = (ocp_qp_hpmpc_args *) c_ptr;
    // ocp_qp_hpmpc_args **args = (ocp_qp_hpmpc_args **) args_;
    int_t N = dims->N;

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

    for (int_t i = 0; i <= N; i++) {
        args->ux0[i] = (real_t *) c_ptr;
        c_ptr += (dims->nu[i] + dims->nx[i]) * sizeof(**(((ocp_qp_hpmpc_args *)0)->ux0));
    }
    for (int_t i = 1; i <= N; i++) {
        args->pi0[i] = (real_t *) c_ptr;
        c_ptr += dims->nx[i] * sizeof(**(((ocp_qp_hpmpc_args *)0)->pi0));
    }
    for (int_t i = 0; i <= N; i++) {
        args->lam0[i] = (real_t *) c_ptr;
        c_ptr += (2 * dims->nb[i] + 2 * dims->ng[i]) * sizeof(**(((ocp_qp_hpmpc_args *)0)->lam0));
    }
    for (int_t i = 0; i <= N; i++) {
        args->t0[i] = (real_t *) c_ptr;
        c_ptr += (2 * dims->nb[i] + 2 * dims->ng[i]) * sizeof(**(((ocp_qp_hpmpc_args *)0)->t0));
    }
    args->inf_norm_res = (real_t *) c_ptr;
    c_ptr += 5 * sizeof(*(((ocp_qp_hpmpc_args *)0)->inf_norm_res));

    // Update submodules pointer
    *submodules_ = NULL;

    return (void *)args;
}



void *ocp_qp_hpmpc_copy_args(ocp_qp_dims *dims, void *raw_memory, void *source_)
{
    ocp_qp_hpmpc_args *source = (ocp_qp_hpmpc_args *)source_;
    ocp_qp_hpmpc_args *dest;

    void *submodules;

    dest = ocp_qp_hpmpc_assign_args(dims, &submodules, raw_memory);

    dest->tol = source->tol;
    dest->max_iter = source->max_iter;
    dest->mu0 = source->mu0;
    dest->warm_start = source->warm_start;
    dest->N2 = source->N2;
    dest->out_iter = source->out_iter;
    dest->sigma_mu = source->sigma_mu;
    dest->N = source->N;
    dest->M = source->M;

    int_t N = dims->N;
    int_t sz;
    for (int_t i = 0; i <= N; i++) {
        sz = (dims->nu[i] + dims->nx[i]) * sizeof(**(((ocp_qp_hpmpc_args *)0)->ux0));
        memcpy(dest->ux0[i], source->ux0[i], sz);
    }
    for (int_t i = 1; i <= N; i++) {
        sz = dims->nx[i] * sizeof(**(((ocp_qp_hpmpc_args *)0)->pi0));
        memcpy(dest->pi0[i], source->pi0[i], sz);
    }
    for (int_t i = 0; i <= N; i++) {
        sz = (2 * dims->nb[i] + 2 * dims->ng[i]) * sizeof(**(((ocp_qp_hpmpc_args *)0)->lam0));
        memcpy(dest->lam0[i], source->lam0[i], sz);
    }
    for (int_t i = 0; i <= N; i++) {
        sz = (2 * dims->nb[i] + 2 * dims->ng[i]) * sizeof(**(((ocp_qp_hpmpc_args *)0)->t0));
        memcpy(dest->t0[i], source->t0[i], sz);
    }
    sz = 5 * sizeof(*(((ocp_qp_hpmpc_args *)0)->inf_norm_res));
    memcpy(dest->inf_norm_res, source->inf_norm_res, sz);

    return (void *)dest;
}



void ocp_qp_hpmpc_initialize_default_args(void *args_)
{
    ocp_qp_hpmpc_args *args = (ocp_qp_hpmpc_args *)args_;
    args->tol = 1e-8;
    args->max_iter = 100;
    args->mu0 = 1e3;
    args->warm_start = 0;

    return;

}



int ocp_qp_hpmpc_calculate_memory_size(ocp_qp_dims *dims, void *args_)
{
    ocp_qp_hpmpc_args *args = (ocp_qp_hpmpc_args*) args_;

    int N = dims->N;
    int *nx = (int *) dims->nx;
    int *nu = (int *) dims->nu;
    int *nb = (int *) dims->nb;
    int *ng = (int *) dims->ng;
    int_t N2 = args->N2;

    int_t ws_size;

    int ii;
    int_t max_ip_iter = args->max_iter;
    ws_size = 8 + 5*max_ip_iter*sizeof(double);
    //        ws_size += 1 * (N + 1) * sizeof(int *);  // hidxb_rev
    for (ii = 0; ii <= N; ii++) {
        ws_size += nb[ii]*sizeof(int);  // hidxb_rev
    }
    ws_size += hpmpc_d_ip_ocp_hard_tv_work_space_size_bytes_noidxb(N, nx, nu, nb, dims->nbx, dims->nbu, ng, N2);

    return ws_size;

}



void *ocp_qp_hpmpc_assign_memory(ocp_qp_dims *dims, void *args_, void *raw_memory)
{
    return raw_memory;
}



int ocp_qp_hpmpc_calculate_workspace_size(ocp_qp_dims *dims, void *args_)
{
    return 0;
}



int ocp_qp_hpmpc(ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *args_, void *mem_, void *work_)
{
    ocp_qp_hpmpc_args *hpmpc_args = (ocp_qp_hpmpc_args*) args_;
    ocp_qp_info *info = (ocp_qp_info *) qp_out->misc;
    acados_timer tot_timer, qp_timer, interface_timer;
    acados_tic(&tot_timer);

    // loop index
    int ii, jj;

    // extract input struct members
    int N = qp_in->dim->N;
    int *nx = (int *) qp_in->dim->nx;
    int *nu = (int *) qp_in->dim->nu;
    int *nb = (int *) qp_in->dim->nb;
    int **hsidxb = (int **) qp_in->idxb;
    int *ng = (int *) qp_in->dim->ng;

    int hpmpc_status = -1;

    int_t M = hpmpc_args->M;
    char *ptr_memory = (char *) mem_;

    // extract args struct members
    double mu_tol = hpmpc_args->tol;
    int k_max = hpmpc_args->max_iter;
    double mu0 = hpmpc_args->mu0;
    int warm_start = hpmpc_args->warm_start;

    //  other solver arguments
    int kk = -1;  // actual number of iterations

    // IPM constants
    // int kk_avg;
    double alpha_min = 1e-8;
    double *stat = (double*)ptr_memory;
    ptr_memory+=sizeof(double)*k_max*6;
	align_char_to(64, &ptr_memory);
	int compute_mult = 1;


    struct blasfeo_dmat *hsmatdummy = NULL;
    struct blasfeo_dvec *hsvecdummy = NULL;

    struct blasfeo_dmat hsBAbt[N+1];
    struct blasfeo_dvec hsb[N+1];
    struct blasfeo_dmat hsRSQrq[N+1];
    struct blasfeo_dvec hsQx[N+1];
    struct blasfeo_dvec hsqx[N+1];
    struct blasfeo_dvec hstinv[N+1];
    struct blasfeo_dvec hsrq[N+1];
    struct blasfeo_dmat hsDCt[N+1];
    struct blasfeo_dvec hsd[N+1];
    struct blasfeo_dvec hsux[N+1];
    struct blasfeo_dvec hsdux[N+1];
    struct blasfeo_dvec hspi[N+1];
    struct blasfeo_dvec hslam[N+1];
    struct blasfeo_dvec hst[N+1];
    struct blasfeo_dvec hsPb[N+1];
    struct blasfeo_dmat hsL[N+1];
    //    struct blasfeo_dmat hsLxt[N+1];
    struct blasfeo_dmat hsric_work_mat[2];

    struct blasfeo_dvec hsdlam[N+1];  // to be checked
    struct blasfeo_dvec hsdt[N+1];
    struct blasfeo_dvec hslamt[N+1];  // to be checked

    acados_tic(&interface_timer);

	// qp_in loop
    for ( ii = 0; ii < N; ii++ ) {

        blasfeo_create_dmat(nu[ii]+nx[ii]+1, nx[ii+1], &hsBAbt[ii], qp_in->BAbt[ii].pA);
        blasfeo_create_dvec(nx[ii+1], &hsb[ii], qp_in->b[ii].pa);
        blasfeo_create_dmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], &hsRSQrq[ii], qp_in->RSQrq[ii].pA);
        blasfeo_create_dvec(nu[ii]+nx[ii], &hsrq[ii], qp_in->rq[ii].pa);
        blasfeo_create_dmat(nu[ii]+nx[ii], ng[ii], &hsDCt[ii], qp_in->DCt[ii].pA);
        blasfeo_create_dvec(2*nb[ii]+2*ng[ii], &hsd[ii], qp_in->d[ii].pa);
		// temporarily invert sign of upper bounds
		blasfeo_dvecsc(nb[ii], -1.0, &hsd[ii], nb[ii] + ng[ii]);
		blasfeo_dvecsc(ng[ii], -1.0, &hsd[ii], 2*nb[ii] + ng[ii]);
    }

    ii = N;
	// ocp_in loop
    blasfeo_create_dmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], &hsRSQrq[ii], qp_in->RSQrq[ii].pA);
    blasfeo_create_dvec(nu[ii]+nx[ii], &hsrq[ii], qp_in->rq[ii].pa);
    blasfeo_create_dmat(nu[ii]+nx[ii], ng[ii], &hsDCt[ii], qp_in->DCt[ii].pA);
    blasfeo_create_dvec(2*nb[ii]+2*ng[ii], &hsd[ii], qp_in->d[ii].pa);
	// temporarily invert sign of upper bounds
	blasfeo_dvecsc(nb[ii], -1.0, &hsd[ii], nb[ii] + ng[ii]);
	blasfeo_dvecsc(ng[ii], -1.0, &hsd[ii], 2*nb[ii] + ng[ii]);

	// dmat loop
	for ( ii = 0; ii < N; ii++ )
		assign_blasfeo_dmat_mem(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], &hsL[ii], &ptr_memory);

    ii = N;
	// dmat loop
	assign_blasfeo_dmat_mem(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], &hsL[ii], &ptr_memory);

	// dvec loop
	for ( ii = 0; ii < N; ii++ ) {

        // initialize hsdux to primal input later usx will be subtracted
        blasfeo_create_dvec(nu[ii]+nx[ii], &hsdux[ii], ptr_memory);
        blasfeo_pack_dvec(nu[ii]+nx[ii], hpmpc_args->ux0[ii], &hsdux[ii], 0);
        ptr_memory += (&hsdux[ii])->memsize;
        blasfeo_create_dvec(nu[ii]+nx[ii], &hsux[ii], ptr_memory);
        blasfeo_pack_dvec(nu[ii]+nx[ii], hpmpc_args->ux0[ii], &hsux[ii], 0);
        ptr_memory += (&hsux[ii])->memsize;

        blasfeo_create_dvec(nx[ii+1], &hspi[ii], ptr_memory);
        ptr_memory += (&hspi[ii])->memsize;

        blasfeo_create_dvec(2*nb[ii]+2*ng[ii], &hslam[ii], ptr_memory);
        // copy multipliers from hpmpc_args
        blasfeo_pack_dvec(2*nb[ii]+2*ng[ii], hpmpc_args->lam0[ii], &hslam[ii], 0);
        ptr_memory += (&hslam[ii])->memsize;

        blasfeo_create_dvec(2*nb[ii]+2*ng[ii], &hst[ii], ptr_memory);
        // copy slacks from hpmpc_args
        blasfeo_pack_dvec(2*nb[ii]+2*ng[ii], hpmpc_args->t0[ii], &hst[ii], 0);
        ptr_memory += (&hst[ii])->memsize;

        blasfeo_create_dvec(nx[ii+1], &hsPb[ii+1], ptr_memory);
        ptr_memory += (&hsPb[ii+1])->memsize;

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
    }

    ii = N;
	// dvec loop
	// initialize hsdux to primal input later usx will be subtracted
    blasfeo_create_dvec(nu[ii]+nx[ii], &hsdux[ii], ptr_memory);
    blasfeo_pack_dvec(nu[ii]+nx[ii], hpmpc_args->ux0[ii], &hsdux[ii], 0);
    ptr_memory += (&hsdux[ii])->memsize;
    blasfeo_create_dvec(nx[ii], &hsux[ii], ptr_memory);
    blasfeo_pack_dvec(nx[ii], hpmpc_args->ux0[ii], &hsux[ii], 0);
    ptr_memory += (&hsux[ii])->memsize;

    blasfeo_create_dvec(nx[ii], &hspi[ii], ptr_memory);  // Andrea: bug?
    ptr_memory += (&hspi[ii])->memsize;

    blasfeo_create_dvec(2*nb[ii]+2*ng[ii], &hslam[ii], ptr_memory);
    blasfeo_pack_dvec(2*nb[ii]+2*ng[ii], hpmpc_args->lam0[ii], &hslam[ii], 0);
    ptr_memory += (&hslam[ii])->memsize;

    blasfeo_create_dvec(2*nb[ii]+2*ng[ii], &hst[ii], ptr_memory);
    blasfeo_pack_dvec(2*nb[ii]+2*ng[ii], hpmpc_args->t0[ii], &hst[ii], 0);
    ptr_memory += (&hst[ii])->memsize;

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

    blasfeo_create_dmat(nx[M]+1, nx[M], &sLxM, ptr_memory);
    ptr_memory += (&sLxM)->memsize;
    blasfeo_create_dmat(nx[M]+1, nx[M], &sPpM, ptr_memory);
    ptr_memory += (&sPpM)->memsize;

    struct blasfeo_dmat hstmpmat0;

    // riccati work space
    void *work_ric;

    work_ric = ptr_memory;
    ptr_memory+=d_back_ric_rec_work_space_size_bytes_libstr(N, nx, nu, nb, ng);

    info->interface_time = acados_toc(&interface_timer);
    acados_tic(&qp_timer);

    if (M < N){
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
        for (int_t i = M; i <= N; i++) {
          // hsdux is initialized to be equal to hpmpc_args.ux0
          temp_p1 = hsdux[i].pa;
          temp_p2 = hsux[i].pa;  // hsux[i].pa;
          for (int_t j = 0; j < nx[i]+nu[i]; j++) temp_p1[j]= - temp_p1[j] + temp_p2[j];
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
    } else {
        // IPM at the beginning
        hpmpc_status = d_ip2_res_mpc_hard_libstr(&kk, k_max, mu0, mu_tol, alpha_min,
          warm_start, stat, N, nx, nu, nb, hsidxb, ng, hsBAbt, hsRSQrq, hsDCt,
          hsd, hsux, compute_mult, hspi, hslam, hst, ptr_memory);  // recover original stage M
    }

    info->solve_QP_time = acados_toc(&qp_timer);
    acados_tic(&interface_timer);

    // copy result to qp_out
    for ( ii = 0; ii < N; ii++ ) {
        blasfeo_dveccp(nx[ii] + nu[ii], &hsux[ii], 0, &qp_out->ux[ii], 0);
        blasfeo_dveccp(nx[ii], &hspi[ii], 0, &qp_out->pi[ii], 0);
        blasfeo_dveccp(2*nb[ii]+2*ng[ii], &hslam[ii], 0, &qp_out->lam[ii], 0);
        blasfeo_dveccp(2*nb[ii]+2*ng[ii], &hst[ii], 0, &qp_out->t[ii], 0);

    }

    ii = N;
    blasfeo_dveccp(nx[ii] + nu[ii], &hsux[ii], 0, &qp_out->ux[ii], 0);
    blasfeo_dveccp(2*nb[ii]+2*ng[ii], &hslam[ii], 0, &qp_out->lam[ii], 0);
    blasfeo_dveccp(2*nb[ii]+2*ng[ii], &hst[ii], 0, &qp_out->t[ii], 0);

    hpmpc_args->out_iter = kk;

	// restore sign of upper bounds
	for(int jj = 0; jj <=N; jj++) {
		blasfeo_dvecsc(nb[jj], -1.0, &hsd[jj], nb[jj] + ng[jj]);
		blasfeo_dvecsc(ng[jj], -1.0, &hsd[jj], 2*nb[jj] + ng[jj]);
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
