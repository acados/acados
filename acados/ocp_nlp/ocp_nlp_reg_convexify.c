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


#include "acados/ocp_nlp/ocp_nlp_reg_convexify.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "acados/utils/math.h"
#include "acados/utils/mem.h"

#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_blas.h"



/************************************************
 * opts
 ************************************************/

acados_size_t ocp_nlp_reg_convexify_opts_calculate_size(void)
{
    return sizeof(ocp_nlp_reg_convexify_opts);
}



void *ocp_nlp_reg_convexify_opts_assign(void *raw_memory)
{
    return raw_memory;
}



void ocp_nlp_reg_convexify_opts_initialize_default(void *config_, ocp_nlp_reg_dims *dims, void *opts_)
{
    ocp_nlp_reg_convexify_opts *opts = opts_;

    opts->delta = 1e-4;
    opts->epsilon = 1e-4;

    return;
}



void ocp_nlp_reg_convexify_opts_set(void *config_, ocp_nlp_reg_dims *dims, void *opts_, char *field, void* value)
{

    ocp_nlp_reg_convexify_opts *opts = opts_;

    if (!strcmp(field, "delta"))
    {
        double *d_ptr = value;
        opts->delta = *d_ptr;
    }
    else if (!strcmp(field, "epsilon"))
    {
        double *d_ptr = value;
        opts->epsilon = *d_ptr;
    }
    else
    {
        printf("\nerror: field %s not available in ocp_nlp_reg_convexify_opts_set\n", field);
        exit(1);
    }

    return;
}



/************************************************
 * memory
 ************************************************/

acados_size_t ocp_nlp_reg_convexify_calculate_memory_size(void *config_, ocp_nlp_reg_dims *dims, void *opts_)
{

    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *nbx = dims->nbx;
    int *nbu = dims->nbu;
    int *ng = dims->ng;

    int ii;

    int nuM = nu[0];
    for(ii=1; ii<=N; ii++)
    {
        nuM = nu[ii]>nuM ? nu[ii] : nuM;
    }

    int nxM = nx[0];
    for(ii=1; ii<=N; ii++)
    {
        nxM = nx[ii]>nxM ? nx[ii] : nxM;
    }

    int nuxM = nu[0]+nx[0];
    for(ii=1; ii<=N; ii++)
    {
        nuxM = nu[ii]+nx[ii]>nuxM ? nu[ii]+nx[ii] : nuxM;
    }

    int nbgM = nbu[0]+nbx[0]+ng[0];
    for(ii=1; ii<=N; ii++)
    {
        nbgM = nbu[ii]+nbx[ii]+ng[ii]>nbgM ? nbu[ii]+nbx[ii]+ng[ii] : nbgM;
    }

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_reg_convexify_memory);

    size += nuM*nuM*sizeof(double);     // R
    size += nuxM*nuxM*sizeof(double);   // V
    size += 2*nuxM*sizeof(double);      // d e
    size += nuxM*nuxM*sizeof(double);   // reg_hess

    size += (N+1)*sizeof(struct blasfeo_dmat); // original_RSQrq

    size += 1 * 64;

    size += 3*blasfeo_memsize_dmat(nxM, nxM);     // Q_tilde Q_bar delta_eye
    size += blasfeo_memsize_dmat(nuxM, nxM);    // BAQ
    size += blasfeo_memsize_dmat(nuM, nuM);     // L
    size += blasfeo_memsize_dmat(nxM, nuM);     // St_copy

    size += blasfeo_memsize_dvec(nuxM);     // tmp_nuxM
    size += blasfeo_memsize_dvec(nbgM);     // tmp_nbgM

    for (ii=0; ii<=N; ii++)
    {
        size += 2*blasfeo_memsize_dmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii]); // original_RSQrq
    }
	size += blasfeo_memsize_dmat(nuxM, nuxM); // tmp_RSQ

//    size += 2*blasfeo_memsize_dvec(nxM); // grad b2

    // giaf's
    size += (2*(N+1)+N)*sizeof(struct blasfeo_dmat *); // RSQrq BAbt DCt
    size += (3*(N+1)+2*N)*sizeof(struct blasfeo_dvec *); // rq b ux pi lam
	size += (N+1)*sizeof(int *); // idxb

    return size;
}



void *ocp_nlp_reg_convexify_assign_memory(void *config_, ocp_nlp_reg_dims *dims, void *opts_, void *raw_memory)
{

    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *nbx = dims->nbx;
    int *nbu = dims->nbu;
    int *ng = dims->ng;

    int ii;

    int nuM = nu[0];
    for(ii=1; ii<=N; ii++)
    {
        nuM = nu[ii]>nuM ? nu[ii] : nuM;
    }

    int nxM = nx[0];
    for(ii=1; ii<=N; ii++)
    {
        nxM = nx[ii]>nxM ? nx[ii] : nxM;
    }

    int nuxM = nu[0]+nx[0];
    for(ii=1; ii<=N; ii++)
    {
        nuxM = nu[ii]+nx[ii]>nuxM ? nu[ii]+nx[ii] : nuxM;
    }

    int nbgM = nbu[0]+nbx[0]+ng[0];
    for(ii=1; ii<=N; ii++)
    {
        nbgM = nbu[ii]+nbx[ii]+ng[ii]>nbgM ? nbu[ii]+nbx[ii]+ng[ii] : nbgM;
    }


    char *c_ptr = (char *) raw_memory;

    ocp_nlp_reg_convexify_memory *mem = (ocp_nlp_reg_convexify_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_reg_convexify_memory);

    mem->R = (double *) c_ptr;
    c_ptr += nuM*nuM*sizeof(double);

    mem->V = (double *) c_ptr;
    c_ptr += nuxM*nuxM*sizeof(double);

    mem->d = (double *) c_ptr;
    c_ptr += nuxM*sizeof(double);

    mem->e = (double *) c_ptr;
    c_ptr += nuxM*sizeof(double);

    mem->reg_hess = (double *) c_ptr;
    c_ptr += nuxM*nuxM*sizeof(double);

    mem->original_RSQrq = (struct blasfeo_dmat *) c_ptr;
    c_ptr += (N+1)*sizeof(struct blasfeo_dmat);

    mem->RSQrq = (struct blasfeo_dmat **) c_ptr;
    c_ptr += (N+1)*sizeof(struct blasfeo_dmat *);

    mem->BAbt = (struct blasfeo_dmat **) c_ptr;
    c_ptr += N*sizeof(struct blasfeo_dmat *);

    mem->DCt = (struct blasfeo_dmat **) c_ptr;
    c_ptr += (N+1)*sizeof(struct blasfeo_dmat *);

    mem->rq = (struct blasfeo_dvec **) c_ptr;
    c_ptr += (N+1)*sizeof(struct blasfeo_dvec *);

    mem->b = (struct blasfeo_dvec **) c_ptr;
    c_ptr += N*sizeof(struct blasfeo_dvec *);

    mem->ux = (struct blasfeo_dvec **) c_ptr;
    c_ptr += (N+1)*sizeof(struct blasfeo_dvec *);

    mem->pi = (struct blasfeo_dvec **) c_ptr;
    c_ptr += N*sizeof(struct blasfeo_dvec *);

    mem->lam = (struct blasfeo_dvec **) c_ptr;
    c_ptr += (N+1)*sizeof(struct blasfeo_dvec *);

	mem->idxb = (int **) c_ptr;
	c_ptr += (N+1)*sizeof(int *);

    align_char_to(64, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nxM, nxM, &mem->Q_tilde, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nxM, nxM, &mem->Q_bar, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nuxM, nxM, &mem->BAQ, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nuM, nuM, &mem->L, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nxM, nxM, &mem->delta_eye, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nxM, nuM, &mem->St_copy, &c_ptr);

    for (ii=0; ii<=N; ii++)
        assign_and_advance_blasfeo_dmat_mem(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], &mem->original_RSQrq[ii], &c_ptr);


    assign_and_advance_blasfeo_dmat_mem(nuxM, nuxM, &mem->tmp_RSQ, &c_ptr);

    assign_and_advance_blasfeo_dvec_mem(nuxM, &mem->tmp_nuxM, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nbgM, &mem->tmp_nbgM, &c_ptr);

//    assign_and_advance_blasfeo_dvec_mem(nxM, &mem->grad, &c_ptr);
//    assign_and_advance_blasfeo_dvec_mem(nxM, &mem->b2, &c_ptr);

    assert((char *)mem + ocp_nlp_reg_convexify_calculate_memory_size(config_, dims, opts_) >= c_ptr);

    return mem;
}



void ocp_nlp_reg_convexify_memory_set_RSQrq_ptr(ocp_nlp_reg_dims *dims, struct blasfeo_dmat *RSQrq, void *memory_)
{
    ocp_nlp_reg_convexify_memory *memory = memory_;

    int ii;

    int N = dims->N;

    for(ii=0; ii<=N; ii++)
    {
        memory->RSQrq[ii] = RSQrq+ii;
//        blasfeo_print_dmat(nu[ii]+nx[ii], nu[ii]+nx[ii], memory->RSQrq[ii], 0, 0);
    }

    return;
}



void ocp_nlp_reg_convexify_memory_set_rq_ptr(ocp_nlp_reg_dims *dims, struct blasfeo_dvec *rq, void *memory_)
{
    ocp_nlp_reg_convexify_memory *memory = memory_;

    int ii;

    int N = dims->N;
    // int *nx = dims->nx;
    // int *nu = dims->nu;

    for(ii=0; ii<=N; ii++)
    {
        memory->rq[ii] = rq+ii;
//        blasfeo_print_dvec(nu[ii]+nx[ii], memory->rq[ii], 0);
    }

    return;
}



void ocp_nlp_reg_convexify_memory_set_BAbt_ptr(ocp_nlp_reg_dims *dims, struct blasfeo_dmat *BAbt, void *memory_)
{
    ocp_nlp_reg_convexify_memory *memory = memory_;

    int ii;

    int N = dims->N;
    // int *nx = dims->nx;
    // int *nu = dims->nu;

    for(ii=0; ii<N; ii++)
    {
        memory->BAbt[ii] = BAbt+ii;
//        blasfeo_print_dmat(nu[ii]+nx[ii]+1, nx[ii+1], memory->BAbt[ii], 0, 0);
    }

    return;
}



void ocp_nlp_reg_convexify_memory_set_b_ptr(ocp_nlp_reg_dims *dims, struct blasfeo_dvec *b, void *memory_)
{
    ocp_nlp_reg_convexify_memory *memory = memory_;

    int ii;

    int N = dims->N;
    // int *nx = dims->nx;
    // int *nu = dims->nu;

    for(ii=0; ii<N; ii++)
    {
        memory->b[ii] = b+ii;
//        blasfeo_print_dvec(nx[ii=1], memory->b[ii], 0);
    }

    return;
}



void ocp_nlp_reg_convexify_memory_set_idxb_ptr(ocp_nlp_reg_dims *dims, int **idxb, void *memory_)
{
    ocp_nlp_reg_convexify_memory *memory = memory_;

    int ii;

    int N = dims->N;

    for(ii=0; ii<=N; ii++)
    {
        memory->idxb[ii] = idxb[ii];
    }

    return;
}



void ocp_nlp_reg_convexify_memory_set_DCt_ptr(ocp_nlp_reg_dims *dims, struct blasfeo_dmat *DCt, void *memory_)
{
    ocp_nlp_reg_convexify_memory *memory = memory_;

    int ii;

    int N = dims->N;

    for(ii=0; ii<=N; ii++)
    {
        memory->DCt[ii] = DCt+ii;
//        blasfeo_print_dmat(nu[ii]+nx[ii]+1, ng[ii], memory->DCt[ii], 0, 0);
    }

    return;
}



void ocp_nlp_reg_convexify_memory_set_ux_ptr(ocp_nlp_reg_dims *dims, struct blasfeo_dvec *ux, void *memory_)
{
    ocp_nlp_reg_convexify_memory *memory = memory_;

    int ii;

    int N = dims->N;
    // int *nx = dims->nx;
    // int *nu = dims->nu;

    for(ii=0; ii<=N; ii++)
    {
        memory->ux[ii] = ux+ii;
//        blasfeo_print_dvec(nu[ii]+nx[ii], memory->ux[ii], 0);
    }

    return;
}



void ocp_nlp_reg_convexify_memory_set_pi_ptr(ocp_nlp_reg_dims *dims, struct blasfeo_dvec *pi, void *memory_)
{
    ocp_nlp_reg_convexify_memory *memory = memory_;

    int ii;

    int N = dims->N;
    // int *nx = dims->nx;
    // int *nu = dims->nu;

    for(ii=0; ii<N; ii++)
    {
        memory->pi[ii] = pi+ii;
//        blasfeo_print_dvec(nx[ii+1], memory->pi[ii], 0);
    }

    return;
}



void ocp_nlp_reg_convexify_memory_set_lam_ptr(ocp_nlp_reg_dims *dims, struct blasfeo_dvec *lam, void *memory_)
{
    ocp_nlp_reg_convexify_memory *memory = memory_;

    int ii;

    int N = dims->N;
    // int *nbu = dims->nbu;
    // int *nbx = dims->nbx;
    // int *ng = dims->ng;

    for(ii=0; ii<=N; ii++)
    {
        memory->lam[ii] = lam+ii;
//        blasfeo_print_dvec(2*nbu[ii]+2*nbx[ii]+2*ng[ii], memory->lam[ii], 0);
    }

    return;
}



void ocp_nlp_reg_convexify_memory_set(void *config_, ocp_nlp_reg_dims *dims, void *memory_, char *field, void *value)
{

    if(!strcmp(field, "RSQrq_ptr"))
    {
        struct blasfeo_dmat *RSQrq = value;
        ocp_nlp_reg_convexify_memory_set_RSQrq_ptr(dims, RSQrq, memory_);
    }
    else if(!strcmp(field, "rq_ptr"))
    {
        struct blasfeo_dvec *rq = value;
        ocp_nlp_reg_convexify_memory_set_rq_ptr(dims, rq, memory_);
    }
    else if(!strcmp(field, "BAbt_ptr"))
    {
        struct blasfeo_dmat *BAbt = value;
        ocp_nlp_reg_convexify_memory_set_BAbt_ptr(dims, BAbt, memory_);
    }
    else if(!strcmp(field, "b_ptr"))
    {
        struct blasfeo_dvec *b = value;
        ocp_nlp_reg_convexify_memory_set_b_ptr(dims, b, memory_);
    }
    else if(!strcmp(field, "idxb_ptr"))
    {
        int **idxb = value;
        ocp_nlp_reg_convexify_memory_set_idxb_ptr(dims, idxb, memory_);
    }
    else if(!strcmp(field, "DCt_ptr"))
    {
        struct blasfeo_dmat *DCt = value;
        ocp_nlp_reg_convexify_memory_set_DCt_ptr(dims, DCt, memory_);
    }
    else if(!strcmp(field, "ux_ptr"))
    {
        struct blasfeo_dvec *ux = value;
        ocp_nlp_reg_convexify_memory_set_ux_ptr(dims, ux, memory_);
    }
    else if(!strcmp(field, "pi_ptr"))
    {
        struct blasfeo_dvec *pi = value;
        ocp_nlp_reg_convexify_memory_set_pi_ptr(dims, pi, memory_);
    }
    else if(!strcmp(field, "lam_ptr"))
    {
        struct blasfeo_dvec *lam = value;
        ocp_nlp_reg_convexify_memory_set_lam_ptr(dims, lam, memory_);
    }
    else
    {
        printf("\nerror: field %s not available in ocp_nlp_reg_convexify_set\n", field);
        exit(1);
    }

    return;
}



/************************************************
 * functions
 ************************************************/

// NOTE this only considers the case of (dynamcs) equality constraints (no inequality constraints)
// TODO inequality constraints case
void ocp_nlp_reg_convexify_regularize_hessian(void *config, ocp_nlp_reg_dims *dims, void *opts_, void *mem_)
{
    ocp_nlp_reg_convexify_memory *mem = mem_;
    ocp_nlp_reg_convexify_opts *opts = opts_;

    int ii, jj;

    int *nx = dims->nx;
    int *nu = dims->nu;
    int N = dims->N;

#if 0
	for(ii=0; ii<=N; ii++)
	{
        blasfeo_drowin(nu[ii]+nx[ii], 1.0, mem->rq[ii], 0, mem->RSQrq[ii], nu[ii]+nx[ii], 0);
        blasfeo_dgecp(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], mem->RSQrq[ii], 0, 0, &mem->original_RSQrq[ii], 0, 0);
	}
	return;
#endif

    double delta = opts->delta;

    // Algorithm 6 from Verschueren2017

	blasfeo_drowin(nu[N]+nx[N], 1.0, mem->rq[N], 0, mem->RSQrq[N], nu[N]+nx[N], 0);

    blasfeo_dgecp(nu[N]+nx[N]+1, nu[N]+nx[N], mem->RSQrq[N], 0, 0, &mem->original_RSQrq[N], 0, 0);

    // TODO regularize R at last stage if needed !!!
	// TODO fix for nu[N]>0 !!!!!!!!!!
    blasfeo_dgese(nx[N], nx[N], 0.0, &mem->delta_eye, 0, 0);
    blasfeo_ddiare(nx[N], delta, &mem->delta_eye, 0, 0);
    blasfeo_dgecp(nx[N], nx[N], &mem->delta_eye, 0, 0, &mem->Q_tilde, 0, 0);
    blasfeo_dgecp(nx[N], nx[N], mem->RSQrq[N], nu[N], nu[N], &mem->Q_bar, 0, 0);
    blasfeo_dgecp(nx[N], nx[N], &mem->Q_tilde, 0, 0, mem->RSQrq[N], nu[N], nu[N]);
    blasfeo_dgead(nx[N], nx[N], -1.0, &mem->Q_tilde, 0, 0, &mem->Q_bar, 0, 0);
    blasfeo_dtrtr_l(nx[N], &mem->Q_bar, 0, 0, &mem->Q_bar, 0, 0);

    for (ii = N-1; ii >= 0; --ii)
    {
        blasfeo_drowin(nx[ii+1], 1.0, mem->b[ii], 0, mem->BAbt[ii], nu[ii]+nx[ii], 0);
        blasfeo_drowin(nu[ii]+nx[ii], 1.0, mem->rq[ii], 0, mem->RSQrq[ii], nu[ii]+nx[ii], 0);

        blasfeo_dgecp(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], mem->RSQrq[ii], 0, 0, &mem->original_RSQrq[ii], 0, 0);

        // printf("----------------\n");
        // printf("--- stage %d ---\n", i);
        // printf("----------------\n");

        // printf("QSR\n");
        // blasfeo_print_dmat(nx+nu+1, nx+nu, &work->qp_in->RSQrq[i], 0, 0);

        // printf("Q_bar\n");
        // blasfeo_print_dmat(nx, nx, &Q_bar, 0, 0);

        // printf("BAbt\n");
        // blasfeo_print_dmat(nx+nu, nx, &work->qp_in->BAbt[i], 0, 0);

        // TODO implement using cholesky
        blasfeo_dgemm_nt(nu[ii]+nx[ii], nx[ii], nx[ii+1], 1.0, mem->BAbt[ii], 0, 0, &mem->Q_bar, 0, 0, 0.0, &mem->BAQ, 0, 0, &mem->BAQ, 0, 0);
        blasfeo_dsyrk_ln_mn(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], nx[ii+1], 1.0, mem->BAbt[ii], 0, 0, &mem->BAQ, 0, 0, 1.0, mem->RSQrq[ii], 0, 0, mem->RSQrq[ii], 0, 0);

		blasfeo_drowex(nu[ii]+nx[ii], 1.0, mem->RSQrq[ii], nu[ii]+nx[ii], 0, mem->rq[ii], 0);

        // printf("BAQ\n");
        // blasfeo_print_dmat(nx+nu, nx, &BAQ, 0, 0);

        blasfeo_unpack_dmat(nu[ii], nu[ii], mem->RSQrq[ii], 0, 0, mem->R, nu[ii]);
        acados_eigen_decomposition(nu[ii], mem->R, mem->V, mem->d, mem->e);

        bool needs_regularization = false;
        for (jj = 0; jj < nu[ii]; jj++)
            if (mem->d[jj] < 1e-10)
                needs_regularization = true;

        if (needs_regularization)
        {
			blasfeo_dgecp(nu[ii]+nx[ii], nu[ii]+nx[ii], mem->RSQrq[ii], 0, 0, &mem->tmp_RSQ, 0, 0);
			// TODO project only nu instead ???????????
			// TODO compute correction as a separate matrix, and apply to original_RSQrq too (TODO change this name then)
//            acados_mirror(nu[ii]+nx[ii], mem->reg_hess, mem->V, mem->d, mem->e, 1e-4);
#if 1
			blasfeo_unpack_dmat(nu[ii]+nx[ii], nu[ii]+nx[ii], mem->RSQrq[ii], 0, 0, mem->reg_hess, nu[ii]+nx[ii]);
			acados_project(nu[ii]+nx[ii], mem->reg_hess, mem->V, mem->d, mem->e, 1e-4);
			blasfeo_pack_dmat(nu[ii]+nx[ii], nu[ii]+nx[ii], mem->reg_hess, nu[ii]+nx[ii], mem->RSQrq[ii], 0, 0);
#else
			blasfeo_unpack_dmat(nu[ii], nu[ii], mem->RSQrq[ii], 0, 0, mem->reg_hess, nu[ii]);
			acados_project(nu[ii], mem->reg_hess, mem->V, mem->d, mem->e, 1e-4);
			blasfeo_pack_dmat(nu[ii], nu[ii], mem->reg_hess, nu[ii], mem->RSQrq[ii], 0, 0);
#endif
//            blasfeo_unpack_dmat(nu[ii], nu[ii], mem->RSQrq[ii], 0, 0, mem->reg_hess, nu[ii]);
//            acados_mirror(nu[ii], mem->reg_hess, mem->V, mem->d, mem->e, 1e-4);
//            acados_project(nu[ii], mem->reg_hess, mem->V, mem->d, mem->e, 1e-4);
//            blasfeo_pack_dmat(nu[ii], nu[ii], mem->reg_hess, nu[ii], mem->RSQrq[ii], 0, 0);
			blasfeo_dgead(nu[ii]+nx[ii], nu[ii]+nx[ii], -1.0, mem->RSQrq[ii], 0, 0, &mem->tmp_RSQ, 0, 0);
			blasfeo_dgead(nu[ii]+nx[ii], nu[ii]+nx[ii], -1.0, &mem->tmp_RSQ, 0, 0, &mem->original_RSQrq[ii], 0, 0);
        }

        // printf("QSR_hat\n");
        // blasfeo_print_dmat(nx+nu+1, nx+nu, &work->qp_in->RSQrq[i], 0, 0);


		// backup Q
        blasfeo_dgecp(nx[ii], nx[ii], mem->RSQrq[ii], nu[ii], nu[ii], &mem->Q_bar, 0, 0);

        // R = L * L^T
        blasfeo_dpotrf_l(nu[ii], mem->RSQrq[ii], 0, 0, &mem->L, 0, 0);
        // Q = S^T * L^-T
        blasfeo_dgecp(nx[ii], nu[ii], mem->RSQrq[ii], nu[ii], 0, &mem->St_copy, 0, 0);
//        blasfeo_dtrsm_rltn(nx[ii], nu[ii], 1.0, &mem->L, 0, 0, &mem->St_copy, 0, 0, &mem->Q_tilde, 0, 0);
        blasfeo_dtrsm_rltn(nx[ii], nu[ii], 1.0, &mem->L, 0, 0, &mem->St_copy, 0, 0, &mem->St_copy, 0, 0);

        // Q = S^T * R^-1 * S + delta*I
		blasfeo_dgese(nx[ii], nx[ii], 0.0, &mem->delta_eye, 0, 0);
		blasfeo_ddiare(nx[ii], delta, &mem->delta_eye, 0, 0);
//        blasfeo_dsyrk_ln(nx[ii], nx[ii], 1.0, &mem->Q_tilde, 0, 0, &mem->Q_tilde, 0, 0, 1.0, &mem->delta_eye, 0, 0, mem->RSQrq[ii], nu[ii], nu[ii]);
        blasfeo_dsyrk_ln(nx[ii], nu[ii], 1.0, &mem->St_copy, 0, 0, &mem->St_copy, 0, 0, 1.0, &mem->delta_eye, 0, 0, mem->RSQrq[ii], nu[ii], nu[ii]);

        // printf("H_tilde\n");
        // blasfeo_print_dmat(nu+nx, nu+nx, &work->qp_in->RSQrq[i], 0, 0);

        // make symmetric
//        blasfeo_dtrtr_l(nx[ii], &mem->Q_bar, 0, 0, &mem->Q_bar, 0, 0);

        // TODO take from b !!!!!!
//        for (jj = 0; jj < nx[ii+1]; jj++)
//            BLASFEO_DVECEL(&mem->b2, jj) = BLASFEO_DMATEL(mem->BAbt[ii], nu[ii]+nx[ii], jj);

        // TODO nx stage is not consistent with above !!!!!!!
//        blasfeo_dgemv_n(nx[ii+1], nx[ii+1], 1.0, &mem->Q_bar, 0, 0, &mem->b2, 0, 0.0, &mem->grad, 0, &mem->grad, 0);
//        blasfeo_dgemv_n(nu[ii]+nx[ii], nx[ii+1], 1.0, mem->BAbt[ii], 0, 0, &mem->grad, 0, 0.0, &mem->b2, 0, &mem->b2, 0);

//        for (jj = 0; jj < nu[ii]+nx[ii]; jj++)
            // TODO maybe 'b' is a bad naming...
//            BLASFEO_DMATEL(mem->RSQrq[ii], nu[ii]+nx[ii], jj) = BLASFEO_DMATEL(mem->RSQrq[ii], nu[ii]+nx[ii], jj) + BLASFEO_DVECEL(&mem->b2, jj);

        blasfeo_dgead(nx[ii], nx[ii], -1.0, mem->RSQrq[ii], nu[ii], nu[ii], &mem->Q_bar, 0, 0);

        // make symmetric
        blasfeo_dtrtr_l(nx[ii], &mem->Q_bar, 0, 0, &mem->Q_bar, 0, 0);

    }

    return;
}



// NOTE this only considers the case of (dynamcs) equality constraints (no inequality constraints)
// TODO inequality constraints case
void ocp_nlp_reg_convexify_correct_dual_sol(void *config, ocp_nlp_reg_dims *dims, void *opts_, void *mem_)
{
    ocp_nlp_reg_convexify_memory *mem = mem_;
    // ocp_nlp_reg_convexify_opts *opts = opts_;

    int ii;

    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *nbx = dims->nbx;
    int *nbu = dims->nbu;
    int *ng = dims->ng;

//	printf("\nRSQrq reg\n");
//	for(ii=0; ii<=N; ii++)
//		blasfeo_print_dmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], mem->RSQrq[ii], 0, 0);

	// restore original hessian and gradient
	for(ii=0; ii<=N; ii++)
	{
        blasfeo_dgecp(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], &mem->original_RSQrq[ii], 0, 0, mem->RSQrq[ii], 0, 0);
		blasfeo_drowex(nu[ii]+nx[ii], 1.0, &mem->original_RSQrq[ii], nu[ii]+nx[ii], 0, mem->rq[ii], 0);
	}

//	printf("\nRSQrq orig\n");
//	for(ii=0; ii<=N; ii++)
//		blasfeo_print_dmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], mem->RSQrq[ii], 0, 0);



	// compute multipliers of equality constraints

#if 0
    blasfeo_dgemv_n(nx[N], nu[N]+nx[N], 1.0, mem->RSQrq[N], nu[N], 0, mem->ux[N], 0, 1.0, mem->rq[N], nu[N], mem->pi[N-1], 0);

    for(ii=1; ii<N; ii++)
    {
		ss = N-ii;
        blasfeo_dgemv_n(nx[ss], nu[ss]+nx[ss], 1.0, mem->RSQrq[ss], nu[ss], 0, mem->ux[ss], 0, 1.0, mem->rq[ss], nu[ss], mem->pi[ss-1], 0);
        blasfeo_dgemv_n(nx[ss], nx[ss+1], 1.0, mem->BAbt[ss], nu[ss], 0, mem->pi[ss], 0, 1.0, mem->pi[ss-1], 0, mem->pi[ss-1], 0);
    }

#else

	// last stage
	blasfeo_dveccp(nx[N], mem->rq[N], nu[N], &mem->tmp_nuxM, nu[N]);
	blasfeo_daxpy(nbu[N]+nbx[N]+ng[N], -1.0, mem->lam[N], 0, mem->lam[N], nbu[N]+nbx[N]+ng[N], &mem->tmp_nbgM, 0);
	blasfeo_dvecad_sp(nbu[N]+nbx[N], 1.0, &mem->tmp_nbgM, 0, mem->idxb[N], &mem->tmp_nuxM, 0);
	// TODO avoid to multiply by R ???
	blasfeo_dsymv_l(nu[N]+nx[N], 1.0, mem->RSQrq[N], 0, 0, mem->ux[N], 0, 1.0, &mem->tmp_nuxM, 0, &mem->tmp_nuxM, 0);
	blasfeo_dgemv_n(nx[N], ng[N], 1.0, mem->DCt[N], nu[N], 0, &mem->tmp_nbgM, nbu[N]+nbx[N], 1.0, &mem->tmp_nuxM, nu[N], &mem->tmp_nuxM, nu[N]);
	blasfeo_dveccp(nx[N], &mem->tmp_nuxM, nu[N], mem->pi[N-1], 0);

	// middle stages
	for(ii=0; ii<N-1; ii++)
		{
		blasfeo_dveccp(nx[N-1-ii], mem->rq[N-1-ii], nu[N-1-ii], &mem->tmp_nuxM, nu[N-1-ii]);
		blasfeo_daxpy(nbu[N-1-ii]+nbx[N-1-ii]+ng[N-1-ii], -1.0, mem->lam[N-1-ii], 0, mem->lam[N-1-ii], nbu[N-1-ii]+nbx[N-1-ii]+ng[N-1-ii], &mem->tmp_nbgM, 0);
		blasfeo_dvecad_sp(nbu[N-1-ii]+nbx[N-1-ii], 1.0, &mem->tmp_nbgM, 0, mem->idxb[N-1-ii], &mem->tmp_nuxM, 0);
		// TODO avoid to multiply by R ???
		blasfeo_dsymv_l(nu[N-1-ii]+nx[N-1-ii], 1.0, mem->RSQrq[N-1-ii], 0, 0, mem->ux[N-1-ii], 0, 1.0, &mem->tmp_nuxM, 0, &mem->tmp_nuxM, 0);
		blasfeo_dgemv_n(nx[N-1-ii], nx[N-ii], 1.0, mem->BAbt[N-1-ii], nu[N-1-ii], 0, mem->pi[N-1-ii], 0, 1.0, &mem->tmp_nuxM, nu[N-1-ii], &mem->tmp_nuxM, nu[N-1-ii]);
		blasfeo_dgemv_n(nx[N-1-ii], ng[N-1-ii], 1.0, mem->DCt[N-1-ii], nu[N-1-ii], 0, &mem->tmp_nbgM, nbu[N-1-ii]+nbx[N-1-ii], 1.0, &mem->tmp_nuxM, nu[N-1-ii], &mem->tmp_nuxM, nu[N-1-ii]);
		blasfeo_dveccp(nx[N-1-ii], &mem->tmp_nuxM, nu[N-1-ii], mem->pi[N-2-ii], 0);
		}
#endif

    return;
}



void ocp_nlp_reg_convexify_config_initialize_default(ocp_nlp_reg_config *config)
{
    // dims
    config->dims_calculate_size = &ocp_nlp_reg_dims_calculate_size;
    config->dims_assign = &ocp_nlp_reg_dims_assign;
    config->dims_set = &ocp_nlp_reg_dims_set;
    // opts
    config->opts_calculate_size = &ocp_nlp_reg_convexify_opts_calculate_size;
    config->opts_assign = &ocp_nlp_reg_convexify_opts_assign;
    config->opts_initialize_default = &ocp_nlp_reg_convexify_opts_initialize_default;
    config->opts_set = &ocp_nlp_reg_convexify_opts_set;
    // memory
    config->memory_calculate_size = &ocp_nlp_reg_convexify_calculate_memory_size;
    config->memory_assign = &ocp_nlp_reg_convexify_assign_memory;
    config->memory_set = &ocp_nlp_reg_convexify_memory_set;
    config->memory_set_RSQrq_ptr = &ocp_nlp_reg_convexify_memory_set_RSQrq_ptr;
    config->memory_set_rq_ptr = &ocp_nlp_reg_convexify_memory_set_rq_ptr;
    config->memory_set_BAbt_ptr = &ocp_nlp_reg_convexify_memory_set_BAbt_ptr;
    config->memory_set_b_ptr = &ocp_nlp_reg_convexify_memory_set_b_ptr;
    config->memory_set_idxb_ptr = &ocp_nlp_reg_convexify_memory_set_idxb_ptr;
    config->memory_set_DCt_ptr = &ocp_nlp_reg_convexify_memory_set_DCt_ptr;
    config->memory_set_ux_ptr = &ocp_nlp_reg_convexify_memory_set_ux_ptr;
    config->memory_set_pi_ptr = &ocp_nlp_reg_convexify_memory_set_pi_ptr;
    config->memory_set_lam_ptr = &ocp_nlp_reg_convexify_memory_set_lam_ptr;
    // functions
    config->regularize_hessian = &ocp_nlp_reg_convexify_regularize_hessian;
    config->correct_dual_sol = &ocp_nlp_reg_convexify_correct_dual_sol;
}
