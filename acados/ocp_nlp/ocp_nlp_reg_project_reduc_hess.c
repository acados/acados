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


#include "acados/ocp_nlp/ocp_nlp_reg_project_reduc_hess.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "acados/ocp_nlp/ocp_nlp_reg_common.h"
#include "acados/utils/math.h"
#include "acados/utils/mem.h"

#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_blas.h"



/************************************************
 * opts
 ************************************************/

acados_size_t ocp_nlp_reg_project_reduc_hess_opts_calculate_size(void)
{
    return sizeof(ocp_nlp_reg_project_reduc_hess_opts);
}



void *ocp_nlp_reg_project_reduc_hess_opts_assign(void *raw_memory)
{
    return raw_memory;
}



void ocp_nlp_reg_project_reduc_hess_opts_initialize_default(void *config_, ocp_nlp_reg_dims *dims, void *opts_)
{
    ocp_nlp_reg_project_reduc_hess_opts *opts = opts_;

    opts->thr_eig = 1e-12;
    opts->min_eig = 1e-4;
    opts->min_pivot = 1e-12;
	opts->pivoting = 1;

    return;
}



void ocp_nlp_reg_project_reduc_hess_opts_set(void *config_, ocp_nlp_reg_dims *dims, void *opts_, char *field, void* value)
{

    ocp_nlp_reg_project_reduc_hess_opts *opts = opts_;

    if (!strcmp(field, "thr_eig"))
    {
        double *d_ptr = value;
        opts->thr_eig = *d_ptr;
    }
    else if (!strcmp(field, "min_eig"))
    {
        double *d_ptr = value;
        opts->min_eig = *d_ptr;
    }
    else if (!strcmp(field, "min_pivot"))
    {
        double *d_ptr = value;
        opts->min_pivot = *d_ptr;
    }
    else if (!strcmp(field, "pivoting"))
    {
        int *i_ptr = value;
        opts->pivoting = *i_ptr;
    }
    else
    {
        printf("\nerror: field %s not available in ocp_nlp_reg_project_reduc_hess_opts_set\n", field);
        exit(1);
    }

    return;
}



/************************************************
 * memory
 ************************************************/

acados_size_t ocp_nlp_reg_project_reduc_hess_memory_calculate_size(void *config_, ocp_nlp_reg_dims *dims, void *opts_)
{
    int *nx = dims->nx;
    int *nu = dims->nu;
    int N = dims->N;

    int ii;

    int nuxM = nu[0]+nx[0];
	int nuM = nu[0];
	int nxM = nx[0];
    for(ii=1; ii<=N; ii++)
    {
        nuxM = nu[ii]+nx[ii]>nuxM ? nu[ii]+nx[ii] : nuxM;
        nuM = nu[ii]>nuM ? nu[ii] : nuM;
        nxM = nx[ii]>nxM ? nx[ii] : nxM;
    }

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_reg_project_reduc_hess_memory);

    size += nuxM*nuxM*sizeof(double);  // reg_hess
    size += nuxM*nuxM*sizeof(double);  // V
    size += 2*nuxM*sizeof(double);     // d e
    size += (N+1)*sizeof(struct blasfeo_dmat *); // RSQrq
    size += N*sizeof(struct blasfeo_dmat *); // BAbt

    size += 1 * 64;

    size += blasfeo_memsize_dmat(nuxM, nuxM);     // L
    size += blasfeo_memsize_dmat(nuxM, nuxM);     // L2
    size += blasfeo_memsize_dmat(nuxM, nuxM);     // L3
    size += blasfeo_memsize_dmat(nxM, nuM);     // Ls
    size += blasfeo_memsize_dmat(nxM, nxM);     // P
    size += blasfeo_memsize_dmat(nuxM, nxM);     // AL

    return size;
}



void *ocp_nlp_reg_project_reduc_hess_memory_assign(void *config_, ocp_nlp_reg_dims *dims, void *opts_, void *raw_memory)
{
    int *nx = dims->nx;
    int *nu = dims->nu;
    int N = dims->N;

    int ii;

    int nuxM = nu[0]+nx[0];
	int nuM = nu[0];
	int nxM = nx[0];
    for(ii=1; ii<=N; ii++)
    {
        nuxM = nu[ii]+nx[ii]>nuxM ? nu[ii]+nx[ii] : nuxM;
        nuM = nu[ii]>nuM ? nu[ii] : nuM;
        nxM = nx[ii]>nxM ? nx[ii] : nxM;
    }

    char *c_ptr = (char *) raw_memory;

    ocp_nlp_reg_project_reduc_hess_memory *mem = (ocp_nlp_reg_project_reduc_hess_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_reg_project_reduc_hess_memory);

    mem->reg_hess = (double *) c_ptr;
    c_ptr += nuxM*nuxM*sizeof(double);  // reg_hess

    mem->V = (double *) c_ptr;
    c_ptr += nuxM*nuxM*sizeof(double);  // V

    mem->d = (double *) c_ptr;
    c_ptr += nuxM*sizeof(double); // d

    mem->e = (double *) c_ptr;
    c_ptr += nuxM*sizeof(double); // e

    mem->RSQrq = (struct blasfeo_dmat **) c_ptr;
    c_ptr += (N+1)*sizeof(struct blasfeo_dmat *); // RSQrq

    mem->BAbt = (struct blasfeo_dmat **) c_ptr;
    c_ptr += N*sizeof(struct blasfeo_dmat *); // BAbt

    align_char_to(64, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nuxM, nuxM, &mem->L, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nuxM, nuxM, &mem->L2, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nuxM, nuxM, &mem->L3, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nxM, nuM, &mem->Ls, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nxM, nxM, &mem->P, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nuxM, nxM, &mem->AL, &c_ptr);

    assert((char *) mem + ocp_nlp_reg_project_reduc_hess_memory_calculate_size(config_, dims, opts_) >= c_ptr);

    return mem;
}



void ocp_nlp_reg_project_reduc_hess_memory_set_RSQrq_ptr(ocp_nlp_reg_dims *dims, struct blasfeo_dmat *RSQrq, void *memory_)
{
    ocp_nlp_reg_project_reduc_hess_memory *memory = memory_;

    int ii;

	int N = dims->N;

    for(ii=0; ii<=N; ii++)
    {
        memory->RSQrq[ii] = RSQrq+ii;
    }

    return;
}



void ocp_nlp_reg_project_reduc_hess_memory_set_rq_ptr(ocp_nlp_reg_dims *dims, struct blasfeo_dvec *rq, void *memory_)
{
    return;
}



void ocp_nlp_reg_project_reduc_hess_memory_set_BAbt_ptr(ocp_nlp_reg_dims *dims, struct blasfeo_dmat *BAbt, void *memory_)
{
    ocp_nlp_reg_project_reduc_hess_memory *memory = memory_;

    int ii;

	int N = dims->N;

    for(ii=0; ii<N; ii++)
    {
        memory->BAbt[ii] = BAbt+ii;
    }

    return;
}



void ocp_nlp_reg_project_reduc_hess_memory_set_b_ptr(ocp_nlp_reg_dims *dims, struct blasfeo_dvec *b, void *memory_)
{
    return;
}



void ocp_nlp_reg_project_reduc_hess_memory_set_idxb_ptr(ocp_nlp_reg_dims *dims, int **idxb, void *memory_)
{
    return;
}



void ocp_nlp_reg_project_reduc_hess_memory_set_DCt_ptr(ocp_nlp_reg_dims *dims, struct blasfeo_dmat *DCt, void *memory_)
{
    return;
}



void ocp_nlp_reg_project_reduc_hess_memory_set_ux_ptr(ocp_nlp_reg_dims *dims, struct blasfeo_dvec *ux, void *memory_)
{
    return;
}



void ocp_nlp_reg_project_reduc_hess_memory_set_pi_ptr(ocp_nlp_reg_dims *dims, struct blasfeo_dvec *pi, void *memory_)
{
    return;
}



void ocp_nlp_reg_project_reduc_hess_memory_set_lam_ptr(ocp_nlp_reg_dims *dims, struct blasfeo_dvec *pi, void *memory_)
{
    return;
}



void ocp_nlp_reg_project_reduc_hess_memory_set(void *config_, ocp_nlp_reg_dims *dims, void *memory_, char *field, void *value)
{

    if(!strcmp(field, "RSQrq_ptr"))
    {
        struct blasfeo_dmat *RSQrq = value;
        ocp_nlp_reg_project_reduc_hess_memory_set_RSQrq_ptr(dims, RSQrq, memory_);
    }
    else if(!strcmp(field, "BAbt_ptr"))
    {
        struct blasfeo_dmat *BAbt = value;
        ocp_nlp_reg_project_reduc_hess_memory_set_BAbt_ptr(dims, BAbt, memory_);
    }
    else
    {
        printf("\nerror: field %s not available in ocp_nlp_reg_project_reduc_hess_set\n", field);
        exit(1);
    }

    return;
}



/************************************************
 * functions
 ************************************************/

void ocp_nlp_reg_project_reduc_hess_regularize_hessian(void *config, ocp_nlp_reg_dims *dims, void *opts_, void *mem_)
{
    ocp_nlp_reg_project_reduc_hess_memory *mem = (ocp_nlp_reg_project_reduc_hess_memory *) mem_;
    ocp_nlp_reg_project_reduc_hess_opts *opts = opts_;

    int ii, jj, kk, ll, ss, idx;

    int *nx = dims->nx;
    int *nu = dims->nu;
    int N = dims->N;

	struct blasfeo_dmat *L = &mem->L;
	struct blasfeo_dmat *L2 = &mem->L2;
	struct blasfeo_dmat *L3 = &mem->L3;
	struct blasfeo_dmat *Ls = &mem->Ls;
	struct blasfeo_dmat *P = &mem->P;
	struct blasfeo_dmat *AL = &mem->AL;

	int do_reg = 0;
	double pivot, tmp_el;


	// last + middle stages
	for(ii=-1; ii<N-1; ii++)
	{

		ss = N-ii-1;
		// last stage
		if(ss==N)
		{
			blasfeo_dtrcp_l(nu[ss]+nx[ss], mem->RSQrq[ss], 0, 0, L, 0, 0);
		}
		// middle stages
		else
		{
			blasfeo_dgemm_nt(nu[ss]+nx[ss], nx[ss+1], nx[ss+1], 1.0, mem->BAbt[ss], 0, 0, P, 0, 0, 0.0, AL, 0, 0, AL, 0, 0); // TODO symm
			blasfeo_dsyrk_ln(nu[ss]+nx[ss], nx[ss+1], 1.0, AL, 0, 0, mem->BAbt[ss], 0, 0, 1.0, mem->RSQrq[ss], 0, 0, L, 0, 0);
		}
		blasfeo_dtrtr_l(nu[ss]+nx[ss], L, 0, 0, L, 0, 0); // necessary ???

		// backup L in L3
		blasfeo_dgese(nu[ss]+nx[ss], nu[ss]+nx[ss], 0.0, L3, 0, 0);
		blasfeo_dgecp(nu[ss]+nx[ss], nu[ss], L, 0, 0, L3, 0, 0);

		// project L_R
		blasfeo_unpack_dmat(nu[ss], nu[ss], L, 0, 0, mem->reg_hess, nu[ss]);
		acados_eigen_decomposition(nu[ss], mem->reg_hess, mem->V, mem->d, mem->e);
		do_reg = 0;
		for(jj=0; jj<nu[ss]; jj++)
		{
			if(mem->d[jj]<opts->thr_eig)
			{
				mem->e[jj] = opts->min_eig - mem->d[jj];
				do_reg = 1;
			}
			else
			{
				mem->e[jj] = 0.0;
			}
		}
		if(do_reg)
		{
			acados_reconstruct_A(nu[ss], mem->reg_hess, mem->V, mem->e);
			blasfeo_dgese(nu[ss]+nx[ss], nu[ss]+nx[ss], 0.0, L2, 0, 0);
			blasfeo_pack_dmat(nu[ss], nu[ss], mem->reg_hess, nu[ss], L2, 0, 0);

			// apply reg to R
			blasfeo_dgead(nu[ss], nu[ss], 1.0, L2, 0, 0, mem->RSQrq[ss], 0, 0);
			// apply reg to L
			blasfeo_dgead(nu[ss], nu[ss], 1.0, L2, 0, 0, L, 0, 0);
		}

		// compute reg_schur
		blasfeo_dgecp(nu[ss]+nx[ss], nu[ss], L, 0, 0, L2, 0, 0);
		blasfeo_dpotrf_l_mn(nu[ss]+nx[ss], nu[ss], L2, 0, 0, L2, 0, 0);
		blasfeo_dgecp(nx[ss], nu[ss], L2, nu[ss], 0, Ls, 0, 0);
		blasfeo_dsyrk_ln_mn(nx[ss], nx[ss], nu[ss], -1.0, Ls, 0, 0, Ls, 0, 0, 0.0, L2, nu[ss], nu[ss], L2, nu[ss], nu[ss]);

		// compute true_schur
		if(do_reg)
		{
			for(jj=0; jj<nu[ss]; jj++)
			{
				if(opts->pivoting)
				{
					// find pivot element
					pivot = BLASFEO_DMATEL(L3, jj, jj);
					idx = jj;
					for(kk=jj+1; kk<nu[ss]; kk++)
					{
						tmp_el = BLASFEO_DMATEL(L3, kk, kk);
						if((tmp_el<pivot) & (tmp_el>-pivot))
						{
							pivot = BLASFEO_DMATEL(L3, kk, kk);
							idx = kk;
						}
					}
					// symmetric permute
					if(idx!=jj)
					{
						// top triangle
						for(kk=jj; kk<idx; kk++)
						{
							tmp_el = BLASFEO_DMATEL(L3, kk, jj);
							BLASFEO_DMATEL(L3, kk, jj) = BLASFEO_DMATEL(L3, idx-jj, idx-kk);
							BLASFEO_DMATEL(L3, idx-jj, idx-kk) = tmp_el;
						}
						// bottom rectangle
						for(kk=idx+1; kk<nu[ss]+nx[ss]; kk++)
                            {
							tmp_el = BLASFEO_DMATEL(L3, kk, jj);
							BLASFEO_DMATEL(L3, kk, jj) = BLASFEO_DMATEL(L3, kk, idx);
							BLASFEO_DMATEL(L3, kk, idx) = tmp_el;
						}
					}
				}
					
				pivot = BLASFEO_DMATEL(L3, jj, jj);
				if ((pivot<opts->min_pivot) & (pivot>-opts->min_pivot))
				{
					if(pivot<0.0)
						pivot = opts->min_pivot;
					else
						pivot = - opts->min_pivot;
				}
				pivot = 1.0/pivot;
				for(kk=jj+1; kk<nu[ss]+nx[ss]; kk++)
				{
					tmp_el = pivot * BLASFEO_DMATEL(L3, kk, jj);
					for(ll=kk; ll<nu[ss]+nx[ss]; ll++)
					{
						BLASFEO_DMATEL(L3, ll, kk) -= BLASFEO_DMATEL(L3, ll, jj) * tmp_el;
					}
				}
			}
		}

		// apply shur to P
		blasfeo_dgecp(nx[ss], nx[ss], L, nu[ss], nu[ss], P, 0, 0);
		if(do_reg)
//		if(0)
		{
			// P
			blasfeo_dgead(nx[ss], nx[ss], 1.0, L3, nu[ss], nu[ss], P, 0, 0);
			// Q
			blasfeo_dgead(nx[ss], nx[ss], -1.0, L2, nu[ss], nu[ss], mem->RSQrq[ss], nu[ss], nu[ss]);
			blasfeo_dgead(nx[ss], nx[ss],  1.0, L3, nu[ss], nu[ss], mem->RSQrq[ss], nu[ss], nu[ss]);
		}
		else
		{
			// P
			blasfeo_dgead(nx[ss], nx[ss], 1.0, L2, nu[ss], nu[ss], P, 0, 0);
		}
		blasfeo_dtrtr_l(nx[ss], P, 0, 0, P, 0, 0);

	}


	// first stage: factorize P in L too
	ss = 0;
	blasfeo_dgemm_nt(nu[ss]+nx[ss], nx[ss+1], nx[ss+1], 1.0, mem->BAbt[ss], 0, 0, P, 0, 0, 0.0, AL, 0, 0, AL, 0, 0); // TODO symm
	blasfeo_dsyrk_ln(nu[ss]+nx[ss], nx[ss+1], 1.0, AL, 0, 0, mem->BAbt[ss], 0, 0, 1.0, mem->RSQrq[ss], 0, 0, L, 0, 0);
	blasfeo_dtrtr_l(nu[ss]+nx[ss], L, 0, 0, L, 0, 0); // necessary ???
	blasfeo_unpack_dmat(nu[ss]+nx[ss], nu[ss]+nx[ss], L, 0, 0, mem->reg_hess, nu[ss]+nx[ss]);
	acados_eigen_decomposition(nu[ss]+nx[ss], mem->reg_hess, mem->V, mem->d, mem->e);
	for(jj=0; jj<nu[ss]+nx[ss]; jj++)
	{
		if(mem->d[jj]<opts->thr_eig)
			mem->e[jj] = opts->min_eig - mem->d[jj];
		else
			mem->e[jj] = 0.0;
	}
	acados_reconstruct_A(nu[ss]+nx[ss], mem->reg_hess, mem->V, mem->e);
	blasfeo_pack_dmat(nu[ss]+nx[ss], nu[ss]+nx[ss], mem->reg_hess, nu[ss]+nx[ss], L2, 0, 0);
	blasfeo_dgead(nu[ss]+nx[ss], nu[ss]+nx[ss], 1.0, L2, 0, 0, mem->RSQrq[ss], 0, 0);


//	printf("\nhessian after\n");
//	for(ii=0; ii<=N; ii++)
//	{
//		printf("\nii = %d\n", ii);
//		blasfeo_print_dmat(nu[ii]+nx[ii], nu[ii]+nx[ii], mem->RSQrq[ii], 0, 0);
//		blasfeo_unpack_dmat(nu[ii]+nx[ii], nu[ii]+nx[ii], mem->RSQrq[ii], 0, 0, mem->reg_hess, nu[ii]+nx[ii]);
//		acados_eigen_decomposition(nu[ii]+nx[ii], mem->reg_hess, mem->V, mem->d, mem->e);
//		d_print_mat(1, nu[ii]+nx[ii], mem->d, 1);
//	}
//	exit(1);

	return;
}



void ocp_nlp_reg_project_reduc_hess_correct_dual_sol(void *config, ocp_nlp_reg_dims *dims, void *opts_, void *mem_)
{
    return;
}



void ocp_nlp_reg_project_reduc_hess_config_initialize_default(ocp_nlp_reg_config *config)
{
    // dims
    config->dims_calculate_size = &ocp_nlp_reg_dims_calculate_size;
    config->dims_assign = &ocp_nlp_reg_dims_assign;
    config->dims_set = &ocp_nlp_reg_dims_set;
    // opts
    config->opts_calculate_size = &ocp_nlp_reg_project_reduc_hess_opts_calculate_size;
    config->opts_assign = &ocp_nlp_reg_project_reduc_hess_opts_assign;
    config->opts_initialize_default = &ocp_nlp_reg_project_reduc_hess_opts_initialize_default;
    config->opts_set = &ocp_nlp_reg_project_reduc_hess_opts_set;
    // memory
    config->memory_calculate_size = &ocp_nlp_reg_project_reduc_hess_memory_calculate_size;
    config->memory_assign = &ocp_nlp_reg_project_reduc_hess_memory_assign;
    config->memory_set = &ocp_nlp_reg_project_reduc_hess_memory_set;
    config->memory_set_RSQrq_ptr = &ocp_nlp_reg_project_reduc_hess_memory_set_RSQrq_ptr;
    config->memory_set_rq_ptr = &ocp_nlp_reg_project_reduc_hess_memory_set_rq_ptr;
    config->memory_set_BAbt_ptr = &ocp_nlp_reg_project_reduc_hess_memory_set_BAbt_ptr;
    config->memory_set_b_ptr = &ocp_nlp_reg_project_reduc_hess_memory_set_b_ptr;
    config->memory_set_idxb_ptr = &ocp_nlp_reg_project_reduc_hess_memory_set_idxb_ptr;
    config->memory_set_DCt_ptr = &ocp_nlp_reg_project_reduc_hess_memory_set_DCt_ptr;
    config->memory_set_ux_ptr = &ocp_nlp_reg_project_reduc_hess_memory_set_ux_ptr;
    config->memory_set_pi_ptr = &ocp_nlp_reg_project_reduc_hess_memory_set_pi_ptr;
    config->memory_set_lam_ptr = &ocp_nlp_reg_project_reduc_hess_memory_set_lam_ptr;
    // functions
    config->regularize_hessian = &ocp_nlp_reg_project_reduc_hess_regularize_hessian;
    config->correct_dual_sol = &ocp_nlp_reg_project_reduc_hess_correct_dual_sol;
}

