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

int ocp_nlp_reg_project_reduc_hess_opts_calculate_size(void)
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

    opts->epsilon = 1e-4;

    return;
}



void ocp_nlp_reg_project_reduc_hess_opts_set(void *config_, ocp_nlp_reg_dims *dims, void *opts_, char *field, void* value)
{

    ocp_nlp_reg_project_reduc_hess_opts *opts = opts_;

    if (!strcmp(field, "epsilon"))
    {
        double *d_ptr = value;
        opts->epsilon = *d_ptr;
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

int ocp_nlp_reg_project_reduc_hess_memory_calculate_size(void *config_, ocp_nlp_reg_dims *dims, void *opts_)
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

    int size = 0;

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
    int *nx = dims->nx;
    int *nu = dims->nu;

    for(ii=0; ii<=N; ii++)
    {
        memory->RSQrq[ii] = RSQrq+ii;
//        blasfeo_print_dmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], memory->RSQrq[ii], 0, 0);
    }

    return;
}



void ocp_nlp_reg_project_reduc_hess_memory_set_rq_ptr(ocp_nlp_reg_dims *dims, struct blasfeo_dvec *rq, void *memory_)
{
#if 0
    ocp_nlp_reg_project_reduc_hess_memory *memory = memory_;

    int ii;

    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;

    for(ii=0; ii<=N; ii++)
    {
        memory->rq[ii] = rq+ii;
//        blasfeo_print_dvec(nu[ii]+nx[ii], memory->rq[ii], 0);
    }
#endif

    return;
}



void ocp_nlp_reg_project_reduc_hess_memory_set_BAbt_ptr(ocp_nlp_reg_dims *dims, struct blasfeo_dmat *BAbt, void *memory_)
{
    ocp_nlp_reg_project_reduc_hess_memory *memory = memory_;

    int ii;

    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;

    for(ii=0; ii<N; ii++)
    {
        memory->BAbt[ii] = BAbt+ii;
//        blasfeo_print_dmat(nu[ii]+nx[ii]+1, nx[ii+1], memory->BAbt[ii], 0, 0);
    }

    return;
}



void ocp_nlp_reg_project_reduc_hess_memory_set_b_ptr(ocp_nlp_reg_dims *dims, struct blasfeo_dvec *b, void *memory_)
{
#if 0
    ocp_nlp_reg_project_reduc_hess_memory *memory = memory_;

    int ii;

    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;

    for(ii=0; ii<N; ii++)
    {
        memory->b[ii] = b+ii;
//        blasfeo_print_dvec(nx[ii=1], memory->b[ii], 0);
    }
#endif

    return;
}



void ocp_nlp_reg_project_reduc_hess_memory_set_ux_ptr(ocp_nlp_reg_dims *dims, struct blasfeo_dvec *ux, void *memory_)
{
#if 0
    ocp_nlp_reg_project_reduc_hess_memory *memory = memory_;

    int ii;

    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;

    for(ii=0; ii<=N; ii++)
    {
        memory->ux[ii] = ux+ii;
//        blasfeo_print_dvec(nu[ii]+nx[ii], memory->ux[ii], 0);
    }
#endif

    return;
}



void ocp_nlp_reg_project_reduc_hess_memory_set_pi_ptr(ocp_nlp_reg_dims *dims, struct blasfeo_dvec *pi, void *memory_)
{
#if 0
    ocp_nlp_reg_project_reduc_hess_memory *memory = memory_;

    int ii;

    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;

    for(ii=0; ii<N; ii++)
    {
        memory->pi[ii] = pi+ii;
//        blasfeo_print_dvec(nx[ii+1], memory->pi[ii], 0);
    }
#endif

    return;
}



void ocp_nlp_reg_project_reduc_hess_memory_set(void *config_, ocp_nlp_reg_dims *dims, void *memory_, char *field, void *value)
{

    if(!strcmp(field, "RSQrq_ptr"))
    {
        struct blasfeo_dmat *RSQrq = value;
        ocp_nlp_reg_project_reduc_hess_memory_set_RSQrq_ptr(dims, RSQrq, memory_);
    }
//    else if(!strcmp(field, "rq_ptr"))
//    {
//        struct blasfeo_dvec *rq = value;
//        ocp_nlp_reg_project_reduc_hess_memory_set_rq_ptr(dims, rq, memory_);
//    }
    else if(!strcmp(field, "BAbt_ptr"))
    {
        struct blasfeo_dmat *BAbt = value;
        ocp_nlp_reg_project_reduc_hess_memory_set_BAbt_ptr(dims, BAbt, memory_);
    }
//    else if(!strcmp(field, "b_ptr"))
//    {
//        struct blasfeo_dvec *b = value;
//        ocp_nlp_reg_project_reduc_hess_memory_set_b_ptr(dims, b, memory_);
//    }
//    else if(!strcmp(field, "ux_ptr"))
//    {
//        struct blasfeo_dvec *ux = value;
//        ocp_nlp_reg_project_reduc_hess_memory_set_ux_ptr(dims, ux, memory_);
//    }
//    else if(!strcmp(field, "pi_ptr"))
//    {
//        struct blasfeo_dvec *pi = value;
//        ocp_nlp_reg_project_reduc_hess_memory_set_pi_ptr(dims, pi, memory_);
//    }
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

    int ii, jj, kk, ll, ss;

//printf("\nhola\n");
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

//	for(ii=0; ii<=N; ii++)
//		blasfeo_print_dmat(nu[ii]+nx[ii], nu[ii]+nx[ii], mem->RSQrq[ii], 0, 0);
//	exit(1);

//printf("\nin project\n");
	// last stage
	ss = N;
	blasfeo_dtrtr_l(nu[ss]+nx[ss], mem->RSQrq[ss], 0, 0, mem->RSQrq[ss], 0, 0); // necessary ???
	// TODO !!!!!!!!!!!!!!!!!!
//	blasfeo_unpack_dmat(nu[ss], nu[ss], mem->RSQrq[ss], 0, 0, mem->reg_hess, nu[ss]);
//	acados_project(nu[ss], mem->reg_hess, mem->V, mem->d, mem->e, opts->epsilon);
//	blasfeo_pack_dmat(nu[ss], nu[ss], mem->reg_hess, nu[ss], mem->RSQrq[ss], 0, 0);
	blasfeo_dpotrf_l_mn(nu[ss]+nx[ss], nu[ss], mem->RSQrq[ss], 0, 0, L, 0, 0);
//printf("\nii = %d\n", ss);
//blasfeo_print_dmat(nu[ss]+nx[ss], nu[ss]+nx[ss], L, 0, 0);
	blasfeo_dgecp(nx[ss], nu[ss], L, nu[ss], 0, Ls, 0, 0);
	blasfeo_dsyrk_ln_mn(nx[ss], nx[ss], nu[ss], -1.0, Ls, 0, 0, Ls, 0, 0, 1.0, mem->RSQrq[ss], nu[ss], nu[ss], P, 0, 0);
	blasfeo_dtrtr_l(nx[ss], P, 0, 0, P, 0, 0);

	// middle stages
	for(ii=0; ii<N-1; ii++)
	{
		ss = N-ii-1;
//printf("\nss = %d\n", ss);
		blasfeo_dgemm_nt(nu[ss]+nx[ss], nx[ss+1], nx[ss+1], 1.0, mem->BAbt[ss], 0, 0, P, 0, 0, 0.0, AL, 0, 0, AL, 0, 0); // TODO symm
		blasfeo_dsyrk_ln(nu[ss]+nx[ss], nx[ss+1], 1.0, AL, 0, 0, mem->BAbt[ss], 0, 0, 1.0, mem->RSQrq[ss], 0, 0, L, 0, 0);
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
			if(mem->d[jj]<opts->epsilon)
			{
				mem->e[jj] = opts->epsilon - mem->d[jj];
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
//printf("\nL2\n");
//blasfeo_print_dmat(nu[ss]+nx[ss], nu[ss]+nx[ss], L2, 0, 0);

		// compute true_schur
		if(do_reg)
		{
			for(jj=0; jj<nu[ss]; jj++)
			{
				// TODO check for singularity of pivot !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				pivot = BLASFEO_DMATEL(L3, jj, jj);
//				printf("\n%f\n", pivot);
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
//printf("\nL3\n");
//blasfeo_print_dmat(nu[ss]+nx[ss], nu[ss]+nx[ss], L3, 0, 0);
		}

		// apply shur to P
		blasfeo_dgecp(nx[ss], nx[ss], L, nu[ss], nu[ss], P, 0, 0);
		if(do_reg)
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
//printf("\nP\n");
//blasfeo_print_dmat(nx[ss], nx[ss], P, 0, 0);

	}

	// first stage: factorize P in L too
	ss = 0;
//printf("\nss %d\n", ss);
	blasfeo_dgemm_nt(nu[ss]+nx[ss], nx[ss+1], nx[ss+1], 1.0, mem->BAbt[ss], 0, 0, P, 0, 0, 0.0, AL, 0, 0, AL, 0, 0); // TODO symm
	blasfeo_dsyrk_ln(nu[ss]+nx[ss], nx[ss+1], 1.0, AL, 0, 0, mem->BAbt[ss], 0, 0, 1.0, mem->RSQrq[ss], 0, 0, L, 0, 0);
	blasfeo_dtrtr_l(nu[ss]+nx[ss], L, 0, 0, L, 0, 0); // necessary ???
	// TODO regularize RSQ (also SQ !!!)
//blasfeo_print_dmat(nu[ss]+nx[ss], nu[ss]+nx[ss], L, 0, 0);
	blasfeo_unpack_dmat(nu[ss]+nx[ss], nu[ss]+nx[ss], L, 0, 0, mem->reg_hess, nu[ss]+nx[ss]);
	acados_eigen_decomposition(nu[ss]+nx[ss], mem->reg_hess, mem->V, mem->d, mem->e);
	for(jj=0; jj<nu[ss]+nx[ss]; jj++)
	{
		if(mem->d[jj]<opts->epsilon)
			mem->e[jj] = opts->epsilon - mem->d[jj];
		else
			mem->e[jj] = 0.0;
	}
//d_print_mat(nu[ss]+nx[ss], nu[ss]+nx[ss], mem->V, nu[ss]+nx[ss]);
//d_print_mat(1, nu[ss]+nx[ss], mem->d, 1);
//d_print_mat(1, nu[ss]+nx[ss], mem->e, 1);
	acados_reconstruct_A(nu[ss]+nx[ss], mem->reg_hess, mem->V, mem->e);
//d_print_mat(nu[ss]+nx[ss], nu[ss]+nx[ss], mem->reg_hess, nu[ss]+nx[ss]);
	blasfeo_pack_dmat(nu[ss]+nx[ss], nu[ss]+nx[ss], mem->reg_hess, nu[ss]+nx[ss], L2, 0, 0);
	blasfeo_dgead(nu[ss]+nx[ss], nu[ss]+nx[ss], 1.0, L2, 0, 0, mem->RSQrq[ss], 0, 0);
//blasfeo_print_dmat(nu[ss]+nx[ss], nu[ss]+nx[ss], mem->RSQrq[ss], 0, 0);
	// TODO till here
//blasfeo_dgead(nu[ss]+nx[ss], nu[ss]+nx[ss], 1.0, L2, 0, 0, L, 0, 0);
//blasfeo_dpotrf_l(nu[ss]+nx[ss], L, 0, 0, L, 0, 0);
//printf("\nL0\n");
//blasfeo_print_dmat(nu[ss]+nx[ss], nu[ss]+nx[ss], L, 0, 0);

//	blasfeo_print_dmat(nx[ii], nx[ii], P, 0, 0);

//	exit(1);

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
    config->memory_set_ux_ptr = &ocp_nlp_reg_project_reduc_hess_memory_set_ux_ptr;
    config->memory_set_pi_ptr = &ocp_nlp_reg_project_reduc_hess_memory_set_pi_ptr;
    // functions
    config->regularize_hessian = &ocp_nlp_reg_project_reduc_hess_regularize_hessian;
    config->correct_dual_sol = &ocp_nlp_reg_project_reduc_hess_correct_dual_sol;
}

