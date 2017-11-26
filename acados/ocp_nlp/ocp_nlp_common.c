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

#include "acados/ocp_nlp/ocp_nlp_common.h"

#include <assert.h>

#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "acados/utils/mem.h"
#include "hpipm/include/hpipm_d_ocp_qp_dim.h"

#ifndef YT
// TODO(dimitris): TEMP!!! REMOVE AFTER tmp_allocate_ocp_nlp_in_sim_solver IS GONE!
#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"
#endif

int number_of_primal_vars(ocp_nlp_dims *dims)
{
    int num_vars = 0;
    for (int ii = 0; ii <= dims->N; ii++) {
        num_vars += dims->nx[ii] + dims->nu[ii];
    }
    return num_vars;
}



void cast_nlp_dims_to_qp_dims(ocp_qp_dims *qp_dims, ocp_nlp_dims *nlp_dims)
{
    qp_dims->N = nlp_dims->N;
    qp_dims->nx = nlp_dims->nx;
    qp_dims->nu = nlp_dims->nu;
    qp_dims->nb = nlp_dims->nb;
    qp_dims->nbx = nlp_dims->nbx;
    qp_dims->nbu = nlp_dims->nbu;
    qp_dims->ng = nlp_dims->ng;
    qp_dims->ns = nlp_dims->ns;

    // TODO(dimitris): probably redundant (can also remove hpipm header)
    qp_dims->memsize = d_memsize_ocp_qp_dim(qp_dims->N);
}



// TODO(dimitris): fix order of funs
int ocp_nlp_in_calculate_size(ocp_nlp_dims *dims)
{
    int N = dims->N;

    int size = sizeof(ocp_nlp_in);

    size += sizeof(ocp_nlp_dims);
    size += 8*sizeof(int)*(N+1);  // dims

    size += sizeof(int *)*(N+1);  // idxb
    size += 2*sizeof(double *)*(N+1);  // lb, ub

    size += 4*sizeof(double *)*(N+1);  // lg, ug, Cx, Cu

    size += 2*sizeof(double *)*(N+1);  // lh, uh

    // TODO(dimitris): check arguments for cost type
    size += sizeof(ocp_nlp_ls_cost);
    size += 2*sizeof(double *)*(N+1);  // W, yref

    for (int ii = 0; ii < N+1; ii++)
    {
        size += sizeof(int)*(dims->nbx[ii]+dims->nbu[ii]);  // idxb
        size += 2*sizeof(double)*(dims->nbx[ii]+dims->nbu[ii]);  // lb, ub

        size += 2*sizeof(double)*dims->ng[ii];  // lg, ug
        size += sizeof(double)*dims->nx[ii]*dims->ng[ii];  // Cx
        size += sizeof(double)*dims->nu[ii]*dims->ng[ii];  // Cu

        size += 2*sizeof(double)*dims->nh[ii];  // lh, uh

        size += sizeof(double)*(dims->nx[ii]+dims->nu[ii])*(dims->nx[ii]+dims->nu[ii]);  // W
        size += sizeof(double)*(dims->nx[ii]+dims->nu[ii]);  // yref
    }

    // TODO(dimitris): add alignment when strvecs/strmats are used
    return size;
}


#ifndef YT
// TEMP!!!
static void tmp_allocate_ocp_nlp_in_sim_solver(int_t N, int_t *nx, int_t *nu, int_t num_stages,
    ocp_nlp_in *const nlp)
{
    nlp->sim = (sim_solver *)calloc(N, sizeof(sim_solver));
    for (int_t i = 0; i < N; i++) {
        int_t nx_i = nx[i];
        int_t nu_i = nu[i];
        nlp->sim[i].in = (sim_in *)malloc(sizeof(sim_in));
        d_zeros(&nlp->sim[i].in->x, nx_i, 1);
        d_zeros(&nlp->sim[i].in->u, nu_i, 1);
        d_zeros(&nlp->sim[i].in->S_forw, nx_i, nx_i + nu_i);
        for (int_t j = 0; j < nx_i; j++)
        nlp->sim[i].in->S_forw[j * (nx_i + 1)] = 1.0;

        d_zeros(&nlp->sim[i].in->S_adj, nx_i + nu_i, 1);
        d_zeros(&nlp->sim[i].in->grad_K, nx_i * num_stages, 1);

        int_t nx_i1 = nx[i + 1];
        nlp->sim[i].out = (sim_out *)malloc(sizeof(sim_out));
        d_zeros(&nlp->sim[i].out->xn, nx_i1, 1);
        d_zeros(&nlp->sim[i].out->S_forw, nx_i1, nx_i + nu_i);
        d_zeros(&nlp->sim[i].out->S_adj, nx_i + nu_i, 1);
        d_zeros(&nlp->sim[i].out->grad, nx_i + nu_i, 1);
        int_t nhess = (nx_i + nu_i + 1) * (nx_i + nu_i) / 2;
        d_zeros(&nlp->sim[i].out->S_hess, nhess, 1);
        nlp->sim[i].out->info = (sim_info *)malloc(sizeof(sim_info));

        nlp->sim[i].mem = NULL;
    }
}


// TEMP!!!
void tmp_free_ocp_nlp_in_sim_solver(ocp_nlp_in *const nlp) {
    for (int_t i = 0; i < nlp->dims->N; i++) {
        free(nlp->sim[i].in->x);
        free(nlp->sim[i].in->u);
        free(nlp->sim[i].in->S_forw);
        free(nlp->sim[i].in->S_adj);
        free(nlp->sim[i].in->grad_K);
        free(nlp->sim[i].in);

        free(nlp->sim[i].out->xn);
        free(nlp->sim[i].out->S_forw);
        free(nlp->sim[i].out->S_adj);
        free(nlp->sim[i].out->info);
        free(nlp->sim[i].out->grad);
        free(nlp->sim[i].out->S_hess);
        free(nlp->sim[i].out);
    }
    free(nlp->sim);
}
#endif // YT


// TODO(dimitris): move num_stages inside args, as nested integrator args
ocp_nlp_in *assign_ocp_nlp_in(ocp_nlp_dims *dims, int num_stages, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    int N = dims->N;

    ocp_nlp_in *in = (ocp_nlp_in *)c_ptr;
    c_ptr += sizeof(ocp_nlp_in);

    in->dims = (ocp_nlp_dims *)c_ptr;
    c_ptr += sizeof(ocp_nlp_dims);

    // TODO(dimitris): check arguments for cost type
    ocp_nlp_ls_cost *cost = (ocp_nlp_ls_cost *)c_ptr;
    in->cost = (void *)cost;
    c_ptr += sizeof(ocp_nlp_ls_cost);

    // double pointers
    assign_int_ptrs(N+1, &in->idxb, &c_ptr);
    assign_double_ptrs(N+1, &in->lb, &c_ptr);
    assign_double_ptrs(N+1, &in->ub, &c_ptr);

    assign_double_ptrs(N+1, &in->lg, &c_ptr);
    assign_double_ptrs(N+1, &in->ug, &c_ptr);
    assign_double_ptrs(N+1, &in->Cx, &c_ptr);
    assign_double_ptrs(N+1, &in->Cu, &c_ptr);

    assign_double_ptrs(N+1, &in->lh, &c_ptr);
    assign_double_ptrs(N+1, &in->uh, &c_ptr);

    assign_double_ptrs(N+1, &cost->W, &c_ptr);
    assign_double_ptrs(N+1, &cost->y_ref, &c_ptr);

    // assign data

    // doubles
    for (int ii = 0; ii < N+1; ii++)
    {
        assign_double(dims->nbx[ii]+dims->nbu[ii], &in->lb[ii], &c_ptr);
        assign_double(dims->nbx[ii]+dims->nbu[ii], &in->ub[ii], &c_ptr);

        assign_double(dims->ng[ii], &in->lg[ii], &c_ptr);
        assign_double(dims->ng[ii], &in->ug[ii], &c_ptr);
        assign_double(dims->ng[ii]*dims->nx[ii], &in->Cx[ii], &c_ptr);
        assign_double(dims->ng[ii]*dims->nu[ii], &in->Cu[ii], &c_ptr);

        assign_double(dims->nh[ii], &in->lh[ii], &c_ptr);
        assign_double(dims->nh[ii], &in->uh[ii], &c_ptr);

        assign_double((dims->nx[ii]+dims->nu[ii])*(dims->nx[ii]+dims->nu[ii]), &cost->W[ii], &c_ptr);
        assign_double(dims->nx[ii]+dims->nu[ii], &cost->y_ref[ii], &c_ptr);
    }
    // ints
    for (int ii = 0; ii < N+1; ii++)
    {
        assign_int(dims->nbx[ii]+dims->nbu[ii], &in->idxb[ii], &c_ptr);
    }

    // assign and copy dimensions
    in->dims->N = N;

    assign_int(N+1, &in->dims->nx, &c_ptr);
    assign_int(N+1, &in->dims->nu, &c_ptr);
    assign_int(N+1, &in->dims->nb, &c_ptr);
    assign_int(N+1, &in->dims->nbx, &c_ptr);
    assign_int(N+1, &in->dims->nbu, &c_ptr);
    assign_int(N+1, &in->dims->ng, &c_ptr);
    assign_int(N+1, &in->dims->nh, &c_ptr);
    assign_int(N+1, &in->dims->ns, &c_ptr);

    for (int ii = 0; ii < N+1; ii++)
    {
        in->dims->nx[ii] = dims->nx[ii];
        in->dims->nu[ii] = dims->nu[ii];
        in->dims->nb[ii] = dims->nbx[ii] + dims->nbu[ii];  // dims->nb[ii];
        in->dims->nbx[ii] = dims->nbx[ii];
        in->dims->nbu[ii] = dims->nbu[ii];
        in->dims->ng[ii] = dims->ng[ii];
        in->dims->nh[ii] = dims->nh[ii];
        in->dims->ns[ii] = dims->ns[ii];
    }

    #ifndef YT
    tmp_allocate_ocp_nlp_in_sim_solver(N, dims->nx, dims->nu, num_stages, in);
    #endif

    // printf("diff = %lld\n", (long long int)(raw_memory + ocp_nlp_in_calculate_size(dims, args)) - (long long int)c_ptr);
    assert((char *) raw_memory + ocp_nlp_in_calculate_size(dims) == c_ptr);

    return in;
}



int ocp_nlp_out_calculate_size(ocp_nlp_dims *dims)
{
    int N = dims->N;

    int size = sizeof(ocp_nlp_out);

    size += sizeof(double *) * (N + 1);  // x
    size += sizeof(double *) * (N + 1);  // u
    size += sizeof(double *) * (N + 1);  // lam
    size += sizeof(double *) * N;  // pi

    for (int ii = 0; ii <= N; ii++)
    {
        size += sizeof(double)*dims->nx[ii];  // x
        size += sizeof(double)*dims->nu[ii];  // u
        size += sizeof(double)*2*(dims->nb[ii] + dims->ng[ii] + dims->nh[ii]);  // lam
        if (ii < N)
        {
            size += sizeof(double)*dims->nx[ii+1];  // pi
        }
    }

    make_int_multiple_of(64, &size);
    size += 1 * 64;

    return size;
}



ocp_nlp_out *assign_ocp_nlp_out(ocp_nlp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    int N = dims->N;

    ocp_nlp_out *out = (ocp_nlp_out *)c_ptr;
    c_ptr += sizeof(ocp_nlp_out);

    // double pointers
    assign_double_ptrs(N+1, &out->x, &c_ptr);
    assign_double_ptrs(N+1, &out->u, &c_ptr);
    assign_double_ptrs(N, &out->pi, &c_ptr);
    assign_double_ptrs(N+1, &out->lam, &c_ptr);

    // doubles
    align_char_to(64, &c_ptr);

    for (int ii = 0; ii <= N; ii++)
    {
        assign_double(dims->nx[ii], &out->x[ii], &c_ptr);
        assign_double(dims->nu[ii], &out->u[ii], &c_ptr);
        if (ii < N)
        {
            assign_double(dims->nx[ii+1], &out->pi[ii], &c_ptr);
        }
        assign_double(2*(dims->nb[ii] + dims->ng[ii] + dims->nh[ii]), &out->lam[ii], &c_ptr);
    }

    assert((char *)raw_memory + ocp_nlp_out_calculate_size(dims) >= c_ptr);

    return out;
}
