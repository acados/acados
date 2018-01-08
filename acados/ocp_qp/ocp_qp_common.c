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

// hpipm
#include "hpipm/include/hpipm_d_ocp_qp.h"
#include "hpipm/include/hpipm_d_ocp_qp_sol.h"
#include "hpipm/include/hpipm_d_ocp_qp_dim.h"
#include "hpipm/include/hpipm_d_ocp_qp_res.h"
#include "hpipm/include/hpipm_d_ocp_qp_ipm.h"
#include "hpipm/include/hpipm_d_ocp_qp_kkt.h"
// acados
#include "acados/utils/types.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_condensing_solver.h"
#include "acados/ocp_qp/ocp_qp_sparse_solver.h"
#include "acados/ocp_qp/ocp_qp_hpipm.h"
#include "acados/ocp_qp/ocp_qp_hpmpc.h"
#include "acados/dense_qp/dense_qp_hpipm.h"
#include "acados/dense_qp/dense_qp_qpoases.h"
#include "acados/dense_qp/dense_qp_qore.h"
// blasfeo
#include "blasfeo_target.h"
#include "blasfeo_common.h"
#include "blasfeo_d_aux.h"
#include "blasfeo_d_blas.h"



int ocp_qp_in_calculate_size(ocp_qp_dims *dims)
{
    int size = sizeof(ocp_qp_in);
    size += d_memsize_ocp_qp(dims);
    size += sizeof(ocp_qp_dims);
    size += d_memsize_ocp_qp_dim(dims->N);
    return size;
}



ocp_qp_in *assign_ocp_qp_in(ocp_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_qp_in *qp_in = (ocp_qp_in *) c_ptr;
    c_ptr += sizeof(ocp_qp_in);

    ocp_qp_dims *dims_copy = (ocp_qp_dims *) c_ptr;
    c_ptr += sizeof(ocp_qp_dims);

    d_create_ocp_qp(dims, qp_in, c_ptr);
    c_ptr += d_memsize_ocp_qp(dims);

    d_create_ocp_qp_dim(dims->N, dims_copy, c_ptr);
    c_ptr += d_memsize_ocp_qp_dim(dims->N);

    dims_copy->N = dims->N;

    for (int ii = 0; ii < dims->N+1; ii++)
    {
        dims_copy->nx[ii] = dims->nx[ii];
        dims_copy->nu[ii] = dims->nu[ii];
        dims_copy->nb[ii] = dims->nb[ii];
        dims_copy->ng[ii] = dims->ng[ii];
        dims_copy->ns[ii] = dims->ns[ii];
        dims_copy->nbu[ii] = dims->nbu[ii];
        dims_copy->nbx[ii] = dims->nbx[ii];
    }

    qp_in->dim = dims_copy;

    assert((char*) raw_memory + ocp_qp_in_calculate_size(dims) == c_ptr);

    return qp_in;
}



int ocp_qp_out_calculate_size(ocp_qp_dims *dims)
{
    int size = sizeof(ocp_qp_out);
    size += d_memsize_ocp_qp_sol(dims);
    size += sizeof(ocp_qp_info);
    return size;
}



ocp_qp_out *assign_ocp_qp_out(ocp_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_qp_out *qp_out = (ocp_qp_out *) c_ptr;
    c_ptr += sizeof(ocp_qp_out);

    d_create_ocp_qp_sol(dims, qp_out, c_ptr);
    c_ptr += d_memsize_ocp_qp_sol(dims);

    qp_out->misc = (void *) c_ptr;
    c_ptr += sizeof(ocp_qp_info);

    assert((char*) raw_memory + ocp_qp_out_calculate_size(dims) == c_ptr);

    return qp_out;
}



int ocp_qp_res_calculate_size(ocp_qp_dims *dims)
{
    int size = sizeof(ocp_qp_res);
    size += d_memsize_ocp_qp_res(dims);
    return size;
}



ocp_qp_res *assign_ocp_qp_res(ocp_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_qp_res *qp_res = (ocp_qp_res *) c_ptr;
    c_ptr += sizeof(ocp_qp_res);

    d_create_ocp_qp_res(dims, qp_res, c_ptr);
    c_ptr += d_memsize_ocp_qp_res(dims);

    assert((char*) raw_memory + ocp_qp_res_calculate_size(dims) == c_ptr);

    return qp_res;
}



int ocp_qp_res_ws_calculate_size(ocp_qp_dims *dims)
{
    int size = sizeof(ocp_qp_res_ws);
    size += d_memsize_ocp_qp_res_workspace(dims);
    return size;
}



ocp_qp_res_ws *assign_ocp_qp_res_ws(ocp_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_qp_res_ws *qp_res_ws = (ocp_qp_res_ws *) c_ptr;
    c_ptr += sizeof(ocp_qp_res_ws);

    d_create_ocp_qp_res_workspace(dims, qp_res_ws, c_ptr);
    c_ptr += d_memsize_ocp_qp_res_workspace(dims);

    assert((char*) raw_memory + ocp_qp_res_ws_calculate_size(dims) == c_ptr);

    return qp_res_ws;
}



void compute_ocp_qp_res(ocp_qp_in *qp_in, ocp_qp_out *qp_out, ocp_qp_res *qp_res, ocp_qp_res_ws *res_ws)
{
    // loop index
	int ii;

	//
	int N = qp_in->dim->N;
	int *nx = qp_in->dim->nx;
	int *nu = qp_in->dim->nu;
	int *nb = qp_in->dim->nb;
    int *ng = qp_in->dim->ng;

    struct d_strmat *DCt = qp_in->DCt;
    struct d_strvec *d = qp_in->d;
    int **idxb = qp_in->idxb;

    struct d_strvec *ux = qp_out->ux;
    struct d_strvec *t = qp_out->t;

    struct d_strvec *tmp_nbgM = res_ws->tmp_nbgM;

    int nx_i, nu_i, nb_i, ng_i;

    for(ii=0; ii<=N; ii++)
    {
        nx_i = nx[ii];
		nu_i = nu[ii];
		nb_i = nb[ii];
		ng_i = ng[ii];

        // set t to zero
        dvecse_libstr(2*nb_i+2*ng_i, 0.0, t+ii, 0);

        // compute slacks for general constraints
        dgemv_t_libstr(nu_i+nx_i, ng_i,  1.0, DCt+ii, 0, 0, ux+ii, 0, -1.0, d+ii, nb_i, t+ii, nb_i);
        dgemv_t_libstr(nu_i+nx_i, ng_i, -1.0, DCt+ii, 0, 0, ux+ii, 0,  1.0, d+ii, 2*nb_i+ng_i, t+ii, 2*nb_i+ng_i);

        // compute slacks for bounds
        dvecex_sp_libstr(nb_i, 1.0, idxb[ii], ux+ii, 0, tmp_nbgM+0, 0);
        daxpy_libstr(nb_i, -1.0, d+ii, 0, t+ii, 0, t+ii, 0);
        daxpy_libstr(nb_i,  1.0, d+ii, nb_i+ng_i, t+ii, nb_i+ng_i, t+ii, nb_i+ng_i);
        daxpy_libstr(nb_i,  1.0, tmp_nbgM+0, 0, t+ii, 0, t+ii, 0);
        daxpy_libstr(nb_i, -1.0, tmp_nbgM+0, 0, t+ii, nb_i+ng_i, t+ii, nb_i+ng_i);
    }

    d_compute_res_ocp_qp(qp_in, qp_out, qp_res, res_ws);
}



void compute_ocp_qp_res_nrm_inf(ocp_qp_res *qp_res, double res[4])
{
    // loop index
	int ii;

	// extract ocp qp size
	int N = qp_res->dim->N;
	int *nx = qp_res->dim->nx;
	int *nu = qp_res->dim->nu;
	int *nb = qp_res->dim->nb;
	int *ng = qp_res->dim->ng;
    int *ns = qp_res->dim->ns;

    // compute size of res_q, res_b, res_d and res_m
    int nvt = 0;
	int net = 0;
    int nct = 0;

    for(ii=0; ii<N; ii++)
	{
		nvt += nx[ii]+nu[ii]+2*ns[ii];
		net += nx[ii+1];
		nct += 2*nb[ii]+2*ng[ii]+2*ns[ii];
	}
	nvt += nx[ii]+nu[ii]+2*ns[ii];
    nct += 2*nb[ii]+2*ng[ii]+2*ns[ii];

    // compute infinity norms of residuals
    dvecnrm_inf_libstr(nvt, qp_res->res_g, 0, &res[0]);
	dvecnrm_inf_libstr(net, qp_res->res_b, 0, &res[1]);
	dvecnrm_inf_libstr(nct, qp_res->res_d, 0, &res[2]);
	dvecnrm_inf_libstr(nct, qp_res->res_m, 0, &res[3]);
}



// void form_nbu_nbx_rev(int N, int *nbu, int *nbx, int *nb, int* nx, int *nu, int **idxb_rev)
// {
//     for (int ii = 0; ii < N+1; ii++)
//     {
//         nbu[ii] = 0;
//         nbx[ii] = 0;
//         for (int jj = 0; jj < nb[ii]; jj++)
//         {
//             if (idxb_rev[ii][jj] < nx[ii])
//             {
//                 nbx[ii]++;
//             }
//             else
//             {
//                 nbu[ii]++;
//             }
//         }
//     }
// }



int set_qp_solver_fun_ptrs(qp_solver_t qp_solver_name, void *qp_solver)
{
    int return_value = ACADOS_SUCCESS;

    switch (qp_solver_name)
    {
        case HPIPM:
            ((ocp_qp_solver *)qp_solver)->calculate_args_size = &ocp_qp_hpipm_calculate_args_size;
            ((ocp_qp_solver *)qp_solver)->assign_args = &ocp_qp_hpipm_assign_args;
            ((ocp_qp_solver *)qp_solver)->initialize_default_args = &ocp_qp_hpipm_initialize_default_args;
            ((ocp_qp_solver *)qp_solver)->calculate_memory_size = &ocp_qp_hpipm_calculate_memory_size;
            ((ocp_qp_solver *)qp_solver)->assign_memory = &ocp_qp_hpipm_assign_memory;
            ((ocp_qp_solver *)qp_solver)->calculate_workspace_size = &ocp_qp_hpipm_calculate_workspace_size;
            ((ocp_qp_solver *)qp_solver)->fun = &ocp_qp_hpipm;
            break;
        case HPMPC:
            ((ocp_qp_solver *)qp_solver)->calculate_args_size = &ocp_qp_hpmpc_calculate_args_size;
            ((ocp_qp_solver *)qp_solver)->assign_args = &ocp_qp_hpmpc_assign_args;
            ((ocp_qp_solver *)qp_solver)->initialize_default_args = &ocp_qp_hpmpc_initialize_default_args;
            ((ocp_qp_solver *)qp_solver)->calculate_memory_size = &ocp_qp_hpmpc_calculate_memory_size;
            ((ocp_qp_solver *)qp_solver)->assign_memory = &ocp_qp_hpmpc_assign_memory;
            ((ocp_qp_solver *)qp_solver)->calculate_workspace_size = &ocp_qp_hpmpc_calculate_workspace_size;
            ((ocp_qp_solver *)qp_solver)->fun = &ocp_qp_hpmpc;
            break;
        case CONDENSING_HPIPM:
            ((dense_qp_solver *)qp_solver)->calculate_args_size = &dense_qp_hpipm_calculate_args_size;
            ((dense_qp_solver *)qp_solver)->assign_args = &dense_qp_hpipm_assign_args;
            ((dense_qp_solver *)qp_solver)->initialize_default_args = &dense_qp_hpipm_initialize_default_args;
            ((dense_qp_solver *)qp_solver)->calculate_memory_size = &dense_qp_hpipm_calculate_memory_size;
            ((dense_qp_solver *)qp_solver)->assign_memory = &dense_qp_hpipm_assign_memory;
            ((dense_qp_solver *)qp_solver)->calculate_workspace_size = &dense_qp_hpipm_calculate_workspace_size;
            ((dense_qp_solver *)qp_solver)->fun = &dense_qp_hpipm;
            break;
        case CONDENSING_QPOASES:
            ((dense_qp_solver *)qp_solver)->calculate_args_size = &dense_qp_qpoases_calculate_args_size;
            ((dense_qp_solver *)qp_solver)->assign_args = &dense_qp_qpoases_assign_args;
            ((dense_qp_solver *)qp_solver)->initialize_default_args = &dense_qp_qpoases_initialize_default_args;
            ((dense_qp_solver *)qp_solver)->calculate_memory_size = &dense_qp_qpoases_calculate_memory_size;
            ((dense_qp_solver *)qp_solver)->assign_memory = &dense_qp_qpoases_assign_memory;
            ((dense_qp_solver *)qp_solver)->calculate_workspace_size = &dense_qp_qpoases_calculate_workspace_size;
            ((dense_qp_solver *)qp_solver)->fun = &dense_qp_qpoases;
            break;
        case CONDENSING_QORE:
            ((dense_qp_solver *)qp_solver)->calculate_args_size = &dense_qp_qore_calculate_args_size;
            ((dense_qp_solver *)qp_solver)->assign_args = &dense_qp_qore_assign_args;
            ((dense_qp_solver *)qp_solver)->initialize_default_args = &dense_qp_qore_initialize_default_args;
            ((dense_qp_solver *)qp_solver)->calculate_memory_size = &dense_qp_qore_calculate_memory_size;
            ((dense_qp_solver *)qp_solver)->assign_memory = &dense_qp_qore_assign_memory;
            ((dense_qp_solver *)qp_solver)->calculate_workspace_size = &dense_qp_qore_calculate_workspace_size;
            ((dense_qp_solver *)qp_solver)->fun = &dense_qp_qore;
        default:
            return_value = ACADOS_FAILURE;
    }
    return return_value;
}



void set_xcond_qp_solver_fun_ptrs(qp_solver_t qp_solver_name, ocp_qp_xcond_solver *qp_solver)
{
    if (qp_solver_name < CONDENSING_HPIPM)
    {
        qp_solver->calculate_args_size = &ocp_qp_sparse_solver_calculate_args_size;
        qp_solver->assign_args = &ocp_qp_sparse_solver_assign_args;
        qp_solver->initialize_default_args = &ocp_qp_sparse_solver_initialize_default_args;
        qp_solver->calculate_memory_size = &ocp_qp_sparse_solver_calculate_memory_size;
        qp_solver->assign_memory = &ocp_qp_sparse_solver_assign_memory;
        qp_solver->calculate_workspace_size = &ocp_qp_sparse_solver_calculate_workspace_size;
        qp_solver->fun = &ocp_qp_sparse_solver;
    } else
    {
        qp_solver->calculate_args_size = &ocp_qp_condensing_solver_calculate_args_size;
        qp_solver->assign_args = &ocp_qp_condensing_solver_assign_args;
        qp_solver->initialize_default_args = &ocp_qp_condensing_solver_initialize_default_args;
        qp_solver->calculate_memory_size = &ocp_qp_condensing_solver_calculate_memory_size;
        qp_solver->assign_memory = &ocp_qp_condensing_solver_assign_memory;
        qp_solver->calculate_workspace_size = &ocp_qp_condensing_solver_calculate_workspace_size;
        qp_solver->fun = &ocp_qp_condensing_solver;
    }
    set_qp_solver_fun_ptrs(qp_solver_name, qp_solver->qp_solver_funs);
}
