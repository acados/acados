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
#include <assert.h>
// hpipm
#include "hpipm/include/hpipm_d_dense_qp.h"
#include "hpipm/include/hpipm_d_dense_qp_sol.h"
#include "hpipm/include/hpipm_d_dense_qp_res.h"
#include "hpipm/include/hpipm_d_dense_qp_ipm.h"
#include "hpipm/include/hpipm_d_dense_qp_kkt.h"
// blasfeo
#include "blasfeo_target.h"
#include "blasfeo_common.h"
#include "blasfeo_d_aux.h"
#include "blasfeo_d_blas.h"
// acados
#include "acados/utils/types.h"
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/dense_qp/dense_qp_hpipm.h"
#include "acados/dense_qp/dense_qp_qpoases.h"




int dense_qp_in_calculate_size(dense_qp_dims *dims)
{
    int size = sizeof(dense_qp_in);
    size += sizeof(dense_qp_dims);
    size += d_memsize_dense_qp(dims);
    return size;
}



dense_qp_in *assign_dense_qp_in(dense_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *)raw_memory;

    dense_qp_in *qp_in = (dense_qp_in *) c_ptr;
    c_ptr += sizeof(dense_qp_in);

    assert((size_t)c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    d_create_dense_qp(dims, qp_in, c_ptr);
    c_ptr += d_memsize_dense_qp(dims);

    qp_in->dim = (dense_qp_dims *) c_ptr;
    c_ptr += sizeof(dense_qp_dims);

    qp_in->dim->nv = dims->nv;
    qp_in->dim->ne = dims->ne;
    qp_in->dim->nb = dims->nb;
    qp_in->dim->ng = dims->ng;
    qp_in->dim->ns = dims->ns;

    assert((char*)raw_memory + dense_qp_in_calculate_size(dims) == c_ptr);

    return qp_in;
}



int dense_qp_out_calculate_size(dense_qp_dims *dims)
{
    int size = sizeof(dense_qp_out);
    size += d_memsize_dense_qp_sol(dims);
    size += sizeof(dense_qp_info);
    return size;
}



dense_qp_out *assign_dense_qp_out(dense_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    dense_qp_out *qp_out = (dense_qp_out *) c_ptr;
    c_ptr += sizeof(dense_qp_out);

    assert((size_t)c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    d_create_dense_qp_sol(dims, qp_out, c_ptr);
    c_ptr += d_memsize_dense_qp_sol(dims);

    qp_out->misc = (void *) c_ptr;
    c_ptr += sizeof(dense_qp_info);

    assert((char*)raw_memory + dense_qp_out_calculate_size(dims) == c_ptr);

    return qp_out;
}



int dense_qp_res_calculate_size(dense_qp_dims *dims)
{
    int size = sizeof(dense_qp_res);
    size += d_memsize_dense_qp_res(dims);
    return size;
}



dense_qp_res *assign_dense_qp_res(dense_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    dense_qp_res *qp_res = (dense_qp_res *) c_ptr;
    c_ptr += sizeof(dense_qp_res);

    d_create_dense_qp_res(dims, qp_res, c_ptr);
    c_ptr += d_memsize_dense_qp_res(dims);

    assert((char*) raw_memory + dense_qp_res_calculate_size(dims) == c_ptr);

    return qp_res;
}



int dense_qp_res_ws_calculate_size(dense_qp_dims *dims)
{
    int size = sizeof(dense_qp_res_ws);
    size += d_memsize_dense_qp_res_workspace(dims);
    return size;
}



dense_qp_res_ws *assign_dense_qp_res_ws(dense_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    dense_qp_res_ws *res_ws = (dense_qp_res_ws *) c_ptr;
    c_ptr += sizeof(dense_qp_res_ws);

    d_create_dense_qp_res_workspace(dims, res_ws, c_ptr);
    c_ptr += d_memsize_dense_qp_res_workspace(dims);

    assert((char*) raw_memory + dense_qp_res_ws_calculate_size(dims) == c_ptr);

    return res_ws;
}



void compute_dense_qp_res(dense_qp_in *qp_in, dense_qp_out *qp_out, dense_qp_res *qp_res, dense_qp_res_ws *res_ws)
{
    int nvd = qp_in->dim->nv;
    int ned = qp_in->dim->ne;
    int ngd = qp_in->dim->ng;
    int nbd = qp_in->dim->nb;

    int *idxb = qp_in->idxb;

    struct d_strvec *tmp_nbg = res_ws->tmp_nbg;

 	// set t to zero
    dvecse_libstr(2*nbd+2*ngd, 0.0, qp_out->t, 0);

    // compute slacks for general constraints
    dgemv_t_libstr(nvd, ngd, 1.0, qp_in->Ct, 0, 0, qp_out->v, 0, -1.0, qp_in->d, nbd, qp_out->t, nbd);
    dgemv_t_libstr(nvd, ngd, -1.0, qp_in->Ct, 0, 0, qp_out->v, 0, 1.0, qp_in->d, 2*nbd+ngd, qp_out->t, 2*nbd+ngd);

    // compute slacks for bounds
    dvecex_sp_libstr(nbd, 1.0, idxb, qp_out->v, 0, tmp_nbg+0, 0);
    daxpy_libstr(nbd, -1.0, qp_in->d, 0, qp_out->t, 0, qp_out->t, 0);
    daxpy_libstr(nbd, 1.0, qp_in->d, nbd+ngd, qp_out->t, nbd+ngd, qp_out->t, nbd+ngd);
    daxpy_libstr(nbd, 1.0, tmp_nbg+0, 0, qp_out->t, 0, qp_out->t, 0);
    daxpy_libstr(nbd, -1.0, tmp_nbg+0, 0, qp_out->t, nbd+ngd, qp_out->t, nbd+ngd);

    // compute residuals
    d_compute_res_dense_qp(qp_in, qp_out, qp_res, res_ws);
}



void compute_dense_qp_res_nrm_inf(dense_qp_res *qp_res, double res[4])
{
    int nv = qp_res->dim->nv;
    int nb = qp_res->dim->nb;
    int ne = qp_res->dim->ne;
    int ng = qp_res->dim->ng;
    int ns = qp_res->dim->ns;

    dvecnrm_inf_libstr(nv+2*ns, qp_res->res_g, 0, &res[0]);
    dvecnrm_inf_libstr(ne, qp_res->res_b, 0, &res[1]);
    dvecnrm_inf_libstr(2*nb+2*ng+2*ns, qp_res->res_d, 0, &res[2]);
    dvecnrm_inf_libstr(2*nb+2*ng+2*ns, qp_res->res_m, 0, &res[3]);
}
