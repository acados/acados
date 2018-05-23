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

// blasfeo
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_blas.h"

// hpipm
#include "hpipm/include/hpipm_d_ocp_qp.h"
#include "hpipm/include/hpipm_d_ocp_qp_dim.h"
#include "hpipm/include/hpipm_d_ocp_qp_ipm.h"
#include "hpipm/include/hpipm_d_ocp_qp_kkt.h"
#include "hpipm/include/hpipm_d_ocp_qp_res.h"
#include "hpipm/include/hpipm_d_ocp_qp_sol.h"

// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/types.h"

/************************************************
 * config
 ************************************************/

int ocp_qp_solver_config_calculate_size()
{
    int size = 0;

    size += sizeof(qp_solver_config);

    return size;
}

qp_solver_config *ocp_qp_solver_config_assign(void *raw_memory)
{
    char *c_ptr = raw_memory;

    qp_solver_config *config = (qp_solver_config *) c_ptr;
    c_ptr += sizeof(qp_solver_config);

    return config;
}

int ocp_qp_xcond_solver_config_calculate_size()
{
    int size = 0;

    size += sizeof(ocp_qp_xcond_solver_config);

    size += ocp_qp_solver_config_calculate_size();  // qp solver

    return size;
}

ocp_qp_xcond_solver_config *ocp_qp_xcond_solver_config_assign(void *raw_memory)
{
    char *c_ptr = raw_memory;

    ocp_qp_xcond_solver_config *config = (ocp_qp_xcond_solver_config *) c_ptr;
    c_ptr += sizeof(ocp_qp_xcond_solver_config);

    config->qp_solver = ocp_qp_solver_config_assign(c_ptr);
    c_ptr += ocp_qp_solver_config_calculate_size();

    return config;
}

int ocp_qp_condensing_config_calculate_size()
{
    int size = 0;

    size += sizeof(ocp_qp_condensing_config);

    return size;
}

ocp_qp_condensing_config *ocp_qp_condensing_config_assign(void *raw_memory)
{
    char *c_ptr = raw_memory;

    ocp_qp_condensing_config *config = (ocp_qp_condensing_config *) c_ptr;
    c_ptr += sizeof(ocp_qp_condensing_config);

    return config;
}

/************************************************
 * dims
 ************************************************/

int ocp_qp_dims_calculate_size(int N)
{
    int size = sizeof(ocp_qp_dims);

    size += d_memsize_ocp_qp_dim(N);

    return size;
}

ocp_qp_dims *ocp_qp_dims_assign(int N, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_qp_dims *dims = (ocp_qp_dims *) c_ptr;
    c_ptr += sizeof(ocp_qp_dims);

    d_create_ocp_qp_dim(N, dims, c_ptr);
    c_ptr += d_memsize_ocp_qp_dim(N);

    dims->N = N;

    assert((char *) raw_memory + ocp_qp_dims_calculate_size(N) == c_ptr);

    return dims;
}

/************************************************
 * in
 ************************************************/

int ocp_qp_in_calculate_size(void *config, ocp_qp_dims *dims)
{
    int size = sizeof(ocp_qp_in);
    size += d_memsize_ocp_qp(dims);
    size += ocp_qp_dims_calculate_size(dims->N);  // TODO(all): remove !!!
    return size;
}

ocp_qp_in *ocp_qp_in_assign(void *config, ocp_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_qp_in *qp_in = (ocp_qp_in *) c_ptr;
    c_ptr += sizeof(ocp_qp_in);

    d_create_ocp_qp(dims, qp_in, c_ptr);
    c_ptr += d_memsize_ocp_qp(dims);

    ocp_qp_dims *dims_copy = ocp_qp_dims_assign(dims->N, c_ptr);  // TODO(all): remove !!!
    c_ptr += ocp_qp_dims_calculate_size(dims->N);                 // TODO(all): remove !!!

    dims_copy->N = dims->N;

    for (int ii = 0; ii < dims->N + 1; ii++)
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

    assert((char *) raw_memory + ocp_qp_in_calculate_size(config, dims) == c_ptr);

    return qp_in;
}

/************************************************
 * out
 ************************************************/

int ocp_qp_out_calculate_size(void *config, ocp_qp_dims *dims)
{
    int size = sizeof(ocp_qp_out);
    size += d_memsize_ocp_qp_sol(dims);
    size += ocp_qp_dims_calculate_size(dims->N);  // TODO(all): remove !!!
    size += sizeof(ocp_qp_info);
    return size;
}

ocp_qp_out *ocp_qp_out_assign(void *config, ocp_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_qp_out *qp_out = (ocp_qp_out *) c_ptr;
    c_ptr += sizeof(ocp_qp_out);

    d_create_ocp_qp_sol(dims, qp_out, c_ptr);
    c_ptr += d_memsize_ocp_qp_sol(dims);

    qp_out->misc = (void *) c_ptr;
    c_ptr += sizeof(ocp_qp_info);

    ocp_qp_dims *dims_copy = ocp_qp_dims_assign(dims->N, c_ptr);  // TODO(all): remove !!!
    c_ptr += ocp_qp_dims_calculate_size(dims->N);                 // TODO(all): remove !!!

    dims_copy->N = dims->N;

    for (int ii = 0; ii < dims->N + 1; ii++)
    {
        dims_copy->nx[ii] = dims->nx[ii];
        dims_copy->nu[ii] = dims->nu[ii];
        dims_copy->nb[ii] = dims->nb[ii];
        dims_copy->ng[ii] = dims->ng[ii];
        dims_copy->ns[ii] = dims->ns[ii];
        dims_copy->nbu[ii] = dims->nbu[ii];
        dims_copy->nbx[ii] = dims->nbx[ii];
    }

    qp_out->dim = dims_copy;

    assert((char *) raw_memory + ocp_qp_out_calculate_size(config, dims) == c_ptr);

    return qp_out;
}

/************************************************
 * res
 ************************************************/

// TODO(all): add config !!!

int ocp_qp_res_calculate_size(ocp_qp_dims *dims)
{
    int size = sizeof(ocp_qp_res);
    size += d_memsize_ocp_qp_res(dims);
    return size;
}

ocp_qp_res *ocp_qp_res_assign(ocp_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_qp_res *qp_res = (ocp_qp_res *) c_ptr;
    c_ptr += sizeof(ocp_qp_res);

    d_create_ocp_qp_res(dims, qp_res, c_ptr);
    c_ptr += d_memsize_ocp_qp_res(dims);

    assert((char *) raw_memory + ocp_qp_res_calculate_size(dims) == c_ptr);

    return qp_res;
}

int ocp_qp_res_workspace_calculate_size(ocp_qp_dims *dims)
{
    int size = sizeof(ocp_qp_res_ws);
    size += d_memsize_ocp_qp_res_workspace(dims);
    return size;
}

ocp_qp_res_ws *ocp_qp_res_workspace_assign(ocp_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_qp_res_ws *qp_res_ws = (ocp_qp_res_ws *) c_ptr;
    c_ptr += sizeof(ocp_qp_res_ws);

    d_create_ocp_qp_res_workspace(dims, qp_res_ws, c_ptr);
    c_ptr += d_memsize_ocp_qp_res_workspace(dims);

    assert((char *) raw_memory + ocp_qp_res_workspace_calculate_size(dims) == c_ptr);

    return qp_res_ws;
}

void ocp_qp_res_compute(ocp_qp_in *qp_in, ocp_qp_out *qp_out, ocp_qp_res *qp_res,
                        ocp_qp_res_ws *res_ws)
{
    // loop index
    int ii;

    //
    int N = qp_in->dim->N;
    int *nx = qp_in->dim->nx;
    int *nu = qp_in->dim->nu;
    int *nb = qp_in->dim->nb;
    int *ng = qp_in->dim->ng;

    struct blasfeo_dmat *DCt = qp_in->DCt;
    struct blasfeo_dvec *d = qp_in->d;
    int **idxb = qp_in->idxb;

    struct blasfeo_dvec *ux = qp_out->ux;
    struct blasfeo_dvec *t = qp_out->t;

    struct blasfeo_dvec *tmp_nbgM = res_ws->tmp_nbgM;

    int nx_i, nu_i, nb_i, ng_i;

    for (ii = 0; ii <= N; ii++)
    {
        nx_i = nx[ii];
        nu_i = nu[ii];
        nb_i = nb[ii];
        ng_i = ng[ii];

        // compute slacks for general constraints
        blasfeo_dgemv_t(nu_i + nx_i, ng_i, 1.0, DCt + ii, 0, 0, ux + ii, 0, -1.0, d + ii, nb_i,
                        t + ii, nb_i);
        blasfeo_dgemv_t(nu_i + nx_i, ng_i, -1.0, DCt + ii, 0, 0, ux + ii, 0, -1.0, d + ii,
                        2 * nb_i + ng_i, t + ii, 2 * nb_i + ng_i);

        // compute slacks for bounds
        blasfeo_dvecex_sp(nb_i, 1.0, idxb[ii], ux + ii, 0, tmp_nbgM + 0, 0);
        blasfeo_daxpby(nb_i, 1.0, tmp_nbgM + 0, 0, -1.0, d + ii, 0, t + ii, 0);
        blasfeo_daxpby(nb_i, -1.0, tmp_nbgM + 0, 0, -1.0, d + ii, nb_i + ng_i, t + ii, nb_i + ng_i);
    }

    d_compute_res_ocp_qp(qp_in, qp_out, qp_res, res_ws);

    return;
}

void ocp_qp_res_compute_nrm_inf(ocp_qp_res *qp_res, double res[4])
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

#if 1

    double tmp;

    res[0] = 0.0;
    for (ii = 0; ii <= N; ii++)
    {
        blasfeo_dvecnrm_inf(nx[ii] + nu[ii] + 2 * ns[ii], &qp_res->res_g[ii], 0, &tmp);
        res[0] = tmp > res[0] ? tmp : res[0];
    }

    res[1] = 0.0;
    for (ii = 0; ii < N; ii++)
    {
        blasfeo_dvecnrm_inf(nx[ii + 1], &qp_res->res_b[ii], 0, &tmp);
        res[1] = tmp > res[1] ? tmp : res[1];
    }

    res[2] = 0.0;
    for (ii = 0; ii <= N; ii++)
    {
        blasfeo_dvecnrm_inf(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], &qp_res->res_d[ii], 0, &tmp);
        res[2] = tmp > res[2] ? tmp : res[2];
    }

    res[3] = 0.0;
    for (ii = 0; ii <= N; ii++)
    {
        blasfeo_dvecnrm_inf(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], &qp_res->res_m[ii], 0, &tmp);
        res[3] = tmp > res[3] ? tmp : res[3];
    }

#else

    // XXX this should be avoided, since it does employ strucutre of the HPIPM core that may change
    // !!!

    // compute size of res_q, res_b, res_d and res_m
    int nvt = 0;
    int net = 0;
    int nct = 0;

    for (ii = 0; ii < N; ii++)
    {
        nvt += nx[ii] + nu[ii] + 2 * ns[ii];
        net += nx[ii + 1];
        nct += 2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii];
    }
    nvt += nx[ii] + nu[ii] + 2 * ns[ii];
    nct += 2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii];

    // compute infinity norms of residuals
    blasfeo_dvecnrm_inf(nvt, qp_res->res_g, 0, &res[0]);
    blasfeo_dvecnrm_inf(net, qp_res->res_b, 0, &res[1]);
    blasfeo_dvecnrm_inf(nct, qp_res->res_d, 0, &res[2]);
    blasfeo_dvecnrm_inf(nct, qp_res->res_m, 0, &res[3]);

#endif

    return;
}
