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
#include <string.h>

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


void ocp_qp_dims_set(void *config_, void *dims_, int stage, const char *field, const int* value)
{
    ocp_qp_dims *dims = (ocp_qp_dims *) dims_;
    // char field_copy[MAX_STR_LEN];
    char *field_copy = (char *) field;

    d_set_ocp_qp_dim(field_copy, stage, *value, dims);
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
        dims_copy->nsbu[ii] = dims->nsbu[ii];
        dims_copy->nsbx[ii] = dims->nsbx[ii];
        dims_copy->nsg[ii] = dims->nsg[ii];
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
        dims_copy->nsbx[ii] = dims->nsbx[ii];
        dims_copy->nsbu[ii] = dims->nsbu[ii];
        dims_copy->nsg[ii] = dims->nsg[ii];
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
    ocp_qp_info *info = (ocp_qp_info *) qp_out->misc;

    if (info->t_computed == 0)
    {
        ocp_qp_compute_t(qp_in, qp_out);
        info->t_computed = 1;
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

void ocp_qp_stack_slacks_dims(ocp_qp_dims *in, ocp_qp_dims *out)
{
    int N     = in->N;
    int *nx   = in->nx;
    int *nu   = in->nu;
    int *nb   = in->nb;
    int *nbx  = in->nbx;
    int *nbu  = in->nbu;
    int *ng   = in->ng;
    int *ns   = in->ns;
    int *nsbx = in->nsbx;
    int *nsbu = in->nsbu;
    // int *nsg  = in->nsg;

    out->N = N;

    for (int ii = 0; ii <= N; ii++)
    {
        out->nx[ii] = nx[ii];
        out->nu[ii] = nu[ii] + 2*ns[ii];
        out->nb[ii] = nb[ii] - nsbx[ii] - nsbu[ii] + 2*ns[ii];
        out->nbx[ii] = nbx[ii] - nsbx[ii];
        out->nbu[ii] = nbu[ii] - nsbu[ii] + 2*ns[ii];
        out->ng[ii] = ns[ii] > 0 ? ng[ii] + nsbx[ii] + nsbu[ii] : ng[ii];
        out->ns[ii] = 0;
        out->nsbx[ii] = 0;
        out->nsbu[ii] = 0;
        out->nsg[ii] = 0;
    }
}

void ocp_qp_stack_slacks(ocp_qp_in *in, ocp_qp_in *out)
{
    int N     = in->dim->N;
    int *nx   = in->dim->nx;
    int *nu   = in->dim->nu;
    int *nb   = in->dim->nb;
    // int *nbx  = in->dim->nbx;
    // int *nbu  = in->dim->nbu;
    int *ng   = in->dim->ng;
    int *ns   = in->dim->ns;
    int *nsbx = in->dim->nsbx;
    int *nsbu = in->dim->nsbu;
    // int *nsg  = in->dim->nsg;
    int **idxb = in->idxb;
    int **idxs = in->idxs;

    int *nx2  = out->dim->nx;
    int *nu2  = out->dim->nu;
    int *nb2  = out->dim->nb;
    // int *nbx2 = out->dim->nbx;
    // int *nbu2 = out->dim->nbu;
    int *ng2  = out->dim->ng;

    UNUSED(nsbx);
    UNUSED(nsbu);

    for (int ii = 0; ii <= N; ii++)
    {
        if (ii < N)
        {
            // set matrices to 0.0
            blasfeo_dgese(nu2[ii]+nx2[ii]+1, nx2[ii+1], 0.0, out->BAbt+ii, 0, 0);

            // copy in->BAbt to out->BAbt, out->BAbt = [0; B'; A'; b']
            blasfeo_dgecp(nu[ii]+nx[ii]+1, nx[ii+1], in->BAbt+ii, 0, 0, out->BAbt+ii, 2*ns[ii], 0);

            // copy in->b to out->b
            blasfeo_dveccp(nx2[ii+1], in->b+ii, 0, out->b+ii, 0);
        }

        // set matrices to 0.0
        blasfeo_dgese(nu2[ii]+nx2[ii]+1, nu2[ii]+nx2[ii], 0.0, out->RSQrq+ii, 0, 0);
        blasfeo_dgese(nu2[ii]+nx2[ii], ng2[ii], 0.0, out->DCt+ii, 0, 0);

        // copy in->Z to the main diagonal of out->RSQrq, out->RSQrq = [Z 0 0; 0 0 0; 0 0 0; 0 0 0]
        blasfeo_ddiain(2*ns[ii], 1.0, in->Z+ii, 0, out->RSQrq+ii, 0, 0);

        // copy in->RSQrq to out->RSQrq, out->RSQrq = [Z 0 0; 0 R S; 0 S' Q; 0 r' q']
        blasfeo_dgecp(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], in->RSQrq+ii, 0, 0, out->RSQrq+ii, 2*ns[ii],
            2*ns[ii]);

        // copy in->rqz to out->rqz, out->rqz = [z r q]
        blasfeo_dveccp(2*ns[ii], in->rqz+ii, nu[ii]+nx[ii], out->rqz+ii, 0);
        blasfeo_dveccp(nu[ii]+nx[ii], in->rqz+ii, 0, out->rqz+ii, 2*ns[ii]);

        // copy in->DCt to out->DCt, out->DCt = [0 0; DCt 0]
        blasfeo_dgecp(nu[ii]+nx[ii], ng[ii], in->DCt+ii, 0, 0, out->DCt+ii, 2*ns[ii], 0);

        if (ns[ii] > 0)
        {
            // set flags for non-softened box constraints
            // use out->m temporarily for this
            for (int jj = 0; jj < nb[ii]; jj++)
            {
                // TODO(bnovoselnik): pick up some workspace for this
                BLASFEO_DVECEL(out->m+ii, jj) = 1.0;
            }

            int col_b = ng[ii];
            for (int jj = 0; jj < ns[ii]; jj++)
            {
                int js = idxs[ii][jj];

                int idx_v_ls1 = jj;
                int idx_v_us1 = ns[ii]+jj;

                int idx_d_ls0 = js;
                int idx_d_us0 = nb[ii]+ng[ii]+js;
                int idx_d_ls1;
                int idx_d_us1;

                if (js < nb[ii])  // soft box constraint
                {
                    // index of a soft box constraint
                    int jv = idxb[ii][js]+2*ns[ii];

                    idx_d_ls1 = nb2[ii]+col_b;
                    idx_d_us1 = 2*nb2[ii]+ng2[ii]+col_b;

                    // soft box constraint, set its flag to -1
                    BLASFEO_DVECEL(out->m+ii, js) = -1.0;

                    // insert soft box constraint into out->DCt, lb <= ux + sl - su <= ub
                    BLASFEO_DMATEL(out->DCt+ii, jv, col_b) = 1.0;
                    BLASFEO_DMATEL(out->DCt+ii, idx_v_ls1, col_b) = +1.0;
                    BLASFEO_DMATEL(out->DCt+ii, idx_v_us1, col_b) = -1.0;
                    BLASFEO_DVECEL(out->d+ii, idx_d_ls1) = BLASFEO_DVECEL(in->d+ii, idx_d_ls0);
                    BLASFEO_DVECEL(out->d+ii, idx_d_us1) = -BLASFEO_DVECEL(in->d+ii, idx_d_us0);

                    col_b++;
                }
                else  // soft general constraint
                {
                    // index of a soft general constraint
                    int col_g = js - nb[ii];

                    // soft general constraint, lg <= D u + C x + sl - su <= ug
                    BLASFEO_DMATEL(out->DCt+ii, idx_v_ls1, col_g) = +1.0;
                    BLASFEO_DMATEL(out->DCt+ii, idx_v_us1, col_g) = -1.0;
                }

                // slacks have box constraints
                out->idxb[ii][idx_v_ls1] = idx_v_ls1;
                out->idxb[ii][idx_v_us1] = idx_v_us1;
            }

            int k_nsb = 0;
            for (int jj = 0; jj < nb[ii]; jj++)
            {
                if (BLASFEO_DVECEL(out->m+ii, jj) > 0)
                {
                    // copy nonsoftened box constraint bounds to out->d
                    BLASFEO_DVECEL(out->d+ii, 2*ns[ii]+k_nsb) = BLASFEO_DVECEL(in->d+ii, jj);
                    BLASFEO_DVECEL(out->d+ii, nb2[ii]+ng2[ii]+2*ns[ii]+k_nsb) =
                                                        -BLASFEO_DVECEL(in->d+ii, nb[ii]+ng[ii]+jj);
                    out->idxb[ii][2*ns[ii]+k_nsb] = 2*ns[ii] + jj;
                    k_nsb++;
                }
            }

            assert(k_nsb == nb[ii] - nsbu[ii] - nsbx[ii] && "Dimensions are wrong!");

            // copy ls and us to out->lb
            blasfeo_dveccp(2*ns[ii], in->d+ii, 2*nb[ii]+2*ng[ii], out->d+ii, 0);

            // for slacks out->ub is +INFTY
            blasfeo_dvecse(2*ns[ii], 1.0e6, out->d+ii, nb2[ii]+ng2[ii]);

            // copy in->lg to out->lg
            blasfeo_dveccp(ng[ii], in->d+ii, nb[ii], out->d+ii, nb2[ii]);

            // copy in->ug to out->ug
            blasfeo_dveccpsc(ng[ii], -1.0, in->d+ii, 2*nb[ii]+ng[ii], out->d+ii, 2*nb2[ii]+ng2[ii]);

            // flip signs for ub and ug
            blasfeo_dvecsc(nb2[ii]+ng2[ii], -1.0, out->d+ii, nb2[ii]+ng2[ii]);

            // set out->m to 0.0
            blasfeo_dvecse(2*nb2[ii]+2*ng2[ii], 0.0, out->m+ii, 0);
        }
        else
        {
            blasfeo_dveccp(2*nb[ii]+2*ng[ii], in->d+ii, 0, out->d+ii, 0);
            blasfeo_dveccp(2*nb[ii]+2*ng[ii], in->m+ii, 0, out->m+ii, 0);
            for (int jj = 0; jj < nb[ii]; jj++) out->idxb[ii][jj] = in->idxb[ii][jj];
        }
    }
}

void ocp_qp_compute_t(ocp_qp_in *qp_in, ocp_qp_out *qp_out)
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

    int nx_i, nu_i, nb_i, ng_i;

    for (ii = 0; ii <= N; ii++)
    {
        nx_i = nx[ii];
        nu_i = nu[ii];
        nb_i = nb[ii];
        ng_i = ng[ii];

        // compute slacks for bounds
        blasfeo_dvecex_sp(nb_i, 1.0, idxb[ii], ux + ii, 0, t+ii, nb_i + ng_i);
        blasfeo_daxpby(nb_i, 1.0, t + ii, nb_i + ng_i, -1.0, d + ii, 0, t + ii, 0);
        blasfeo_daxpby(nb_i, -1.0, t + ii, nb_i + ng_i, -1.0, d + ii, nb_i + ng_i, t + ii,
            nb_i + ng_i);

        // compute slacks for general constraints
        blasfeo_dgemv_t(nu_i + nx_i, ng_i, 1.0, DCt + ii, 0, 0, ux + ii, 0, -1.0, d + ii, nb_i,
                        t + ii, nb_i);
        blasfeo_dgemv_t(nu_i + nx_i, ng_i, -1.0, DCt + ii, 0, 0, ux + ii, 0, -1.0, d + ii,
                        2 * nb_i + ng_i, t + ii, 2 * nb_i + ng_i);
    }
}
