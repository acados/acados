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
#include <stddef.h>
// hpipm
#include "hpipm/include/hpipm_d_dense_qp.h"
#include "hpipm/include/hpipm_d_dense_qp_ipm.h"
#include "hpipm/include/hpipm_d_dense_qp_kkt.h"
#include "hpipm/include/hpipm_d_dense_qp_res.h"
#include "hpipm/include/hpipm_d_dense_qp_sol.h"
// blasfeo
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_blas.h"
// acados
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/utils/types.h"

/************************************************
 * config
 ************************************************/

int dense_qp_solver_config_calculate_size()
{
    int size = 0;

    size += sizeof(qp_solver_config);

    return size;
}

qp_solver_config *dense_qp_solver_config_assign(void *raw_memory)
{
    char *c_ptr = raw_memory;

    qp_solver_config *config = (qp_solver_config *) c_ptr;
    c_ptr += sizeof(qp_solver_config);

    return config;
}

/************************************************
 * dims
 ************************************************/

int dense_qp_dims_calculate_size()
{
    int size = sizeof(dense_qp_dims);

    size += d_memsize_dense_qp_dim();

    return size;
}

dense_qp_dims *dense_qp_dims_assign(void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    dense_qp_dims *dims = (dense_qp_dims *) c_ptr;
    c_ptr += sizeof(dense_qp_dims);

    d_create_dense_qp_dim(dims, c_ptr);
    c_ptr += d_memsize_dense_qp_dim();

    assert((char *) raw_memory + dense_qp_dims_calculate_size() == c_ptr);

    return dims;
}

/************************************************
 * in
 ************************************************/

int dense_qp_in_calculate_size(void *config, dense_qp_dims *dims)
{
    int size = sizeof(dense_qp_in);
    size += sizeof(dense_qp_dims);
    size += d_memsize_dense_qp(dims);
    return size;
}

dense_qp_in *dense_qp_in_assign(void *config, dense_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    dense_qp_in *qp_in = (dense_qp_in *) c_ptr;
    c_ptr += sizeof(dense_qp_in);

    assert((size_t) c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    d_create_dense_qp(dims, qp_in, c_ptr);
    c_ptr += d_memsize_dense_qp(dims);

    qp_in->dim = (dense_qp_dims *) c_ptr;
    c_ptr += sizeof(dense_qp_dims);

    qp_in->dim->nv = dims->nv;
    qp_in->dim->ne = dims->ne;
    qp_in->dim->nb = dims->nb;
    qp_in->dim->ng = dims->ng;
    qp_in->dim->ns = dims->ns;
    qp_in->dim->nsb = dims->nsb;
    qp_in->dim->nsg = dims->nsg;

    assert((char *) raw_memory + dense_qp_in_calculate_size(config, dims) == c_ptr);

    return qp_in;
}

/************************************************
 * out
 ************************************************/

int dense_qp_out_calculate_size(void *config, dense_qp_dims *dims)
{
    int size = sizeof(dense_qp_out);
    size += d_memsize_dense_qp_sol(dims);
    size += sizeof(dense_qp_info);
    return size;
}

dense_qp_out *dense_qp_out_assign(void *config, dense_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    dense_qp_out *qp_out = (dense_qp_out *) c_ptr;
    c_ptr += sizeof(dense_qp_out);

    assert((size_t) c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    d_create_dense_qp_sol(dims, qp_out, c_ptr);
    c_ptr += d_memsize_dense_qp_sol(dims);

    qp_out->misc = (void *) c_ptr;
    c_ptr += sizeof(dense_qp_info);

    assert((char *) raw_memory + dense_qp_out_calculate_size(config, dims) == c_ptr);

    return qp_out;
}

/************************************************
 * res
 ************************************************/

int dense_qp_res_calculate_size(dense_qp_dims *dims)
{
    int size = sizeof(dense_qp_res);
    size += d_memsize_dense_qp_res(dims);
    return size;
}

dense_qp_res *dense_qp_res_assign(dense_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    dense_qp_res *qp_res = (dense_qp_res *) c_ptr;
    c_ptr += sizeof(dense_qp_res);

    d_create_dense_qp_res(dims, qp_res, c_ptr);
    c_ptr += d_memsize_dense_qp_res(dims);

    assert((char *) raw_memory + dense_qp_res_calculate_size(dims) == c_ptr);

    return qp_res;
}

int dense_qp_res_workspace_calculate_size(dense_qp_dims *dims)
{
    int size = sizeof(dense_qp_res_ws);
    size += d_memsize_dense_qp_res_workspace(dims);
    return size;
}

dense_qp_res_ws *dense_qp_res_workspace_assign(dense_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    dense_qp_res_ws *res_ws = (dense_qp_res_ws *) c_ptr;
    c_ptr += sizeof(dense_qp_res_ws);

    d_create_dense_qp_res_workspace(dims, res_ws, c_ptr);
    c_ptr += d_memsize_dense_qp_res_workspace(dims);

    assert((char *) raw_memory + dense_qp_res_workspace_calculate_size(dims) == c_ptr);

    return res_ws;
}

void dense_qp_compute_t(dense_qp_in *qp_in, dense_qp_out *qp_out)
{
    int nvd = qp_in->dim->nv;
    // int ned = qp_in->dim->ne;
    int ngd = qp_in->dim->ng;
    int nbd = qp_in->dim->nb;
    int nsd = qp_in->dim->ns;

    int *idxb = qp_in->idxb;
    int *idxs = qp_in->idxs;

    // compute slacks for bounds
    blasfeo_dvecex_sp(nbd, 1.0, idxb, qp_out->v, 0, qp_out->t, nbd+ngd);
    blasfeo_daxpby(nbd, 1.0, qp_out->t, nbd+ngd, -1.0, qp_in->d, 0, qp_out->t, 0);
    blasfeo_daxpby(nbd, -1.0, qp_out->t, nbd+ngd, -1.0, qp_in->d, nbd + ngd, qp_out->t, nbd + ngd);

    // compute slacks for general constraints
    blasfeo_dgemv_t(nvd, ngd, 1.0, qp_in->Ct, 0, 0, qp_out->v, 0, -1.0, qp_in->d, nbd, qp_out->t,
                    nbd);
    blasfeo_dgemv_t(nvd, ngd, -1.0, qp_in->Ct, 0, 0, qp_out->v, 0, -1.0, qp_in->d, 2 * nbd + ngd,
                    qp_out->t, 2 * nbd + ngd);

	// soft
	blasfeo_dvecad_sp(nsd, 1.0, qp_out->v, nvd, idxs, qp_out->t, 0);
	blasfeo_dvecad_sp(nsd, 1.0, qp_out->v, nvd+nsd, idxs, qp_out->t, nbd+ngd);
//	blasfeo_daxpy(2*nsd, -1.0, qp_out->v, nvd, qp_out->t, 2*nbd+2*ngd, qp_out->t, 2*nbd+2*ngd);
//	blasfeo_dvecse(2*nsd, 0.0, qp_out->t, 2*nbd+2*ngd);
   	blasfeo_dveccp(2*nsd, qp_out->v, nvd, qp_out->t, 2*nbd+2*ngd);
    blasfeo_daxpy(2*nsd, -1.0, qp_in->d, 2*nbd+2*ngd, qp_out->t, 2*nbd+2*ngd, qp_out->t, 2*nbd+2*ngd);

}

void dense_qp_res_compute(dense_qp_in *qp_in, dense_qp_out *qp_out, dense_qp_res *qp_res,
                          dense_qp_res_ws *res_ws)
{
    dense_qp_info *info = (dense_qp_info *) qp_out->misc;

    if (info->t_computed == 0)
    {
        dense_qp_compute_t(qp_in, qp_out);
        info->t_computed = 1;
    }

    // compute residuals
    d_compute_res_dense_qp(qp_in, qp_out, qp_res, res_ws);
}



void dense_qp_res_compute_nrm_inf(dense_qp_res *qp_res, double res[4])
{
    int nv = qp_res->dim->nv;
    int nb = qp_res->dim->nb;
    int ne = qp_res->dim->ne;
    int ng = qp_res->dim->ng;
    int ns = qp_res->dim->ns;

    blasfeo_dvecnrm_inf(nv + 2 * ns, qp_res->res_g, 0, &res[0]);
    blasfeo_dvecnrm_inf(ne, qp_res->res_b, 0, &res[1]);
    blasfeo_dvecnrm_inf(2 * nb + 2 * ng + 2 * ns, qp_res->res_d, 0, &res[2]);
    blasfeo_dvecnrm_inf(2 * nb + 2 * ng + 2 * ns, qp_res->res_m, 0, &res[3]);
}



void dense_qp_stack_slacks_dims(dense_qp_dims *in, dense_qp_dims *out)
{
    out->nv = in->nv + 2 * in->ns;
    out->ne = in->ne;
    out->nb = in->nb - in->nsb + 2 * in->ns;
    out->ng = in->ns > 0 ? 2 * in->ng + 2 * in->nsb : in->ng;
    out->ns = 0;
    out->nsb = 0;
    out->nsg = 0;
}



void dense_qp_stack_slacks(dense_qp_in *in, dense_qp_in *out)
{
    int nv = in->dim->nv;
    int ne = in->dim->ne;
    int nb = in->dim->nb;
    int ng = in->dim->ng;
    int ns = in->dim->ns;
    int nsb = in->dim->nsb;
    int nsg = in->dim->nsg;
    int *idxs = in->idxs;
    int *idxb = in->idxb;

    int nv2 = out->dim->nv;
    int nb2 = out->dim->nb;
    int ng2 = out->dim->ng;

    assert(nv2 == nv+2*ns && "Dimensions are wrong!");
    assert(nb2 == nb-nsb+2*ns && "Dimensions are wrong!");
    assert(ng2 == 2*ng+2*nsb && "Dimensions are wrong!");

    // set matrice to 0.0
    blasfeo_dgese(nv2, nv2, 0.0, out->Hv, 0, 0);
    blasfeo_dgese(nv2, ng2, 0.0, out->Ct, 0, 0);

    // copy in->Hv to upper left corner of out->Hv, out->Hv = [in->Hv 0; 0 0]
    blasfeo_dgecp(nv, nv, in->Hv, 0, 0, out->Hv, 0, 0);

    // copy in->Z to main diagonal of out->Hv, out->Hv = [in->Hv 0; 0 Z]
    blasfeo_ddiain(2 * ns, 1.0, in->Z, 0, out->Hv, nv, nv);

    // copy in->gz to out->gz
    blasfeo_dveccp(nv + 2 * ns, in->gz, 0, out->gz, 0);

    if (ne > 0)
    {
        // copy in->A to out->A
        blasfeo_dgecp(ne, nv, in->A, 0, 0, out->A, 0, 0);

        // copy in->b to out->b
        blasfeo_dveccp(ne, in->b, 0, out->b, 0);
    }

    // copy in->Ct to out->Ct, out->Ct = [in->Ct 0 0 0; 0 0 0 0]
    blasfeo_dgecp(nv, ng, in->Ct, 0, 0, out->Ct, 0, 0);

    if (ns > 0)
    {
        // copy -in->Ct to out->Ct, out->Ct = [in->Ct -in->Ct 0 0; 0 0 0 0]
        blasfeo_dgecpsc(nv, ng, -1.0, in->Ct, 0, 0, out->Ct, 0, ng);

        // set flags for non-softened box constraints
        // use out->m temporarily for this
        for (int ii = 0; ii < nb; ii++)
        {
            BLASFEO_DVECEL(out->m, ii) = 1.0;
        }

        int k_s = 0, col_b = 2*ng;
        for (int ii = 0; ii < ns; ii++)
        {
            int js = idxs[ii];

            if (js < nb)
            {
                // index of a soft box constraint
                int jv = idxb[js];

                // softened box constraint, set its flag to -1
                BLASFEO_DVECEL(out->m, js) = -1.0;

                // insert softened box constraint into out->Ct, x_i - su_i <= ub_i
                BLASFEO_DMATEL(out->Ct, jv, col_b) = 1.0;
                BLASFEO_DMATEL(out->Ct, nv+ns+k_s, col_b) = -1.0;
                BLASFEO_DVECEL(out->d, 2*nb2+ng2+col_b) = -BLASFEO_DVECEL(in->d, nb+ng+js);

                // insert softened box constraint into out->Ct, -x_i - sl_i <= -lb_i
                BLASFEO_DMATEL(out->Ct, jv, col_b+nsb) = -1.0;
                BLASFEO_DMATEL(out->Ct, nv+k_s, col_b+nsb) = -1.0;
                BLASFEO_DVECEL(out->d, 2*nb2+ng2+col_b+nsb) = -BLASFEO_DVECEL(in->d, js);

                col_b++;
            }
            else
            {
                // index of a soft general constraint
                int col_g = js - nb;

                // C_i x - su_i <= ug_i
                BLASFEO_DMATEL(out->Ct, nv + ns + k_s, col_g) = -1.0;

                // -C_i x - sl_i <= -lg_i
                BLASFEO_DMATEL(out->Ct, nv + k_s, col_g+ng) = -1.0;
            }

            k_s++;

            // slack variables have box constraints
            out->idxb[nb-nsb+ii] = ii + nv;
            out->idxb[nb-nsb+ns+ii] = ii + nv + ns;
        }

        int k_nsb = 0;
        for (int ii = 0; ii < nb; ii++)
        {
            if (BLASFEO_DVECEL(out->m, ii) > 0)
            {
                // copy nonsoftened box constraint bounds to out->d
                BLASFEO_DVECEL(out->d, k_nsb) = BLASFEO_DVECEL(in->d, ii);
                BLASFEO_DVECEL(out->d, nb2+ng2+k_nsb) = -BLASFEO_DVECEL(in->d, nb+ng+ii);
                out->idxb[k_nsb] = ii;
                k_nsb++;
            }
        }

        assert(k_nsb == nb-nsb && "Dimensions are wrong!");

        // copy ls and us to out->lb, jump over nonsoftened box constraints
        blasfeo_dveccp(2*ns, in->d, 2*nb+2*ng, out->d, k_nsb);

        // for slack variables out->ub is +INFTY
        blasfeo_dvecse(2*ns, 1.0e6, out->d, nb2+ng2+k_nsb);

        // set out->lg to -INFTY
        blasfeo_dvecse(ng2, -1.0e6, out->d, nb2);

        // copy in->ug and -in->lg to out->ug
        blasfeo_dveccpsc(ng, -1.0, in->d, 2*nb+ng, out->d, 2*nb2+ng2);
        blasfeo_dveccpsc(ng, -1.0, in->d, nb, out->d, 2*nb2+ng2+ng);

        // flip signs for ub and ug
        blasfeo_dvecsc(nb2+ng2, -1.0, out->d, nb2+ng2);

        // set out->m to 0.0
        blasfeo_dvecse(2*nb2+2*ng2, 0.0, out->m, 0);
    }
    else
    {
        blasfeo_dveccp(2*nb+2*ng, in->d, 0, out->d, 0);
        blasfeo_dveccp(2*nb+2*ng, in->m, 0, out->m, 0);
    }
}
