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
#include <math.h>
#include <string.h>
// blasfeo
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
// acados
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/dense_qp/dense_qp_qore.h"
#include "acados/utils/mem.h"
#include "acados/utils/timing.h"



/************************************************
 * opts
 ************************************************/

int dense_qp_qore_opts_calculate_size(void *config_, dense_qp_dims *dims)
{
    int size = 0;
    size += sizeof(dense_qp_qore_opts);

    return size;
}



void *dense_qp_qore_opts_assign(void *config_, dense_qp_dims *dims, void *raw_memory)
{
    dense_qp_qore_opts *opts;

    char *c_ptr = (char *) raw_memory;

    opts = (dense_qp_qore_opts *) c_ptr;
    c_ptr += sizeof(dense_qp_qore_opts);

    assert((char *) raw_memory + dense_qp_qore_opts_calculate_size(config_, dims) >= c_ptr);

    return (void *) opts;
}



void dense_qp_qore_opts_initialize_default(void *config_, dense_qp_dims *dims, void *opts_)
{
    dense_qp_qore_opts *opts = (dense_qp_qore_opts *) opts_;

    opts->print_freq = -1;
    opts->warm_start = 0;
    opts->warm_strategy = 0;
    opts->nsmax = 400;
    opts->hot_start = 0;
    opts->max_iter = 1000;
    opts->compute_t = 1;

    return;
}



void dense_qp_qore_opts_update(void *config_, dense_qp_dims *dims, void *opts_)
{
    //    dense_qp_qore_opts *opts = (dense_qp_qore_opts *)opts_;

    return;
}



void dense_qp_qore_opts_set(void *config_, void *opts_, const char *field, void *value)
{
    dense_qp_qore_opts *opts = opts_;

    if (!strcmp(field, "tol_stat"))
    {
		// TODO set solver exit tolerance
    }
    else if (!strcmp(field, "tol_eq"))
    {
		// TODO set solver exit tolerance
    }
    else if (!strcmp(field, "tol_ineq"))
    {
		// TODO set solver exit tolerance
    }
    else if (!strcmp(field, "tol_comp"))
    {
		// TODO set solver exit tolerance
    }
    else if (!strcmp(field, "warm_start"))
    {
		// TODO set solver warm start
    }
	else
	{
		printf("\nerror: dense_qp_qore_opts_set: wrong field: %s\n", field);
		exit(1);
	}

	return;
}



/************************************************
 * memory
 ************************************************/

int dense_qp_qore_memory_calculate_size(void *config_, dense_qp_dims *dims, void *opts_)
{
    dense_qp_qore_opts *opts = (dense_qp_qore_opts *) opts_;
    dense_qp_dims dims_stacked;

    int nv = dims->nv;
    int ne = dims->ne;
    int ng = dims->ng;
    int nb = dims->nb;
    int ns = dims->ns;
    int nsb = dims->nsb;
    // int nsg = dims->nsg;
    int nsmax = (2 * nv >= opts->nsmax) ? opts->nsmax : 2 * nv;

    int nv2 = nv + 2*ns;
    int ng2 = (ns > 0) ? ng + nsb : ng;
    int nb2 = nb - nsb + 2*ns;

    // size in bytes
    int size = sizeof(dense_qp_qore_memory);

    size += 1 * nv * nv * sizeof(double);      // H
    size += 1 * nv2 * nv2 * sizeof(double);    // HH
    size += 1 * nv2 * ne * sizeof(double);     // A
    size += 2 * nv * ng * sizeof(double);      // C, Ct
    size += 2 * nv2 * ng2 * sizeof(double);    // CC, CCt
    size += 1 * nv * sizeof(double);           // g
    size += 3 * nv2 * sizeof(double);          // gg d_lb d_ub
    size += 1 * ne * sizeof(double);           // b
    size += 2 * nb2 * sizeof(double);          // d_lb0 d_ub0
    size += 2 * ng2 * sizeof(double);          // d_lg d_ug
    size += 1 * nb * sizeof(int);              // idxb
    size += 1 * nb2 * sizeof(int);             // idxb_stacked
    size += 1 * ns * sizeof(int);
    size += 2 * (nv2 + ng2) * sizeof(double);  // lb, ub
    size += 1 * (nv2 + ng2) * sizeof(double);  // prim_sol
    size += 1 * (nv2 + ng2) * sizeof(double);  // dual_sol
    size += 6 * ns * sizeof(double);           // Zl, Zu, zl, zu, d_ls, d_us
    size += QPDenseSize(nv2, ng2, nsmax);

    if (ns > 0)
    {
        dense_qp_stack_slacks_dims(dims, &dims_stacked);
        size += dense_qp_in_calculate_size(config_, &dims_stacked);
    }

    make_int_multiple_of(8, &size);

    return size;
}



void *dense_qp_qore_memory_assign(void *config_, dense_qp_dims *dims, void *opts_, void *raw_memory)
{
    dense_qp_qore_memory *mem;
    dense_qp_qore_opts *opts = (dense_qp_qore_opts *) opts_;
    dense_qp_dims dims_stacked;

    int nv = dims->nv;
    int ne = dims->ne;
    int ng = dims->ng;
    int nb = dims->nb;
    int ns = dims->ns;
    int nsb = dims->nsb;
    // int nsg = dims->nsg;
    int nsmax = (2 * nv >= opts->nsmax) ? opts->nsmax : 2 * nv;

    int nv2 = nv + 2*ns;
    int ng2 = (ns > 0) ? ng + nsb : ng;
    int nb2 = nb - nsb + 2*ns;

    // char pointer
    char *c_ptr = (char *) raw_memory;

    mem = (dense_qp_qore_memory *) c_ptr;
    c_ptr += sizeof(dense_qp_qore_memory);

    assert((size_t) c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    if (ns > 0)
    {
        dense_qp_stack_slacks_dims(dims, &dims_stacked);
        mem->qp_stacked = dense_qp_in_assign(config_, &dims_stacked, c_ptr);
        c_ptr += dense_qp_in_calculate_size(config_, &dims_stacked);
    }
    else
    {
        mem->qp_stacked = NULL;
    }

    assert((size_t) c_ptr % 8 == 0 && "double not 8-byte aligned!");

    assign_and_advance_double(nv * nv, &mem->H, &c_ptr);
    assign_and_advance_double(nv2 * nv2, &mem->HH, &c_ptr);
    assign_and_advance_double(nv2 * ne, &mem->A, &c_ptr);
    assign_and_advance_double(nv * ng, &mem->C, &c_ptr);
    assign_and_advance_double(nv * ng, &mem->Ct, &c_ptr);
    assign_and_advance_double(nv2 * ng2, &mem->CC, &c_ptr);
    assign_and_advance_double(nv2 * ng2, &mem->CCt, &c_ptr);
    assign_and_advance_double(nv, &mem->g, &c_ptr);
    assign_and_advance_double(nv2, &mem->gg, &c_ptr);
    assign_and_advance_double(ne, &mem->b, &c_ptr);
    assign_and_advance_double(nb2, &mem->d_lb0, &c_ptr);
    assign_and_advance_double(nb2, &mem->d_ub0, &c_ptr);
    assign_and_advance_double(nv2, &mem->d_lb, &c_ptr);
    assign_and_advance_double(nv2, &mem->d_ub, &c_ptr);
    assign_and_advance_double(ng2, &mem->d_lg, &c_ptr);
    assign_and_advance_double(ng2, &mem->d_ug, &c_ptr);
    assign_and_advance_double(ns, &mem->Zl, &c_ptr);
    assign_and_advance_double(ns, &mem->Zu, &c_ptr);
    assign_and_advance_double(ns, &mem->zl, &c_ptr);
    assign_and_advance_double(ns, &mem->zu, &c_ptr);
    assign_and_advance_double(ns, &mem->d_ls, &c_ptr);
    assign_and_advance_double(ns, &mem->d_us, &c_ptr);
    assign_and_advance_double(nv2 + ng2, &mem->lb, &c_ptr);
    assign_and_advance_double(nv2 + ng2, &mem->ub, &c_ptr);
    assign_and_advance_double(nv2 + ng2, &mem->prim_sol, &c_ptr);
    assign_and_advance_double(nv2 + ng2, &mem->dual_sol, &c_ptr);

    assert((size_t) c_ptr % 8 == 0 && "double not 8-byte aligned!");

    mem->QP = (QoreProblemDense *) c_ptr;
    QPDenseCreate(&mem->QP, nv2, ng2, nsmax, c_ptr);
    c_ptr += QPDenseSize(nv2, ng2, nsmax);

    // int stuff
    assign_and_advance_int(nb, &mem->idxb, &c_ptr);
    assign_and_advance_int(nb2, &mem->idxb_stacked, &c_ptr);
    assign_and_advance_int(ns, &mem->idxs, &c_ptr);

    assert((char *) raw_memory + dense_qp_qore_memory_calculate_size(config_, dims, opts_) >=
           c_ptr);

    return mem;
}



/************************************************
 * workspace
 ************************************************/

int dense_qp_qore_workspace_calculate_size(void *config_, dense_qp_dims *dims, void *opts_)
{
    return 0;
}



/************************************************
 * functions
 ************************************************/

int dense_qp_qore(void *config_, dense_qp_in *qp_in, dense_qp_out *qp_out, void *opts_,
                  void *memory_, void *work_)
{
    dense_qp_info *info = (dense_qp_info *) qp_out->misc;
    acados_timer tot_timer, qp_timer, interface_timer;

    acados_tic(&tot_timer);
    acados_tic(&interface_timer);
    info->t_computed = 0;

    // cast structures
    dense_qp_qore_opts *opts = (dense_qp_qore_opts *) opts_;
    dense_qp_qore_memory *memory = (dense_qp_qore_memory *) memory_;

    // extract qpoases data
    double *H = memory->H;
    double *HH = memory->HH;
    double *A = memory->A;
    double *C = memory->C;
    double *CC = memory->CC;
    double *Ct = memory->Ct;
    double *CCt = memory->CCt;
    // double *g = memory->g;
    double *gg = memory->gg;
    double *b = memory->b;
    double *d_lb0 = memory->d_lb0;
    double *d_ub0 = memory->d_ub0;
    double *d_lb = memory->d_lb;
    double *d_ub = memory->d_ub;
    double *d_lg = memory->d_lg;
    double *d_ug = memory->d_ug;
    double *lb = memory->lb;
    double *ub = memory->ub;
    double *Zl = memory->Zl;
    double *Zu = memory->Zu;
    double *zl = memory->zl;
    double *zu = memory->zu;
    double *d_ls = memory->d_ls;
    double *d_us = memory->d_us;
    int *idxb = memory->idxb;
    int *idxb_stacked = memory->idxb_stacked;
    int *idxs = memory->idxs;
    QoreProblemDense *QP = memory->QP;
    double *prim_sol = memory->prim_sol;
    double *dual_sol = memory->dual_sol;
    dense_qp_in *qp_stacked = memory->qp_stacked;

    // extract dense qp size
    int nv = qp_in->dim->nv;
    int ng = qp_in->dim->ng;
    int nb = qp_in->dim->nb;
    int ns = qp_in->dim->ns;
    int nsb = qp_in->dim->nsb;
    // int nsg = qp_in->dim->nsg;

    int nv2 = nv + 2*ns;
    int ng2 = (ns > 0) ? ng + nsb : ng;
    int nb2 = nb - nsb + 2 * ns;

    // fill in the upper triangular of H in dense_qp
    blasfeo_dtrtr_l(nv, qp_in->Hv, 0, 0, qp_in->Hv, 0, 0);

    // extract data from qp_in in col-major
    d_cvt_dense_qp_to_colmaj(qp_in, H, gg, A, b, idxb, d_lb0, d_ub0, C, d_lg, d_ug,
                             Zl, Zu, zl, zu, idxs, d_ls, d_us);

    // reorder bounds
    for (int ii = 0; ii < nv2; ii++)
    {
        d_lb[ii] = -INFINITY;
        d_ub[ii] = +INFINITY;
    }

    if (ns > 0)
    {
        dense_qp_stack_slacks(qp_in, qp_stacked);
        d_cvt_dense_qp_to_colmaj(qp_stacked, HH, gg, A, b, idxb_stacked, d_lb0, d_ub0, CC, d_lg,
            d_ug, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

        for (int ii = 0; ii < nb2; ii++)
        {
            d_lb[idxb_stacked[ii]] = d_lb0[ii];
            d_ub[idxb_stacked[ii]] = d_ub0[ii];
        }

        // transpose CC as expected by QORE
        for (int j = 0; j < nv2; j++)
        {
            for (int i = 0; i < ng2; i++)
            {
                CCt[j + i * nv2] = CC[i + j * ng2];
            }
        }
    }
    else
    {
        for (int ii = 0; ii < nb; ii++)
        {
            d_lb[idxb[ii]] = d_lb0[ii];
            d_ub[idxb[ii]] = d_ub0[ii];
        }

        // transpose C as expected by QORE
        for (int j = 0; j < nv; j++)
        {
            for (int i = 0; i < ng; i++)
            {
                Ct[j + i * nv] = C[i + j * ng];
            }
        }
    }

    memcpy(lb, d_lb, nv2 * sizeof(double));
    memcpy(lb + nv2, d_lg, ng2 * sizeof(double));

    memcpy(ub, d_ub, nv2 * sizeof(double));
    memcpy(ub + nv2, d_ug, ng2 * sizeof(double));


    info->interface_time = acados_toc(&interface_timer);

    // solve dense qp
    acados_tic(&qp_timer);

    if (opts->warm_start)
    {
        QPDenseSetInt(QP, "warmstrategy", opts->warm_strategy);
        (ns > 0) ? QPDenseUpdateMatrices(QP, nv2, ng2, CCt, HH) :
                   QPDenseUpdateMatrices(QP, nv, ng, Ct, H);
    }
    else if (!opts->hot_start)
    {
        (ns > 0) ? QPDenseSetData(QP, nv2, ng2, CCt, HH) :
                   QPDenseSetData(QP, nv, ng, Ct, H);
    }

    QPDenseSetInt(QP, "maxiter", opts->max_iter);
    QPDenseSetInt(QP, "prtfreq", opts->print_freq);
    QPDenseOptimize(QP, lb, ub, gg, 0, 0);
    int qore_status;
    QPDenseGetInt(QP, "status", &qore_status);

    QPDenseGetDblVector(QP, "primalsol", prim_sol);
    QPDenseGetDblVector(QP, "dualsol", dual_sol);
    int num_iter;
    QPDenseGetInt(QP, "itercount", &num_iter);
    memory->num_iter = num_iter;

    info->solve_QP_time = acados_toc(&qp_timer);
    acados_tic(&interface_timer);

    // copy prim_sol and dual_sol to qpd_sol
    blasfeo_pack_dvec(nv2, prim_sol, qp_out->v, 0);
    for (int ii = 0; ii < 2 * nb + 2 * ng + 2 * ns; ii++) qp_out->lam->pa[ii] = 0.0;
    for (int ii = 0; ii < nb; ii++)
    {
        if (dual_sol[idxb[ii]] >= 0.0)
            qp_out->lam->pa[ii] = dual_sol[idxb[ii]];
        else
            qp_out->lam->pa[nb + ng + ii] = -dual_sol[idxb[ii]];
    }

    for (int ii = 0; ii < ng; ii++)
    {
        if (dual_sol[nv2 + ii] >= 0.0)
            qp_out->lam->pa[nb + ii] = dual_sol[nv2 + ii];
        else
            qp_out->lam->pa[2 * nb + ng + ii] = -dual_sol[nv2 + ii];
    }

    int k = 0;
    for (int ii = 0; ii < ns; ii++)
    {
        int js = idxs[ii];

        double offset_l = 0.0;
        double offset_u = 0.0;

        if (js < nb)
        {
            if (dual_sol[nv2 + ng + k] <= 0.0)  // softened upper box constraints
            {
                qp_out->lam->pa[nb + ng + js] = -dual_sol[nv2 + ng + k];
                offset_u = -dual_sol[nv2 + ng + k];
            }
            else  // softened lower box constraints
            {
                qp_out->lam->pa[js] = dual_sol[nv2 + ng + k];
                offset_l = dual_sol[nv2 + ng + k];
            }

            k++;
        }
        else
        {
            offset_l = qp_out->lam->pa[nb+js-nb];
            offset_u = qp_out->lam->pa[2*nb+ng+js-nb];
        }

        // dual variables for sl >= d_ls
        if (dual_sol[nv + ii] >= 0)
            qp_out->lam->pa[2*nb + 2*ng + ii] = dual_sol[nv + ii] - offset_u;

        // dual variables for su >= d_us
        if (dual_sol[nv + ns + ii] >= 0)
            qp_out->lam->pa[2*nb + 2*ng + ns + ii] = dual_sol[nv + ns + ii] - offset_l;
    }

    info->interface_time += acados_toc(&interface_timer);
    info->total_time = acados_toc(&tot_timer);
    info->num_iter = num_iter;

    // compute slacks
    if (opts->compute_t)
    {
        dense_qp_compute_t(qp_in, qp_out);
        info->t_computed = 1;
    }

    int acados_status = qore_status;
    if (qore_status == QPSOLVER_DENSE_OPTIMAL) acados_status = ACADOS_SUCCESS;
    if (qore_status == QPSOLVER_DENSE_ITER_LIMIT) acados_status = ACADOS_MAXITER;
    return acados_status;
}



void dense_qp_qore_config_initialize_default(void *config_)
{
    qp_solver_config *config = config_;

    config->opts_calculate_size = (int (*)(void *, void *)) & dense_qp_qore_opts_calculate_size;
    config->opts_assign = (void *(*) (void *, void *, void *) ) & dense_qp_qore_opts_assign;
    config->opts_initialize_default =
        (void (*)(void *, void *, void *)) & dense_qp_qore_opts_initialize_default;
    config->opts_update = (void (*)(void *, void *, void *)) & dense_qp_qore_opts_update;
    config->opts_set = &dense_qp_qore_opts_set;
    config->memory_calculate_size =
        (int (*)(void *, void *, void *)) & dense_qp_qore_memory_calculate_size;
    config->memory_assign =
        (void *(*) (void *, void *, void *, void *) ) & dense_qp_qore_memory_assign;
    config->workspace_calculate_size =
        (int (*)(void *, void *, void *)) & dense_qp_qore_workspace_calculate_size;
    config->evaluate = (int (*)(void *, void *, void *, void *, void *, void *)) & dense_qp_qore;

    return;
}
