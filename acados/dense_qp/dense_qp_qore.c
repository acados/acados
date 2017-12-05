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
#include <math.h>
#include <string.h>
// blasfeo
#include "blasfeo_target.h"
#include "blasfeo_common.h"
#include "blasfeo_d_aux.h"
#include "blasfeo_d_aux_ext_dep.h"
// acados
#include "acados/dense_qp/dense_qp_qore.h"
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/utils/mem.h"
#include "acados/utils/timing.h"


int dense_qp_qore_calculate_args_size(dense_qp_dims *dims)
{
    int size = 0;
    size += sizeof(dense_qp_qore_args);

    return size;
}



void *dense_qp_qore_assign_args(dense_qp_dims *dims, void *raw_memory)
{
    dense_qp_qore_args *args;

    char *c_ptr = (char *) raw_memory;

    args = (dense_qp_qore_args *) c_ptr;
    c_ptr += sizeof(dense_qp_qore_args);

    assert((char*)raw_memory + dense_qp_qore_calculate_args_size(dims) == c_ptr);

    return (void *)args;
}



void dense_qp_qore_initialize_default_args(void *args_)
{
    dense_qp_qore_args *args = (dense_qp_qore_args *)args_;

    args->prtfreq = -1;
    args->warm_start = 0;
    args->warm_strategy = 0;
    args->nsmax = 400;
    args->hot_start = 0;
}



int dense_qp_qore_calculate_memory_size(dense_qp_dims *dims, void *args_)
{
    dense_qp_qore_args *args = (dense_qp_qore_args *) args_;

    int nvd = dims->nv;
    int ned = dims->ne;
    int ngd = dims->ng;
    int nbd = dims->nb;
    int nsmax = (2*nvd >= args->nsmax) ? args->nsmax : 2*nvd;

    // size in bytes
    int size = sizeof(dense_qp_qore_memory);

    size += 1 * nvd * nvd * sizeof(double);        // H
    size += 1 * nvd * ned * sizeof(double);        // A
    size += 2 * nvd * ngd * sizeof(double);        // C, Ct
    size += 3 * nvd * sizeof(double);              // g d_lb d_ub
    size += 1 * ned * sizeof(double);              // b
    size += 2 * nbd * sizeof(double);              // d_lb0 d_ub0
    size += 2 * ngd * sizeof(double);              // d_lg d_ug
    size += 1 * nbd * sizeof(int);                 // idxb
    size += 2 * (nvd + ngd) * sizeof(double);      // lb, ub
    size += (nvd + ngd) * sizeof(double);          // prim_sol
    size += (nvd + ngd) * sizeof(double);          // dual_sol
    size += QPDenseSize(nvd,ngd,nsmax);

    make_int_multiple_of(8, &size);

    return size;
}



void *dense_qp_qore_assign_memory(dense_qp_dims *dims, void *args_, void *raw_memory)
{
    dense_qp_qore_memory *mem;
    dense_qp_qore_args *args = (dense_qp_qore_args *) args_;

    int nvd = dims->nv;
    int ned = dims->ne;
    int ngd = dims->ng;
    int nbd = dims->nb;
    int nsmax = (2*nvd >= args->nsmax) ? args->nsmax : 2*nvd;

    // char pointer
    char *c_ptr = (char *)raw_memory;

    mem = (dense_qp_qore_memory *) c_ptr;
    c_ptr += sizeof(dense_qp_qore_memory);

    assert((size_t)c_ptr % 8 == 0 && "double not 8-byte aligned!");

    assign_double(nvd*nvd, &mem->H, &c_ptr);
    assign_double(nvd*ned, &mem->A, &c_ptr);
    assign_double(nvd*ngd, &mem->C, &c_ptr);
    assign_double(nvd*ngd, &mem->Ct, &c_ptr);
    assign_double(nvd, &mem->g, &c_ptr);
    assign_double(ned, &mem->b, &c_ptr);
    assign_double(nbd, &mem->d_lb0, &c_ptr);
    assign_double(nbd, &mem->d_ub0, &c_ptr);
    assign_double(nvd, &mem->d_lb, &c_ptr);
    assign_double(nvd, &mem->d_ub, &c_ptr);
    assign_double(ngd, &mem->d_lg, &c_ptr);
    assign_double(ngd, &mem->d_ug, &c_ptr);
    assign_double(nvd+ngd, &mem->lb, &c_ptr);
    assign_double(nvd+ngd, &mem->ub, &c_ptr);
    assign_double(nvd+ngd, &mem->prim_sol, &c_ptr);
    assign_double(nvd+ngd, &mem->dual_sol, &c_ptr);

    assert((size_t)c_ptr % 8 == 0 && "double not 8-byte aligned!");

    mem->QP = (QoreProblemDense *) c_ptr;
    QPDenseCreate(&mem->QP, nvd, ngd, nsmax, c_ptr);
    c_ptr += QPDenseSize(nvd,ngd,nsmax);

    // int stuff
    assign_int(nbd, &mem->idxb, &c_ptr);

    assert((char *)raw_memory + dense_qp_qore_calculate_memory_size(dims, args_) >= c_ptr);

    return mem;
}



int dense_qp_qore_calculate_workspace_size(dense_qp_dims *dims, void *args_)
{
    return 0;
}



int dense_qp_qore(dense_qp_in *qp_in, dense_qp_out *qp_out, void *args_, void *memory_, void *work_)
{
    dense_qp_info *info = (dense_qp_info *) qp_out->misc;
    acados_timer tot_timer, qp_timer, interface_timer;

    acados_tic(&tot_timer);
    acados_tic(&interface_timer);

    // cast structures
    dense_qp_qore_args *args = (dense_qp_qore_args *)args_;
    dense_qp_qore_memory *memory = (dense_qp_qore_memory *)memory_;

    // initialize return code
    int acados_status = ACADOS_SUCCESS;

    // extract qpoases data
    double *H = memory->H;
    double *A = memory->A;
    double *C = memory->C;
    double *Ct = memory->Ct;
    double *g = memory->g;
    double *b = memory->b;
    double *d_lb0 = memory->d_lb0;
    double *d_ub0 = memory->d_ub0;
    double *d_lb = memory->d_lb;
    double *d_ub = memory->d_ub;
    double *d_lg = memory->d_lg;
    double *d_ug = memory->d_ug;
    double *lb = memory->lb;
    double *ub = memory->ub;
    int *idxb = memory->idxb;
    double *prim_sol = memory->prim_sol;
    double *dual_sol = memory->dual_sol;
    QoreProblemDense *QP = memory->QP;
    int num_iter;

    // extract dense qp size
    int nvd = qp_in->dim->nv;
    int ned = qp_in->dim->ne;
    int ngd = qp_in->dim->ng;
    int nbd = qp_in->dim->nb;

    assert(ned == 0 && "ned != 0 not supported yet");

    // fill in the upper triangular of H in dense_qp
    dtrtr_l_libstr(nvd, qp_in->Hv, 0, 0, qp_in->Hv, 0, 0);

    // dense qp row-major
    d_cvt_dense_qp_to_colmaj(qp_in, H, g, A, b, idxb, d_lb0, d_ub0, C, d_lg, d_ug,
        NULL, NULL, NULL, NULL, NULL);

    // reorder bounds
    for (int ii = 0; ii < nvd; ii++) {
        d_lb[ii] = -INFINITY;
        d_ub[ii] = +INFINITY;
    }
    for (int ii = 0; ii < nbd; ii++) {
        d_lb[idxb[ii]] = d_lb0[ii];
        d_ub[idxb[ii]] = d_ub0[ii];
    }

    memcpy(lb, d_lb, nvd*sizeof(double));
    memcpy(lb+nvd, d_lg, ngd*sizeof(double));

    memcpy(ub, d_ub, nvd*sizeof(double));
    memcpy(ub+nvd, d_ug, ngd*sizeof(double));

    /* transpose C as expected by QORE */
	int i, j;
	for (j=0; j<nvd; j++) {
		for (i=0; i<ngd; i++) {
			Ct[j+i*nvd] = C[i+j*ngd];
		}
    }

    // solve dense qp
    int prtfreq = args->prtfreq;
    int warm_start = args->warm_start;
    int warm_strategy = args->warm_strategy;
    int hot_start = args->hot_start;
    int return_flag = 0;

    info->interface_time = acados_toc(&interface_timer);
    acados_tic(&qp_timer);

    if (warm_start)
    {
        QPDenseSetInt(QP, "warmstrategy", warm_strategy);
        QPDenseUpdateMatrices(QP, nvd, ngd, Ct, H);
    }
    else if (!hot_start)
    {
        QPDenseSetData(QP, nvd, ngd, Ct, H);
    }

    QPDenseSetInt(QP, "prtfreq", prtfreq);
    return_flag = QPDenseOptimize( QP, lb, ub, g, 0, 0 );

    QPDenseGetDblVector( QP, "primalsol", prim_sol );
    QPDenseGetDblVector( QP, "dualsol", dual_sol );
    QPDenseGetInt(QP, "itercount", &num_iter);
    memory->num_iter = num_iter;

#if 0
    d_print_mat(1, nvd, prim_sol, 1);
    exit(1);
#endif

    info->solve_QP_time = acados_toc(&qp_timer);
    acados_tic(&interface_timer);

    // copy prim_sol and dual_sol to qpd_sol
    d_cvt_vec2strvec(nvd, prim_sol, qp_out->v, 0);
    for (int ii = 0; ii < 2*nbd+2*ngd; ii++)
        qp_out->lam->pa[ii] = 0.0;
    for (int ii = 0; ii < nbd; ii++) {
        if (dual_sol[idxb[ii]] >= 0.0)
            qp_out->lam->pa[ii] = dual_sol[idxb[ii]];
        else
            qp_out->lam->pa[nbd+ngd+ii] = - dual_sol[idxb[ii]];
    }
    for (int ii = 0; ii < ngd; ii++) {
        if (dual_sol[nvd+ii] >= 0.0)
            qp_out->lam->pa[nbd+ii] =   dual_sol[nvd+ii];
        else
            qp_out->lam->pa[2*nbd+ngd+ii] = - dual_sol[nvd+ii];
    }

    info->interface_time += acados_toc(&interface_timer);
    info->total_time = acados_toc(&tot_timer);

    // return
    // TODO(bnovoselnik): cast qore return to acados return
    acados_status = return_flag;
    return acados_status;
}
