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
#include "blasfeo_target.h"
#include "blasfeo_common.h"
#include "blasfeo_d_aux.h"
// qpoases
#include "qpOASES_e.h"
// acados
#include "acados/dense_qp/dense_qp_qpoases.h"
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/utils/mem.h"
#include "acados/utils/timing.h"


int dense_qp_qpoases_calculate_args_size(dense_qp_dims *dims)
{
    int size = 0;
    size += sizeof(dense_qp_qpoases_args);

    return size;
}



void *dense_qp_qpoases_assign_args(dense_qp_dims *dims, void *raw_memory)
{
    dense_qp_qpoases_args *args;

    char *c_ptr = (char *) raw_memory;

    args = (dense_qp_qpoases_args *) c_ptr;
    c_ptr += sizeof(dense_qp_qpoases_args);

    assert((char*)raw_memory + dense_qp_qpoases_calculate_args_size(dims) == c_ptr);

    return (void *)args;
}



void dense_qp_qpoases_initialize_default_args(void *args_)
{
    dense_qp_qpoases_args *args = (dense_qp_qpoases_args *)args_;

    args->max_cputime = 1000.0;
    args->warm_start = 0;
    args->max_nwsr = 1000;
}



int dense_qp_qpoases_calculate_memory_size(dense_qp_dims *dims, void *args_)
{
    dense_qp_qpoases_args *args = (dense_qp_qpoases_args *) args_;

    int nvd = dims->nv;
    int ned = dims->ne;
    int ngd = dims->ng;
    int nbd = dims->nb;

    // size in bytes
    int size = sizeof(dense_qp_qpoases_memory);

    size += 1 * nvd * nvd * sizeof(double);        // H
    size += 1 * nvd * ned * sizeof(double);        // A
    size += 1 * nvd * ngd * sizeof(double);        // C
    size += 3 * nvd * sizeof(double);              // g d_lb d_ub
    size += 1 * ned * sizeof(double);              // b
    size += 2 * nbd * sizeof(double);              // d_lb0 d_ub0
    size += 2 * ngd * sizeof(double);              // d_lg d_ug
    size += 1 * nbd * sizeof(int);                 // idxb
    size += 1 * nvd * sizeof(double);              // prim_sol
    size += (nvd+ngd) * sizeof(double);  // dual_sol

    if (ngd > 0)  // QProblem
        size += QProblem_calculateMemorySize(nvd, ngd);
    else  // QProblemB
        size += QProblemB_calculateMemorySize(nvd);

    make_int_multiple_of(8, &size);

    return size;
}



void *dense_qp_qpoases_assign_memory(dense_qp_dims *dims, void *args_, void *raw_memory)
{
    dense_qp_qpoases_memory *mem;
    dense_qp_qpoases_args *args = (dense_qp_qpoases_args *) args_;

    int nvd = dims->nv;
    int ned = dims->ne;
    int ngd = dims->ng;
    int nbd = dims->nb;

    // char pointer
    char *c_ptr = (char *)raw_memory;

    mem = (dense_qp_qpoases_memory *) c_ptr;
    c_ptr += sizeof(dense_qp_qpoases_memory);

    assert((size_t)c_ptr % 8 == 0 && "double not 8-byte aligned!");

    assign_double(nvd*nvd, &mem->H, &c_ptr);
    assign_double(nvd*ned, &mem->A, &c_ptr);
    assign_double(nvd*ngd, &mem->C, &c_ptr);
    assign_double(nvd, &mem->g, &c_ptr);
    assign_double(ned, &mem->b, &c_ptr);
    assign_double(nbd, &mem->d_lb0, &c_ptr);
    assign_double(nbd, &mem->d_ub0, &c_ptr);
    assign_double(nvd, &mem->d_lb, &c_ptr);
    assign_double(nvd, &mem->d_ub, &c_ptr);
    assign_double(ngd, &mem->d_lg, &c_ptr);
    assign_double(ngd, &mem->d_ug, &c_ptr);
    assign_double(nvd, &mem->prim_sol, &c_ptr);
    assign_double(nvd+ngd, &mem->dual_sol, &c_ptr);

    // TODO(dimitris): update assign syntax in qpOASES
    assert((size_t)c_ptr % 8 == 0 && "double not 8-byte aligned!");

    if (ngd > 0) {  // QProblem
        QProblem_assignMemory(nvd, ngd, (QProblem **) &(mem->QP), c_ptr);
        c_ptr += QProblem_calculateMemorySize(nvd, ngd);
    } else {  // QProblemB
        QProblemB_assignMemory(nvd, (QProblemB **) &(mem->QPB), c_ptr);
        c_ptr += QProblemB_calculateMemorySize(nvd);
    }

    assign_int(nbd, &mem->idxb, &c_ptr);

    assert((char *)raw_memory + dense_qp_qpoases_calculate_memory_size(dims, args_) >= c_ptr);

    return mem;
}



int dense_qp_qpoases_calculate_workspace_size(dense_qp_dims *dims, void *args_)
{
    return 0;
}



int dense_qp_qpoases(dense_qp_in *qp_in, dense_qp_out *qp_out, void *args_, void *memory_, void *work_)
{
    dense_qp_info *info = (dense_qp_info *) qp_out->misc;
    acados_timer tot_timer, qp_timer, interface_timer;

    acados_tic(&tot_timer);
    acados_tic(&interface_timer);

    // cast structures
    dense_qp_qpoases_args *args = (dense_qp_qpoases_args *)args_;
    dense_qp_qpoases_memory *memory = (dense_qp_qpoases_memory *)memory_;

    // initialize return code
    int acados_status = ACADOS_SUCCESS;

    // extract qpoases data
    double *H = memory->H;
    double *A = memory->A;
    double *C = memory->C;
    double *g = memory->g;
    double *b = memory->b;
    double *d_lb0 = memory->d_lb0;
    double *d_ub0 = memory->d_ub0;
    double *d_lb = memory->d_lb;
    double *d_ub = memory->d_ub;
    double *d_lg = memory->d_lg;
    double *d_ug = memory->d_ug;
    int *idxb = memory->idxb;
    double *prim_sol = memory->prim_sol;
    double *dual_sol = memory->dual_sol;
    QProblemB *QPB = memory->QPB;
    QProblem *QP = memory->QP;

    // extract dense qp size
    int nvd = qp_in->dim->nv;
    int ned = qp_in->dim->ne;
    int ngd = qp_in->dim->ng;
    int nbd = qp_in->dim->nb;

    assert(ned == 0 && "ned != 0 not supported yet");

    // fill in the upper triangular of H in dense_qp
    dtrtr_l_libstr(nvd, qp_in->Hv, 0, 0, qp_in->Hv, 0, 0);

    // dense qp row-major
    d_cvt_dense_qp_to_rowmaj(qp_in, H, g, A, b, idxb, d_lb0, d_ub0, C, d_lg, d_ug,
        NULL, NULL, NULL, NULL, NULL);

    // reorder bounds
    for (int ii = 0; ii < nvd; ii++) {
        d_lb[ii] = -QPOASES_INFTY;
        d_ub[ii] = +QPOASES_INFTY;
    }
    for (int ii = 0; ii < nbd; ii++) {
        d_lb[idxb[ii]] = d_lb0[ii];
        d_ub[idxb[ii]] = d_ub0[ii];
    }

    // cholesky factorization of H
    // dpotrf_l_libstr(nvd, qpd->Hv, 0, 0, sR, 0, 0);

    // fill in upper triangular of R
    // dtrtr_l_libstr(nvd, sR, 0, 0, sR, 0, 0);

    // extract R
    // d_cvt_strmat2mat(nvd, nvd, sR, 0, 0, R, nvd);

#if 0
#endif

    // cold start the dual solution with no active constraints
    int warm_start = args->warm_start;
    if (!warm_start) {
        for (int ii = 0; ii < nvd + ngd; ii++)
            dual_sol[ii] = 0;
    }

    info->interface_time = acados_toc(&interface_timer);
    acados_tic(&qp_timer);

    // solve dense qp
    int nwsr = args->max_nwsr;
    double cputime = args->max_cputime;
    int return_flag = 0;
    if (ngd > 0) {  // QProblem
        QProblemCON(QP, nvd, ngd, HST_POSDEF);
        QProblem_setPrintLevel(QP, PL_MEDIUM);
        QProblem_printProperties(QP);
        return_flag = QProblem_initW(QP, H, g, C, d_lb, d_ub, d_lg, d_ug, &nwsr, &cputime,
            NULL, dual_sol, NULL, NULL, NULL);  // NULL or 0
        //            NULL, NULL, NULL, NULL);
        //            NULL, NULL, NULL, R);  // to provide Cholesky factor
        QProblem_getPrimalSolution(QP, prim_sol);
        QProblem_getDualSolution(QP, dual_sol);
    } else {  // QProblemB
        QProblemBCON(QPB, nvd, HST_POSDEF);
        QProblemB_setPrintLevel(QPB, PL_MEDIUM);
        QProblemB_printProperties(QPB);
        return_flag = QProblemB_initW(QPB, H, g, d_lb, d_ub, &nwsr, &cputime,
            NULL, dual_sol, NULL, NULL);  // NULL or 0
        QProblemB_getPrimalSolution(QPB, prim_sol);
        QProblemB_getDualSolution(QPB, dual_sol);
    }

    // save solution statistics to memory
    memory->cputime = cputime;
    memory->nwsr = nwsr;

#if 0
    double cmpl = 0.0, feas = 0.0, stat = 0.0;
    qpOASES_getKktViolation(nvd, ngd, H, g, C, d_lb, d_ub, d_lg, d_ug, prim_sol, dual_sol, &stat, &feas, &cmpl);
    printf("\nstat=%e, feas=%e, cmpl=%e\n", stat, feas, cmpl);
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

    // return
    // TODO(dimitris): cast qpoases return to acados return
    acados_status = return_flag;

    info->interface_time += acados_toc(&interface_timer);
    info->total_time = acados_toc(&tot_timer);

    // printf("total time = \t\t\t%f\n", 1000*info->total_time);
    // printf("interface time = \t\t%f\n", 1000*info->interface_time);
    // printf("qp time = \t\t\t%f\n", 1000*info->solve_QP_time);
    // printf("total time from qpOASES = \t%f\n", 1000*cputime);  // does not include getSolution
    // printf("**************\n");

    return acados_status;
}
