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
#if defined(RUNTIME_CHECKS)
#include <assert.h>
#endif
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


int dense_qp_qpoases_calculate_args_size(dense_qp_dims *dims)
{
    int size = 0;
    size += sizeof(dense_qp_qpoases_args);

    make_int_multiple_of(8, &size);
    size += 1 * 8;

    return size;
}



dense_qp_qpoases_args *dense_qp_qpoases_assign_args(dense_qp_dims *dims, void *raw_memory)
{
    dense_qp_qpoases_args *args;

    char *c_ptr = (char *) raw_memory;

    args = (dense_qp_qpoases_args *) c_ptr;
    c_ptr += sizeof(dense_qp_qpoases_args);

    assert((char*)raw_memory + dense_qp_qpoases_calculate_args_size(dims) >= c_ptr);

    return args;
}



void dense_qp_qpoases_initialize_default_args(dense_qp_qpoases_args *args)
{
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
    size += (2 * nvd + 2 * ngd) * sizeof(double);  // dual_sol

    if (ngd > 0)  // QProblem
        size += QProblem_calculateMemorySize(nvd, ngd);
    else  // QProblemB
        size += QProblemB_calculateMemorySize(nvd);

    make_int_multiple_of(8, &size);
    size += 1 * 8;

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

    align_char_to(8, &c_ptr);

    //
    mem->H = (double *)c_ptr;
    c_ptr += nvd * nvd * sizeof(double);
    //
    mem->A = (double *)c_ptr;
    c_ptr += nvd * ned * sizeof(double);
    //
    mem->C = (double *)c_ptr;
    c_ptr += nvd * ngd * sizeof(double);
    //
    mem->g = (double *)c_ptr;
    c_ptr += nvd * sizeof(double);
    // TODO(dimitris): use this instead
    // assign_double(nvd, &mem->g, &c_ptr);
    //
    mem->b = (double *)c_ptr;
    c_ptr += ned * sizeof(double);
    //
    mem->d_lb0 = (double *)c_ptr;
    c_ptr += nbd * sizeof(double);
    //
    mem->d_ub0 = (double *)c_ptr;
    c_ptr += nbd * sizeof(double);
    //
    mem->d_lb = (double *)c_ptr;
    c_ptr += nvd * sizeof(double);
    //
    mem->d_ub = (double *)c_ptr;
    c_ptr += nvd * sizeof(double);
    //
    mem->d_lg = (double *)c_ptr;
    c_ptr += ngd * sizeof(double);
    //
    mem->d_ug = (double *)c_ptr;
    c_ptr += ngd * sizeof(double);
    //
    mem->prim_sol = (double *)c_ptr;
    c_ptr += nvd * sizeof(double);
    //
    mem->dual_sol = (double *)c_ptr;
    c_ptr += (2 * nvd + 2 * ngd) * sizeof(double);


    // TODO(dimitris): update syntax in qpOASES

    if (ngd > 0) {  // QProblem
        c_ptr = QProblem_assignMemory(nvd, ngd, (QProblem **) &(mem->QP), c_ptr);
    } else {  // QProblemB
        c_ptr = QProblemB_assignMemory(nvd, (QProblemB **) &(mem->QPB), c_ptr);
    }

    // int stuff
    mem->idxb = (int *)c_ptr;
    c_ptr += nbd * sizeof(int);

    assert((char *)raw_memory + dense_qp_qpoases_calculate_memory_size(dims, args_) >= c_ptr);

    return mem;
}



int dense_qp_qpoases(dense_qp_in *qp_in, dense_qp_out *qp_out, void *args_, void *memory_)
{
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

    // fill in the upper triangular of H in dense_qp
    dtrtr_l_libstr(nvd, qp_in->Hg, 0, 0, qp_in->Hg, 0, 0);

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
    //  dpotrf_l_libstr(nvd, qpd->Hg, 0, 0, sR, 0, 0);

    //  fill in upper triangular of R
    //  dtrtr_l_libstr(nvd, sR, 0, 0, sR, 0, 0);

    //  extract R
    //  d_cvt_strmat2mat(nvd, nvd, sR, 0, 0, R, nvd);

#if 0
    d_print_mat(nvd, nvd, H, nvd);
    d_print_mat(nvd, nvd, R, nvd);
    exit(1);
#endif

    // cold start the dual solution with no active constraints
    int warm_start = args->warm_start;
    if (!warm_start) {
        for (int ii = 0; ii < 2 * nvd + 2 * ngd; ii++)
            dual_sol[ii] = 0;
    }

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
    d_print_mat(1, nvd, prim_sol, 1);
    exit(1);
#endif

    // copy prim_sol and dual_sol to qpd_sol
    d_cvt_vec2strvec(nvd, prim_sol, qp_out->v, 0);
    for (int ii = 0; ii < 2*nbd+2*ngd; ii++)
        qp_out->lam->pa[ii] = 0.0;
    for (int ii = 0; ii < nbd; ii++) {
        if (dual_sol[ii] >= 0.0)
            qp_out->lam->pa[ii] = dual_sol[ii];
        else
            qp_out->lam->pa[nbd+ngd+ii] = - dual_sol[ii];
    }
    for (int ii = 0; ii < ngd; ii++) {
        if (dual_sol[nbd+ii] >= 0.0)
            qp_out->lam->pa[nbd+ii] =   dual_sol[nbd+ii];
        else
            qp_out->lam->pa[2*nbd+ngd+ii] = - dual_sol[nbd+ii];
    }

    // return
    // TODO(dimitris): cast qpoases return to acados return
    acados_status = return_flag;
    return acados_status;
}
