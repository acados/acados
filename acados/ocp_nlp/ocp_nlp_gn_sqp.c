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

#include "acados/ocp_nlp/ocp_nlp_gn_sqp.h"

// external
#include <assert.h>
#include <stdio.h>
#if defined(EXT_DEPS)
#include <stdlib.h>
#endif
// blasfeo
#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_blas.h"
// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_common_ext_dep.h"
#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/sim/sim_common.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"
#include "acados/utils/mem.h"

// TEMP
// #include "acados/ocp_qp/ocp_qp_condensing_qpoases.h"


int ocp_nlp_gn_sqp_calculate_args_size(ocp_nlp_dims *dims, ocp_qp_solver *qp_solver)
{
    int size = 0;

    ocp_qp_dims qp_dims;
    cast_nlp_dims_to_qp_dims(&qp_dims, dims);

    size += sizeof(ocp_nlp_gn_sqp_args);
    size += sizeof(ocp_nlp_args);  // TODO(dimitris): REPLACE WITH CALCULATE SIZE?
    size += qp_solver->calculate_args_size(&qp_dims);

    return size;
}



ocp_nlp_gn_sqp_args *ocp_nlp_gn_sqp_assign_args(ocp_nlp_dims *dims, ocp_qp_solver *qp_solver, void *raw_memory)
{
    ocp_nlp_gn_sqp_args *args;

    ocp_qp_dims qp_dims;
    cast_nlp_dims_to_qp_dims(&qp_dims, dims);

    char *c_ptr = (char *) raw_memory;

    args = (ocp_nlp_gn_sqp_args *) c_ptr;
    c_ptr += sizeof(ocp_nlp_gn_sqp_args);

    args->common = (ocp_nlp_args *) c_ptr;
    c_ptr += sizeof(ocp_nlp_args);  // TODO(dimitris): REPLACE WITH ASSIGN?

    args->qp_solver_args = qp_solver->assign_args(&qp_dims, c_ptr);
    c_ptr += qp_solver->calculate_args_size(&qp_dims);  // TODO(dimitris): replace with memsize?

    assert((char*)raw_memory + ocp_nlp_gn_sqp_calculate_args_size(dims, qp_solver) >= c_ptr);

    return args;
}



int ocp_nlp_gn_sqp_calculate_memory_size(ocp_nlp_dims *dims, ocp_qp_solver *qp_solver, ocp_nlp_gn_sqp_args *args)
{
    int size = 0;

    ocp_qp_dims qp_dims;
    cast_nlp_dims_to_qp_dims(&qp_dims, dims);

    size+= sizeof(ocp_nlp_gn_sqp_memory);
    size+= ocp_nlp_calculate_memory_size(dims, args->common);
    size+= qp_solver->calculate_memory_size(&qp_dims, args->qp_solver_args);

    make_int_multiple_of(8, &size);
    size += 2 * 8;

    return size;
}



ocp_nlp_gn_sqp_memory *ocp_nlp_gn_sqp_assign_memory(ocp_nlp_dims *dims, ocp_qp_solver *qp_solver, ocp_nlp_gn_sqp_args *args, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_qp_dims qp_dims;
    cast_nlp_dims_to_qp_dims(&qp_dims, dims);

    ocp_nlp_gn_sqp_memory *mem = (ocp_nlp_gn_sqp_memory *)c_ptr;
    c_ptr += sizeof(ocp_nlp_gn_sqp_memory);

    align_char_to(8, &c_ptr);
    mem->common = ocp_nlp_assign_memory(dims, args->common, c_ptr);
    c_ptr += ocp_nlp_calculate_memory_size(dims, args->common);

    align_char_to(8, &c_ptr);
    mem->qp_mem = qp_solver->assign_memory(&qp_dims, args->qp_solver_args, c_ptr);
    c_ptr += qp_solver->calculate_memory_size(&qp_dims, args->qp_solver_args);

    mem->qp_solver = qp_solver;
    mem->qp_solver->mem = mem->qp_mem;
    mem->dims = dims;

    assert((char *)raw_memory + ocp_nlp_gn_sqp_calculate_memory_size(dims, qp_solver, args) >= c_ptr);

    return mem;
}


#if defined(EXT_DEPS)

ocp_nlp_gn_sqp_args *ocp_nlp_gn_sqp_create_args(ocp_nlp_dims *dims, ocp_qp_solver *qp_solver)
{
    int size = ocp_nlp_gn_sqp_calculate_args_size(dims, qp_solver);
    void *ptr = malloc(size);
    ocp_nlp_gn_sqp_args *args = ocp_nlp_gn_sqp_assign_args(dims, qp_solver, ptr);

    // TODO(dimitris): nest in initialize default args of each module
    qp_solver->initialize_default_args(args->qp_solver_args);
    args->common->maxIter = 30;

    return args;
}



ocp_nlp_gn_sqp_memory *new_ocp_nlp_gn_sqp_create_memory(ocp_nlp_dims *dims, ocp_qp_solver *qp_solver, ocp_nlp_gn_sqp_args *args)
{
    int size = ocp_nlp_gn_sqp_calculate_memory_size(dims, qp_solver, args);
    void *ptr = malloc(size);
    ocp_nlp_gn_sqp_memory *mem = ocp_nlp_gn_sqp_assign_memory(dims, qp_solver, args, ptr);
    return mem;
}

#endif


int ocp_nlp_gn_sqp_calculate_workspace_size(ocp_nlp_dims *dims, ocp_qp_solver *qp_solver, ocp_nlp_gn_sqp_args *args)
{
    int size = 0;

    ocp_qp_dims qp_dims;
    cast_nlp_dims_to_qp_dims(&qp_dims, dims);

    size += sizeof(ocp_nlp_gn_sqp_work);
    size += ocp_nlp_calculate_workspace_size(dims, args->common);
    size += ocp_qp_in_calculate_size(&qp_dims);
    size += ocp_qp_out_calculate_size(&qp_dims);

    size += (dims->N+1)*sizeof(struct d_strvec);  // tmp_vecs

    for (int ii = 0; ii < dims->N+1; ii++)
    {
        size += d_size_strvec(qp_dims.nx[ii] + qp_dims.nu[ii]);
    }

    make_int_multiple_of(64, &size);
    size += 1 * 64;

    return size;
}



static void ocp_nlp_gn_sqp_cast_workspace(ocp_nlp_gn_sqp_work *work, ocp_nlp_gn_sqp_memory *mem, ocp_nlp_gn_sqp_args *args)
{
    char *c_ptr = (char *)work;
    c_ptr += sizeof(ocp_nlp_gn_sqp_work);

    int N = mem->dims->N;
    int *nx = mem->dims->nx;
    int *nu = mem->dims->nu;

    ocp_qp_dims qp_dims;
    cast_nlp_dims_to_qp_dims(&qp_dims, mem->dims);

    // set up common nlp workspace
    work->common = (ocp_nlp_work *)c_ptr;
    ocp_nlp_cast_workspace(work->common, mem->common);
    c_ptr += ocp_nlp_calculate_workspace_size(mem->dims, args->common);

    // set up local SQP data
    work->tmp_vecs = (struct d_strvec *)c_ptr;
    c_ptr += (N+1)*sizeof(struct d_strvec);

    align_char_to(64, &c_ptr);
    for (int ii = 0; ii < N+1; ii++)
    {
        d_create_strvec(nx[ii]+nu[ii], &work->tmp_vecs[ii], c_ptr);
        c_ptr += work->tmp_vecs[ii].memory_size;
    }

    // set up QP solver
    work->qp_in = ocp_qp_in_assign(&qp_dims, c_ptr);
    c_ptr += ocp_qp_in_calculate_size(&qp_dims);
    work->qp_out = ocp_qp_out_assign(&qp_dims, c_ptr);
    c_ptr += ocp_qp_out_calculate_size(&qp_dims);
    mem->qp_solver->qp_in = work->qp_in;
    mem->qp_solver->qp_out = work->qp_out;
    mem->qp_solver->args = args->qp_solver_args;

    assert((char *)work + ocp_nlp_gn_sqp_calculate_workspace_size(mem->dims, mem->qp_solver, args) >= c_ptr);
}



static void initialize_objective(const ocp_nlp_in *nlp_in, ocp_nlp_gn_sqp_memory *gn_sqp_mem, ocp_nlp_gn_sqp_work *work)
{
    int N = nlp_in->dims->N;
    int *nx = nlp_in->dims->nx;
    int *nu = nlp_in->dims->nu;
    ocp_nlp_ls_cost *cost = (ocp_nlp_ls_cost*) nlp_in->cost;

    // real_t **qp_Q = (real_t **) gn_sqp_mem->qp_solver->qp_in->Q;
    // real_t **qp_S = (real_t **) gn_sqp_mem->qp_solver->qp_in->S;
    // real_t **qp_R = (real_t **) gn_sqp_mem->qp_solver->qp_in->R;
	struct d_strmat *sRSQrq = gn_sqp_mem->qp_solver->qp_in->RSQrq;

    // TODO(rien): only for least squares cost with state and control reference atm
    for (int_t i = 0; i <= N; i++) {

        // // TODO(dimitris): DOUBLE CHECK THAT IT'S CORRECT ALSO FOR FULL HESSIANS!
        // cost->W[i][nx[i]] = 66;
        // cost->W[i][nx[i]+1] = 77;

        // copy R
        d_cvt_mat2strmat(nu[i], nu[i], &cost->W[i][nx[i]*(nx[i]+nu[i])+nx[i]], nx[i] + nu[i], &sRSQrq[i], 0, 0);
        // copy Q
        d_cvt_mat2strmat(nx[i], nx[i], &cost->W[i][0], nx[i] + nu[i], &sRSQrq[i], nu[i], nu[i]);
        // copy S
        d_cvt_tran_mat2strmat(nu[i], nx[i], &cost->W[i][nx[i]], nx[i] + nu[i], &sRSQrq[i], nu[i], 0);

        // printf("W = \n");
        // d_print_mat(nx[i]+nu[i], nx[i]+nu[i], cost->W[i], nx[i]+nu[i]);

        // printf("RSQrq=\n");
        // d_print_strmat(nx[i]+nu[i]+1, nx[i]+nu[i], &sRSQrq[i], 0, 0);

        // for (int_t j = 0; j < nx[i]; j++) {
        //     for (int_t k = 0; k < nx[i]; k++) {
        //         qp_Q[i][j * nx[i] + k] = cost->W[i][j * (nx[i] + nu[i]) + k];
        //     }
        //     for (int_t k = 0; k < nu[i]; k++) {
        //         qp_S[i][j * nu[i] + k] =
        //             cost->W[i][j * (nx[i] + nu[i]) + nx[i] + k];
        //     }
        // }
        // for (int_t j = 0; j < nu[i]; j++) {
        //     for (int_t k = 0; k < nu[i]; k++) {
        //         qp_R[i][j * nu[i] + k] =
        //             cost->W[i][(nx[i] + j) * (nx[i] + nu[i]) + nx[i] + k];
        //     }
        // }
    }
}



static void initialize_trajectories(
    const ocp_nlp_in *nlp_in,
    ocp_nlp_gn_sqp_memory *gn_sqp_mem,
    ocp_nlp_gn_sqp_work *work) {

    int N = nlp_in->dims->N;
    int *nx = nlp_in->dims->nx;
    int *nu = nlp_in->dims->nu;
    real_t *w = work->common->w;

    int_t w_idx = 0;
    for (int_t i = 0; i <= N; i++) {
        for (int_t j = 0; j < nx[i]; j++) {
            w[w_idx + j] = gn_sqp_mem->common->x[i][j];
        }
        for (int_t j = 0; j < nu[i]; j++) {
            w[w_idx + nx[i] + j] = gn_sqp_mem->common->u[i][j];
        }
        w_idx += nx[i] + nu[i];
    }
}



static void multiple_shooting(const ocp_nlp_in *nlp, ocp_nlp_gn_sqp_memory *mem, ocp_nlp_gn_sqp_work *work) {

    int N = nlp->dims->N;
    int *nx = nlp->dims->nx;
    int *nu = nlp->dims->nu;
    int *nb = nlp->dims->nb;
    int *ng = nlp->dims->ng;

    real_t *w = work->common->w;
    struct d_strvec *stmp = work->tmp_vecs;

    sim_solver *sim = nlp->sim;
    ocp_nlp_ls_cost *cost = (ocp_nlp_ls_cost *) nlp->cost;
    real_t **y_ref = cost->y_ref;

    struct d_strmat *sBAbt = mem->qp_solver->qp_in->BAbt;
    struct d_strvec *sb = mem->qp_solver->qp_in->b;
    struct d_strmat *sRSQrq = mem->qp_solver->qp_in->RSQrq;
    struct d_strvec *srq = mem->qp_solver->qp_in->rq;
    struct d_strvec *sd = mem->qp_solver->qp_in->d;

    // real_t **qp_A = (real_t **) mem->qp_solver->qp_in->A;
    // real_t **qp_B = (real_t **) mem->qp_solver->qp_in->B;
    // real_t **qp_b = (real_t **) mem->qp_solver->qp_in->b;
    // real_t **qp_q = (real_t **) mem->qp_solver->qp_in->q;
    // real_t **qp_r = (real_t **) mem->qp_solver->qp_in->r;
    // real_t **qp_lb = (real_t **) mem->qp_solver->qp_in->lb;
    // real_t **qp_ub = (real_t **) mem->qp_solver->qp_in->ub;

    int_t w_idx = 0;

    for (int_t i = 0; i < N; i++)
    {
        // Pass state and control to integrator
        for (int_t j = 0; j < nx[i]; j++) sim[i].in->x[j] = w[w_idx+j];
        for (int_t j = 0; j < nu[i]; j++) sim[i].in->u[j] = w[w_idx+nx[i]+j];
        sim[i].fun(sim[i].in, sim[i].out, sim[i].args, sim[i].mem, sim[i].work);

        // TODO(rien): transition functions for changing dimensions not yet implemented!

        // convert b
        d_cvt_vec2strvec(nx[i+1], sim[i].out->xn, &sb[i], 0);
        for (int j = 0; j < nx[i+1]; j++)
            DVECEL_LIBSTR(&sb[i], j) -= w[w_idx+nx[i]+nu[i]+j];
        // copy B
        d_cvt_tran_mat2strmat(nx[i+1], nu[i], &sim[i].out->S_forw[nx[i+1]*nx[i]], nx[i+1], &sBAbt[i], 0, 0);
        // copy A
        d_cvt_tran_mat2strmat(nx[i+1], nx[i], &sim[i].out->S_forw[0], nx[i+1], &sBAbt[i], nu[i], 0);
        // copy b
        drowin_libstr(nx[i+1], 1.0, &sb[i], 0, &sBAbt[i], nu[i]+nx[i], 0);

        // printf("AB = \n");
        // d_print_mat(nx[i+1], nx[i]+nu[i], sim[i].out->S_forw, nx[i+1]);

        // printf("ABbt=\n");
        // d_print_strmat(nx[i]+nu[i]+1, nx[i+1], &sBAbt[i], 0, 0);

        // for (int_t j = 0; j < nx[i]; j++) {
        //     qp_b[i][j] = sim[i].out->xn[j] - w[w_idx+nx[i]+nu[i]+j];
        //     for (int_t k = 0; k < nx[i]; k++)
        //         qp_A[i][j*nx[i]+k] = sim[i].out-> [j*nx[i]+k];
        // }
        // for (int_t j = 0; j < nu[i]; j++)
        //     for (int_t k = 0; k < nx[i]; k++)
        //         qp_B[i][j*nx[i]+k] = sim[i].out->S_forw[(nx[i]+j)*nx[i]+k];

        // Update bounds:
        for (int_t j = 0; j < nb[i]; j++) {
// #ifdef FLIP_BOUNDS

            // NOTE!!!! ATM IDXB OF NLP IS FLIPPED WRT QP
            DVECEL_LIBSTR(&sd[i], j) = nlp->lb[i][j] - w[w_idx+nlp->idxb[i][j]];
            DVECEL_LIBSTR(&sd[i], j+nb[i]+ng[i]) = nlp->ub[i][j] - w[w_idx+nlp->idxb[i][j]];

            // if (nlp->idxb[i][j] < nu[i]) {
            //     DVECEL_LIBSTR(&sd[i], j) = nlp->lb[i][j] - w[w_idx + nx[i] + nlp->idxb[i][j]];
            //     DVECEL_LIBSTR(&sd[i], j+nb[i]+ng[i]) = nlp->ub[i][j] - w[w_idx + nx[i] + nlp->idxb[i][j]];
            //     // qp_lb[i][j] = nlp->lb[i][j] - w[w_idx + nx[i] + nlp->idxb[i][j]];
            //     // qp_ub[i][j] = nlp->ub[i][j] - w[w_idx + nx[i] + nlp->idxb[i][j]];
            // } else {
            //     DVECEL_LIBSTR(&sd[i], j) = nlp->lb[i][j] - w[w_idx - nu[i] + nlp->idxb[i][j]];
            //     DVECEL_LIBSTR(&sd[i], j+nb[i]+ng[i]) = nlp->ub[i][j] - w[w_idx - nu[i] + nlp->idxb[i][j]];
            //     // qp_lb[i][j] = nlp->lb[i][j] - w[w_idx - nu[i] + nlp->idxb[i][j]];
            //     // qp_ub[i][j] = nlp->ub[i][j] - w[w_idx - nu[i] + nlp->idxb[i][j]];
            // }

// #else
//             qp_lb[i][j] = nlp->lb[i][j] - w[w_idx+nlp->idxb[i][j]];
//             qp_ub[i][j] = nlp->ub[i][j] - w[w_idx+nlp->idxb[i][j]];
// #endif
        }

        // Update gradients
        // TODO(rien): only for diagonal Q, R matrices atm
        // TODO(rien): only for least squares cost with state and control reference atm
        sim_RK_opts *opts = (sim_RK_opts*) sim[i].args;

        for (int j = 0; j < nx[i]; j++)
            DVECEL_LIBSTR(&stmp[i], nu[i]+j) = w[w_idx+j]-y_ref[i][j];
        for (int j = 0; j < nu[i]; j++)
            DVECEL_LIBSTR(&stmp[i], j) = w[w_idx+nx[i]+j]-y_ref[i][nx[i]+j];

        dsymv_l_libstr(nu[i]+nx[i], nu[i]+nx[i], 1.0, &sRSQrq[i], 0, 0, &stmp[i], 0, 0.0, &srq[i], 0, &srq[i], 0);

        if (opts->scheme.type != exact)
        {
            for (int_t j = 0; j < nx[i]; j++)
                DVECEL_LIBSTR(&srq[i], nu[i]+j) += sim[i].out->grad[j];
            for (int_t j = 0; j < nu[i]; j++)
                DVECEL_LIBSTR(&srq[i], j) += sim[i].out->grad[nx[i]+j];
        }
        drowin_libstr(nu[i]+nx[i], 1.0, &srq[i], 0, &sRSQrq[i], nu[i]+nx[i], 0);

        // for (int_t j = 0; j < nx[i]; j++) {
        //     qp_q[i][j] = cost->W[i][j*(nx[i]+nu[i]+1)]*(w[w_idx+j]-y_ref[i][j]);
        //     // adjoint-based gradient correction:
        //     if (opts->scheme.type != exact) qp_q[i][j] += sim[i].out->grad[j];
        // }
        // for (int_t j = 0; j < nu[i]; j++) {
        //     qp_r[i][j] = cost->W[i][(nx[i]+j)*(nx[i]+nu[i]+1)]*(w[w_idx+nx[i]+j]-y_ref[i][nx[i]+j]);
        //     // adjoint-based gradient correction:
        //     if (opts->scheme.type != exact) qp_r[i][j] += sim[i].out->grad[nx[i]+j];
        // }
        w_idx += nx[i]+nu[i];
    }

    for (int_t j = 0; j < nb[N]; j++) {
// #ifdef FLIP_BOUNDS
        if (nlp->idxb[N][j] < nu[N]) {
            DVECEL_LIBSTR(&sd[N], j) = nlp->lb[N][j] - w[w_idx + nx[N] + nlp->idxb[N][j]];
            DVECEL_LIBSTR(&sd[N], j+nb[N]+ng[N]) = nlp->ub[N][j] - w[w_idx + nx[N] + nlp->idxb[N][j]];
        } else {
            DVECEL_LIBSTR(&sd[N], j) =  nlp->lb[N][j] - w[w_idx - nu[N] + nlp->idxb[N][j]];
            DVECEL_LIBSTR(&sd[N], j+nb[N]+ng[N]) = nlp->ub[N][j] - w[w_idx - nu[N] + nlp->idxb[N][j]];
        }
// #else
//         qp_lb[N][j] = nlp->lb[N][j] - w[w_idx+nlp->idxb[N][j]];
//         qp_ub[N][j] = nlp->ub[N][j] - w[w_idx+nlp->idxb[N][j]];
// #endif

    }

    for (int j=0; j<nx[N]; j++)
        DVECEL_LIBSTR(&stmp[N], j) = w[w_idx+j]-y_ref[N][j];
    dsymv_l_libstr(nx[N], nx[N], 1.0, &sRSQrq[N], 0, 0, &stmp[N], 0, 0.0, &srq[N], 0, &srq[N], 0);

    drowin_libstr(nx[N], 1.0, &srq[N], 0, &sRSQrq[N], nx[N], 0);

    // for (int_t j = 0; j < nx[N]; j++)
    //     qp_q[N][j] = cost->W[N][j*(nx[N]+nu[N]+1)]*(w[w_idx+j]-y_ref[N][j]);
}



static void update_variables(const ocp_nlp_in *nlp, ocp_nlp_gn_sqp_memory *mem, real_t *w) {
    int N = nlp->dims->N;
    int *nx = nlp->dims->nx;
    int *nu = nlp->dims->nu;
    sim_solver *sim = nlp->sim;

    for (int_t i = 0; i < N; i++)
        for (int_t j = 0; j < nx[i+1]; j++)
            sim[i].in->S_adj[j] = -DVECEL_LIBSTR(&mem->qp_solver->qp_out->pi[i], j);
            // sim[i].in->S_adj[j] = -mem->qp_solver->qp_out->pi[i][j];

    int_t w_idx = 0;
    for (int_t i = 0; i <= N; i++) {
        for (int_t j = 0; j < nx[i]; j++) {
            w[w_idx+j] += DVECEL_LIBSTR(&mem->qp_solver->qp_out->ux[i], j + nu[i]);
            // w[w_idx+j] += mem->qp_solver->qp_out->x[i][j];
        }
        for (int_t j = 0; j < nu[i]; j++)
            w[w_idx+nx[i]+j] += DVECEL_LIBSTR(&mem->qp_solver->qp_out->ux[i], j);
            // w[w_idx+nx[i]+j] += mem->qp_solver->qp_out->u[i][j];
        w_idx += nx[i]+nu[i];
    }
}



static void store_trajectories(const ocp_nlp_in *nlp, ocp_nlp_memory *memory, ocp_nlp_out *out,
    real_t *w) {

    int N = nlp->dims->N;
    int *nx = nlp->dims->nx;
    int *nu = nlp->dims->nu;

    int_t w_idx = 0;
    for (int_t i = 0; i <= N; i++) {
        for (int_t j = 0; j < nx[i]; j++) {
            memory->x[i][j] = w[w_idx+j];
            out->x[i][j] = w[w_idx+j];
        }
        for (int_t j = 0; j < nu[i]; j++) {
            memory->u[i][j] = w[w_idx+nx[i]+j];
            out->u[i][j] = w[w_idx+nx[i]+j];
        }
        w_idx += nx[i] + nu[i];
    }
}



// Simple fixed-step Gauss-Newton based SQP routine
int ocp_nlp_gn_sqp(ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out, ocp_nlp_gn_sqp_args *args, ocp_nlp_gn_sqp_memory *mem, void *work_)
{
    ocp_nlp_gn_sqp_work *work = (ocp_nlp_gn_sqp_work*) work_;
    ocp_nlp_gn_sqp_cast_workspace(work, mem, args);

    int N = nlp_in->dims->N;
    int *nx = nlp_in->dims->nx;
    int *nu = nlp_in->dims->nu;
    int *nb = nlp_in->dims->nb;

    initialize_objective(nlp_in, mem, work);

    initialize_trajectories(nlp_in, mem, work);

    // TODO(dimitris): move somewhere else
    int_t **qp_idxb = (int_t **) mem->qp_solver->qp_in->idxb;
    for (int_t i = 0; i <= N; i++) {
        for (int_t j = 0; j < nb[i]; j++) {
            if (nlp_in->idxb[i][j] < nx[i]) {
                // state bound
                qp_idxb[i][j] = nlp_in->idxb[i][j]+nu[i];
            } else {
                // input bound
                qp_idxb[i][j] = nlp_in->idxb[i][j]-nx[i];
            }
        }
    }

    int_t max_sqp_iterations =  args->common->maxIter;

    acados_timer timer;
    real_t total_time = 0;
    acados_tic(&timer);
    for (int_t sqp_iter = 0; sqp_iter < max_sqp_iterations; sqp_iter++)
    {
        multiple_shooting(nlp_in, mem, work);

        int_t qp_status = mem->qp_solver->fun(mem->qp_solver->qp_in, mem->qp_solver->qp_out,
            mem->qp_solver->args, mem->qp_solver->mem);

        // ocp_qp_condensing_qpoases_memory *tmp_mem = mem->qp_mem;
        // printf("SQP iter #%d - QP status %d - nwsr %d - cputime %f\n", sqp_iter, qp_status, tmp_mem->solver_memory->nwsr, tmp_mem->solver_memory->cputime*1e3);

        if (qp_status != 0)
        {
            printf("QP solver returned error status %d\n", qp_status);
            return -1;
        }

        update_variables(nlp_in, mem, work->common->w);

        for (int_t i = 0; i < N; i++)
        {
            sim_RK_opts *opts = (sim_RK_opts*) nlp_in->sim[i].args;
            nlp_in->sim[i].in->sens_adj = (opts->scheme.type != exact);
            if (nlp_in->freezeSens)
            {  // freeze inexact sensitivities after first SQP iteration !!
                opts->scheme.freeze = true;
            }
        }
    }

    total_time += acados_toc(&timer);
    store_trajectories(nlp_in, mem->common, nlp_out, work->common->w);

    return 0;
}
