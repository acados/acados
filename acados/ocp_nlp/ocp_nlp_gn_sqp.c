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
#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/sim/sim_casadi_wrapper.h"
#include "acados/sim/sim_common.h"
#include "acados/sim/sim_collocation.h"
#include "acados/utils/create.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"
#include "acados/utils/mem.h"

static int get_max_sim_workspace_size(ocp_nlp_dims *dims, ocp_nlp_gn_sqp_args *args)
{
    sim_dims sim_dims;
    int sim_work_size;

    int max_sim_work_size = 0;

    for (int ii = 0; ii < dims->N; ii++)
    {
        cast_nlp_dims_to_sim_dims(&sim_dims, dims, ii);
        // sim_in_size = sim_in_calculate_size(&sim_dims);
        // if (sim_in_size > *max_sim_in_size) *max_sim_in_size = sim_in_size;
        // sim_out_size = sim_out_calculate_size(&sim_dims);
        // if (sim_out_size > *max_sim_out_size) *max_sim_out_size = sim_out_size;
        sim_work_size = args->sim_solvers[ii]->calculate_workspace_size(&sim_dims, args->sim_solvers_args[ii]);
        if (sim_work_size > max_sim_work_size) max_sim_work_size = sim_work_size;
    }
    return max_sim_work_size;
}


int ocp_nlp_gn_sqp_calculate_args_size(ocp_nlp_dims *dims, ocp_qp_xcond_solver *qp_solver, sim_solver *sim_solvers)
{
    int size = 0;

    ocp_qp_dims qp_dims;
    cast_nlp_dims_to_qp_dims(&qp_dims, dims);
    size += sizeof(ocp_nlp_gn_sqp_args);
    size += sizeof(ocp_qp_xcond_solver);

    size += qp_solver->calculate_args_size(&qp_dims, qp_solver->qp_solver_funs);

    sim_dims sim_dims;

    size += dims->N*sizeof(sim_solver *);
    size += dims->N*sizeof(void *);  //sim_solvers_args

    int return_value;

    for (int ii = 0; ii < dims->N; ii++)
    {
        cast_nlp_dims_to_sim_dims(&sim_dims, dims, ii);
        size += sizeof(sim_solver);
        size += sim_solvers[ii].calculate_args_size(&sim_dims);
    }

    return size;
}



ocp_nlp_gn_sqp_args *ocp_nlp_gn_sqp_assign_args(ocp_nlp_dims *dims, ocp_qp_xcond_solver *qp_solver, sim_solver *sim_solvers, void *raw_memory)
{
    ocp_nlp_gn_sqp_args *args;

    ocp_qp_dims qp_dims;
    cast_nlp_dims_to_qp_dims(&qp_dims, dims);

    char *c_ptr = (char *) raw_memory;

    args = (ocp_nlp_gn_sqp_args *) c_ptr;
    c_ptr += sizeof(ocp_nlp_gn_sqp_args);

    args->qp_solver = (ocp_qp_xcond_solver*) c_ptr;
    c_ptr += sizeof(ocp_qp_xcond_solver);

    copy_module_pointers_to_args(args->qp_solver, qp_solver);
    args->qp_solver->qp_solver_funs = qp_solver->qp_solver_funs;
    copy_module_pointers_to_args(args->qp_solver->qp_solver_funs, qp_solver->qp_solver_funs);

    args->qp_solver_args = args->qp_solver->assign_args(&qp_dims, qp_solver->qp_solver_funs, c_ptr);
    c_ptr += args->qp_solver->calculate_args_size(&qp_dims, qp_solver->qp_solver_funs);

    sim_dims sim_dims;

    args->sim_solvers = (sim_solver **) c_ptr;
    c_ptr += dims->N*sizeof(sim_solver *);

    args->sim_solvers_args = (void **) c_ptr;
    c_ptr += dims->N*sizeof(void *);

    int return_value, ns;

    for (int ii = 0; ii < dims->N; ii++)
    {
        cast_nlp_dims_to_sim_dims(&sim_dims, dims, ii);

        args->sim_solvers[ii] = (sim_solver *) c_ptr;
        c_ptr += sizeof(sim_solver);

        copy_module_pointers_to_args(args->sim_solvers[ii], &sim_solvers[ii]);

        args->sim_solvers_args[ii] = args->sim_solvers[ii]->assign_args(&sim_dims, c_ptr);
        c_ptr += args->sim_solvers[ii]->calculate_args_size(&sim_dims);
    }

    assert((char*)raw_memory + ocp_nlp_gn_sqp_calculate_args_size(dims, qp_solver, sim_solvers) == c_ptr);

    return args;
}



int ocp_nlp_gn_sqp_calculate_memory_size(ocp_nlp_dims *dims, ocp_nlp_gn_sqp_args *args)
{
    int size = 0;

    ocp_qp_dims qp_dims;
    cast_nlp_dims_to_qp_dims(&qp_dims, dims);

    size+= sizeof(ocp_nlp_gn_sqp_memory);

    int N = dims->N;

    size += sizeof(double *) * (N + 1);  // x
    size += sizeof(double *) * (N + 1);  // u
    size += sizeof(double *) * (N + 1);  // lam
    size += sizeof(double *) * N;  // pi

    for (int ii = 0; ii <= N; ii++)
    {
        size += sizeof(double)*dims->nx[ii];  // x
        size += sizeof(double)*dims->nu[ii];  // u
        size += sizeof(double)*2*(dims->nb[ii] + dims->ng[ii] + dims->nh[ii]);  // lam
        if (ii < N)
        {
            size += sizeof(double)*dims->nx[ii+1];  // pi
        }
    }

    size += args->qp_solver->calculate_memory_size(&qp_dims, args->qp_solver_args);

    sim_dims sim_dims;

    size += N*sizeof(void *);  // sim_solvers_mem

    for (int ii = 0; ii < N; ii++)
    {
        cast_nlp_dims_to_sim_dims(&sim_dims, dims, ii);
        size += args->sim_solvers[ii]->calculate_memory_size(&sim_dims, args->sim_solvers_args[ii]);
    }

    return size;
}



ocp_nlp_gn_sqp_memory *ocp_nlp_gn_sqp_assign_memory(ocp_nlp_dims *dims, ocp_nlp_gn_sqp_args *args, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_qp_dims qp_dims;
    cast_nlp_dims_to_qp_dims(&qp_dims, dims);

    ocp_nlp_gn_sqp_memory *mem = (ocp_nlp_gn_sqp_memory *)c_ptr;
    c_ptr += sizeof(ocp_nlp_gn_sqp_memory);

    assert((size_t)c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    int N = dims->N;

    mem->num_vars = number_of_primal_vars(dims);

    // double pointers
    assign_double_ptrs(N+1, &mem->x, &c_ptr);
    assign_double_ptrs(N+1, &mem->u, &c_ptr);
    assign_double_ptrs(N, &mem->pi, &c_ptr);
    assign_double_ptrs(N+1, &mem->lam, &c_ptr);

    // doubles
    assert((size_t)c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    for (int ii = 0; ii <= N; ii++)
    {
        assign_double(dims->nx[ii], &mem->x[ii], &c_ptr);
        assign_double(dims->nu[ii], &mem->u[ii], &c_ptr);
        if (ii < N)
        {
            assign_double(dims->nx[ii+1], &mem->pi[ii], &c_ptr);
        }
        assign_double(2*(dims->nb[ii] + dims->ng[ii] + dims->nh[ii]), &mem->lam[ii], &c_ptr);
    }

    assert((size_t)c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    // QP solver
    mem->qp_solver_mem = args->qp_solver->assign_memory(&qp_dims, args->qp_solver_args, c_ptr);
    c_ptr += args->qp_solver->calculate_memory_size(&qp_dims, args->qp_solver_args);

    // integrators
    sim_dims sim_dims;

    mem->sim_solvers_mem = (void **) c_ptr;
    c_ptr += N*sizeof(void *);

    for (int ii = 0; ii < N; ii++)
    {
        cast_nlp_dims_to_sim_dims(&sim_dims, dims, ii);
        mem->sim_solvers_mem[ii] = args->sim_solvers[ii]->assign_memory(&sim_dims, args->sim_solvers_args[ii], c_ptr);
        c_ptr += args->sim_solvers[ii]->calculate_memory_size(&sim_dims, args->sim_solvers_args[ii]);
    }

    mem->dims = dims;

    assert((char *)raw_memory + ocp_nlp_gn_sqp_calculate_memory_size(dims, args) == c_ptr);

    return mem;
}



int ocp_nlp_gn_sqp_calculate_workspace_size(ocp_nlp_dims *dims, ocp_nlp_gn_sqp_args *args)
{
    int size = 0;

    ocp_qp_dims qp_dims;
    cast_nlp_dims_to_qp_dims(&qp_dims, dims);

    size += sizeof(ocp_nlp_gn_sqp_work);

    size += number_of_primal_vars(dims) * sizeof(double);  // w

    size += ocp_qp_in_calculate_size(&qp_dims);
    size += ocp_qp_out_calculate_size(&qp_dims);
    size += args->qp_solver->calculate_workspace_size(&qp_dims, args->qp_solver_args);

    sim_dims sim_dims;

    size += dims->N*sizeof(sim_in *);
    size += dims->N*sizeof(sim_out *);
    size += dims->N*sizeof(void *);  // sim_work

    size += get_max_sim_workspace_size(dims, args);

    for (int ii = 0; ii < dims->N; ii++)
    {
        cast_nlp_dims_to_sim_dims(&sim_dims, dims, ii);
        size += sim_in_calculate_size(&sim_dims);
        size += sim_out_calculate_size(&sim_dims);
    }


    size += (dims->N+1)*sizeof(struct d_strvec);  // tmp_vecs

    for (int ii = 0; ii < dims->N+1; ii++)
    {
        size += d_size_strvec(qp_dims.nx[ii] + qp_dims.nu[ii]);
    }

    make_int_multiple_of(64, &size);
    size += 1 * 64;

    return size;
}



void ocp_nlp_gn_sqp_cast_workspace(ocp_nlp_gn_sqp_work *work, ocp_nlp_gn_sqp_memory *mem, ocp_nlp_gn_sqp_args *args)
{
    char *c_ptr = (char *)work;
    c_ptr += sizeof(ocp_nlp_gn_sqp_work);

    int N = mem->dims->N;
    int *nx = mem->dims->nx;
    int *nu = mem->dims->nu;

    ocp_qp_dims qp_dims;
    cast_nlp_dims_to_qp_dims(&qp_dims, mem->dims);

    // set up common nlp workspace
    assign_double(mem->num_vars, &work->w, &c_ptr);

    // set up local SQP data
    assign_strvec_ptrs(N+1, &work->tmp_vecs, &c_ptr);

    align_char_to(64, &c_ptr);

    for (int ii = 0; ii < N+1; ii++)
    {
        assign_strvec(nx[ii]+nu[ii], &work->tmp_vecs[ii], &c_ptr);
    }

    // set up QP solver
    work->qp_in = assign_ocp_qp_in(&qp_dims, c_ptr);
    c_ptr += ocp_qp_in_calculate_size(&qp_dims);
    work->qp_out = assign_ocp_qp_out(&qp_dims, c_ptr);
    c_ptr += ocp_qp_out_calculate_size(&qp_dims);

    work->qp_work = (void *)c_ptr;
    c_ptr += args->qp_solver->calculate_workspace_size(&qp_dims, args->qp_solver_args);

    // set up integrators
    sim_dims sim_dims;

    work->sim_in = (sim_in **) c_ptr;
    c_ptr += mem->dims->N*sizeof(sim_in *);
    work->sim_out = (sim_out **) c_ptr;
    c_ptr += mem->dims->N*sizeof(sim_out *);
    work->sim_solvers_work = (void **) c_ptr;
    c_ptr += mem->dims->N*sizeof(void *);

    int max_sim_work_size = get_max_sim_workspace_size(mem->dims, args);

    work->sim_solvers_work[0] = (void *)c_ptr;
    c_ptr += max_sim_work_size;

    for (int ii = 0; ii < mem->dims->N; ii++)
    {
        cast_nlp_dims_to_sim_dims(&sim_dims, mem->dims, ii);

        work->sim_in[ii] = assign_sim_in(&sim_dims, c_ptr);
        c_ptr += sim_in_calculate_size(&sim_dims);
        work->sim_out[ii] = assign_sim_out(&sim_dims, c_ptr);
        c_ptr += sim_out_calculate_size(&sim_dims);

        if (ii > 0) work->sim_solvers_work[ii] = work->sim_solvers_work[0];
    }

    assert((char *)work + ocp_nlp_gn_sqp_calculate_workspace_size(mem->dims, args) >= c_ptr);
}



static void initialize_objective(const ocp_nlp_in *nlp_in, ocp_nlp_gn_sqp_args *args, ocp_nlp_gn_sqp_memory *gn_sqp_mem, ocp_nlp_gn_sqp_work *work)
{
    int N = nlp_in->dims->N;
    int *nx = nlp_in->dims->nx;
    int *nu = nlp_in->dims->nu;
    ocp_nlp_ls_cost *cost = (ocp_nlp_ls_cost*) nlp_in->cost;

    // real_t **qp_Q = (real_t **) gn_sqp_mem->qp_solver->qp_in->Q;
    // real_t **qp_S = (real_t **) gn_sqp_mem->qp_solver->qp_in->S;
    // real_t **qp_R = (real_t **) gn_sqp_mem->qp_solver->qp_in->R;
	struct d_strmat *sRSQrq = work->qp_in->RSQrq;

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
    real_t *w = work->w;

    int_t w_idx = 0;
    for (int_t i = 0; i <= N; i++) {
        for (int_t j = 0; j < nx[i]; j++) {
            w[w_idx + j] = gn_sqp_mem->x[i][j];
        }
        for (int_t j = 0; j < nu[i]; j++) {
            w[w_idx + nx[i] + j] = gn_sqp_mem->u[i][j];
        }
        w_idx += nx[i] + nu[i];
    }
}



static void multiple_shooting(const ocp_nlp_in *nlp, ocp_nlp_gn_sqp_args *args, ocp_nlp_gn_sqp_memory *mem, ocp_nlp_gn_sqp_work *work) {

    int N = nlp->dims->N;
    int *nx = nlp->dims->nx;
    int *nu = nlp->dims->nu;
    int *nb = nlp->dims->nb;
    int *ng = nlp->dims->ng;

    real_t *w = work->w;
    struct d_strvec *stmp = work->tmp_vecs;

    ocp_nlp_ls_cost *cost = (ocp_nlp_ls_cost *) nlp->cost;
    real_t **y_ref = cost->y_ref;

    struct d_strmat *sBAbt = work->qp_in->BAbt;
    struct d_strvec *sb = work->qp_in->b;
    struct d_strmat *sRSQrq = work->qp_in->RSQrq;
    struct d_strvec *srq = work->qp_in->rq;
    struct d_strvec *sd = work->qp_in->d;

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
        for (int_t j = 0; j < nx[i]; j++) work->sim_in[i]->x[j] = w[w_idx+j];
        for (int_t j = 0; j < nu[i]; j++) work->sim_in[i]->u[j] = w[w_idx+nx[i]+j];
        args->sim_solvers[i]->fun(work->sim_in[i], work->sim_out[i], args->sim_solvers_args[i],
            mem->sim_solvers_mem[i], work->sim_solvers_work[i]);

        // TODO(rien): transition functions for changing dimensions not yet implemented!

        // convert b
        d_cvt_vec2strvec(nx[i+1], work->sim_out[i]->xn, &sb[i], 0);
        for (int j = 0; j < nx[i+1]; j++)
            DVECEL_LIBSTR(&sb[i], j) -= w[w_idx+nx[i]+nu[i]+j];
        // copy B
        d_cvt_tran_mat2strmat(nx[i+1], nu[i], &work->sim_out[i]->S_forw[nx[i+1]*nx[i]], nx[i+1], &sBAbt[i], 0, 0);
        // copy A
        d_cvt_tran_mat2strmat(nx[i+1], nx[i], &work->sim_out[i]->S_forw[0], nx[i+1], &sBAbt[i], nu[i], 0);
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
        sim_rk_opts *opts = (sim_rk_opts*) args->sim_solvers_args[i];

        for (int j = 0; j < nx[i]; j++)
            DVECEL_LIBSTR(&stmp[i], nu[i]+j) = w[w_idx+j]-y_ref[i][j];
        for (int j = 0; j < nu[i]; j++)
            DVECEL_LIBSTR(&stmp[i], j) = w[w_idx+nx[i]+j]-y_ref[i][nx[i]+j];

        dsymv_l_libstr(nu[i]+nx[i], nu[i]+nx[i], 1.0, &sRSQrq[i], 0, 0, &stmp[i], 0, 0.0, &srq[i], 0, &srq[i], 0);

        if (opts->scheme != NULL && opts->scheme->type != exact)
        {
            for (int_t j = 0; j < nx[i]; j++)
                DVECEL_LIBSTR(&srq[i], nu[i]+j) += work->sim_out[i]->grad[j];
            for (int_t j = 0; j < nu[i]; j++)
                DVECEL_LIBSTR(&srq[i], j) += work->sim_out[i]->grad[nx[i]+j];
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



static void update_variables(const ocp_nlp_in *nlp, ocp_nlp_gn_sqp_args *args, ocp_nlp_gn_sqp_memory *mem, ocp_nlp_gn_sqp_work *work) {
    int N = nlp->dims->N;
    int *nx = nlp->dims->nx;
    int *nu = nlp->dims->nu;


    for (int_t i = 0; i < N; i++)
    {
        for (int_t j = 0; j < nx[i+1]; j++)
        {
            work->sim_in[i]->S_adj[j] = -DVECEL_LIBSTR(&work->qp_out->pi[i], j);
            // sim[i].in->S_adj[j] = -mem->qp_solver->qp_out->pi[i][j];
        }
    }
    int_t w_idx = 0;
    for (int_t i = 0; i <= N; i++) {
        for (int_t j = 0; j < nx[i]; j++) {
            work->w[w_idx+j] += DVECEL_LIBSTR(&work->qp_out->ux[i], j + nu[i]);
            // w[w_idx+j] += mem->qp_solver->qp_out->x[i][j];
        }
        for (int_t j = 0; j < nu[i]; j++)
            work->w[w_idx+nx[i]+j] += DVECEL_LIBSTR(&work->qp_out->ux[i], j);
            // w[w_idx+nx[i]+j] += mem->qp_solver->qp_out->u[i][j];
        w_idx += nx[i]+nu[i];
    }
}



static void store_trajectories(const ocp_nlp_in *nlp, ocp_nlp_gn_sqp_memory *memory, ocp_nlp_out *out,
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

    sim_rk_opts *sim_opts;

    // set up integrators
    for (int ii = 0; ii < N; ii++)
    {
        sim_opts = args->sim_solvers_args[ii];
        work->sim_in[ii]->step = sim_opts->interval/sim_opts->num_steps;

        work->sim_in[ii]->vde = nlp_in->vde[ii];
        work->sim_in[ii]->jac = nlp_in->jac[ii];
        work->sim_in[ii]->vde_adj = nlp_in->vde_adj[ii];
        work->sim_in[ii]->forward_vde_wrapper = &vde_fun;
        work->sim_in[ii]->jacobian_wrapper = &jac_fun;
        work->sim_in[ii]->adjoint_vde_wrapper = &vde_hess_fun;

        // TODO(dimitris): REVISE IF THIS IS CORRECT FOR VARYING DIMENSIONS!
        for (int jj = 0; jj < nx[ii+1] * (nx[ii] + nu[ii]); jj++)
            work->sim_in[ii]->S_forw[jj] = 0.0;
        for (int jj = 0; jj < nx[ii+1]; jj++)
            work->sim_in[ii]->S_forw[jj * (nx[ii] + 1)] = 1.0;
        for (int jj = 0; jj < nx[ii] + nu[ii]; jj++)
            work->sim_in[ii]->S_adj[jj] = 0.0;
        // for (int jj = 0; jj < nlp_in->dims->num_stages[ii] * nx[ii+1]; jj++)
            // work->sim_in[ii]->grad_K[jj] = 0.0;
    }

    initialize_objective(nlp_in, args, mem, work);

    initialize_trajectories(nlp_in, mem, work);

    // TODO(dimitris): move somewhere else (not needed after new nlp_in)
    int_t **qp_idxb = (int_t **) work->qp_in->idxb;
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

    int_t max_sqp_iterations =  args->maxIter;

    acados_timer timer;
    real_t total_time = 0;
    acados_tic(&timer);
    for (int_t sqp_iter = 0; sqp_iter < max_sqp_iterations; sqp_iter++)
    {
        multiple_shooting(nlp_in, args, mem, work);

        int_t qp_status = args->qp_solver->fun(work->qp_in, work->qp_out,
            args->qp_solver_args, mem->qp_solver_mem, work->qp_work);

        if (qp_status != 0)
        {
            printf("QP solver returned error status %d\n", qp_status);
            return -1;
        }

        update_variables(nlp_in, args, mem, work);

        for (int_t i = 0; i < N; i++)
        {
            sim_rk_opts *opts = (sim_rk_opts*) args->sim_solvers_args[i];
            if (opts->scheme == NULL)
                continue;
            opts->sens_adj = (opts->scheme->type != exact);
            if (nlp_in->freezeSens) {
                // freeze inexact sensitivities after first SQP iteration !!
                opts->scheme->freeze = true;
            }
        }
    }

    total_time += acados_toc(&timer);
    store_trajectories(nlp_in, mem, nlp_out, work->w);

    return 0;
}
