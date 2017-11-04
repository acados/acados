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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

int_t ocp_nlp_gn_sqp_calculate_workspace_size(const ocp_nlp_in *in, void *args_) {
    ocp_nlp_gn_sqp_args *args = (ocp_nlp_gn_sqp_args*) args_;

    int_t size;

    size = sizeof(ocp_nlp_gn_sqp_work);
    size += ocp_nlp_calculate_workspace_size(in, args->common);
    return size;
}

static void ocp_nlp_gn_sqp_cast_workspace(ocp_nlp_gn_sqp_work *work,
                                          ocp_nlp_gn_sqp_memory *mem) {
    char *ptr = (char *)work;

    ptr += sizeof(ocp_nlp_gn_sqp_work);
    work->common = (ocp_nlp_work *)ptr;
    ocp_nlp_cast_workspace(work->common, mem->common);
}


static void initialize_objective(
    const ocp_nlp_in *nlp_in,
    ocp_nlp_gn_sqp_memory *gn_sqp_mem,
    ocp_nlp_gn_sqp_work *work) {

    const int_t N = nlp_in->N;
    const int_t *nx = nlp_in->nx;
    const int_t *nu = nlp_in->nu;
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

    const int_t N = nlp_in->N;
    const int_t *nx = nlp_in->nx;
    const int_t *nu = nlp_in->nu;
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


static void multiple_shooting(const ocp_nlp_in *nlp, ocp_nlp_gn_sqp_memory *mem, real_t *w) {

    const int_t N = nlp->N;
    const int_t *nx = nlp->nx;
    const int_t *nu = nlp->nu;
    const int_t *nb = nlp->nb;
    const int_t *ng = nlp->nc;

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

    // TODO(dimitris): TEMPORARY HACK TO STORE INTERM. RESULT
    int nvmax = 0;
    for (int i = 0; i < N+1; i++) {
        if (nx[i]+nu[i] > nvmax) {
            nvmax = nx[i]+nu[i];
        }
    }
    struct d_strvec stmp;
    d_allocate_strvec(nvmax, &stmp);

    int_t w_idx = 0;

    for (int_t i = 0; i < N; i++) {
        // Pass state and control to integrator
        for (int_t j = 0; j < nx[i]; j++) sim[i].in->x[j] = w[w_idx+j];
        for (int_t j = 0; j < nu[i]; j++) sim[i].in->u[j] = w[w_idx+nx[i]+j];
        sim[i].fun(sim[i].in, sim[i].out, sim[i].args, sim[i].mem, sim[i].work);

        // TODO(rien): transition functions for changing dimensions not yet implemented!

        d_cvt_vec2strvec(nx[i+1], &w[w_idx+nx[i]+nu[i]], &stmp, 0);

        // copy first part of b
        d_cvt_vec2strvec(nx[i+1], sim[i].out->xn, &sb[i], 0);
        // correct b
        daxpy_libstr(nx[i+1], -1.0, &stmp, 0, &sb[i], 0, &sb[i], 0);
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
        for (int_t j = 0; j < nlp->nb[i]; j++) {
// #ifdef FLIP_BOUNDS
            if (nlp->idxb[i][j] < nu[i]) {
                DVECEL_LIBSTR(&sd[i], j) = nlp->lb[i][j] - w[w_idx + nx[i] + nlp->idxb[i][j]];
                DVECEL_LIBSTR(&sd[i], j+nb[i]+ng[i]) = nlp->ub[i][j] - w[w_idx + nx[i] + nlp->idxb[i][j]];
                // qp_lb[i][j] = nlp->lb[i][j] - w[w_idx + nx[i] + nlp->idxb[i][j]];
                // qp_ub[i][j] = nlp->ub[i][j] - w[w_idx + nx[i] + nlp->idxb[i][j]];
            } else {
                DVECEL_LIBSTR(&sd[i], j) = nlp->lb[i][j] - w[w_idx - nu[i] + nlp->idxb[i][j]];
                DVECEL_LIBSTR(&sd[i], j+nb[i]+ng[i]) = nlp->ub[i][j] - w[w_idx - nu[i] + nlp->idxb[i][j]];
                // qp_lb[i][j] = nlp->lb[i][j] - w[w_idx - nu[i] + nlp->idxb[i][j]];
                // qp_ub[i][j] = nlp->ub[i][j] - w[w_idx - nu[i] + nlp->idxb[i][j]];
            }
// #else
//             qp_lb[i][j] = nlp->lb[i][j] - w[w_idx+nlp->idxb[i][j]];
//             qp_ub[i][j] = nlp->ub[i][j] - w[w_idx+nlp->idxb[i][j]];
// #endif
        }



        // Update gradients
        // TODO(rien): only for diagonal Q, R matrices atm
        // TODO(rien): only for least squares cost with state and control reference atm
        sim_RK_opts *opts = (sim_RK_opts*) sim[i].args;

        for (int j=0; j<nx[i]; j++)
            DVECEL_LIBSTR(&stmp, nu[i]+j) = w[w_idx+j]-y_ref[i][j];
        for (int j=0; j<nu[i]; j++)
            DVECEL_LIBSTR(&stmp, j) = w[w_idx+nx[i]+j]-y_ref[i][nx[i]+j];
        dsymv_l_libstr(nu[i]+nx[i], nu[i]+nx[i], 1.0, &sRSQrq[i], 0, 0, &stmp, 0, 0.0, &srq[i], 0, &srq[i], 0);
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

    for (int_t j = 0; j < nlp->nb[N]; j++) {
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
        DVECEL_LIBSTR(&stmp, j) = w[w_idx+j]-y_ref[N][j];
    dsymv_l_libstr(nx[N], nx[N], 1.0, &sRSQrq[N], 0, 0, &stmp, 0, 0.0, &srq[N], 0, &srq[N], 0);

    drowin_libstr(nx[N], 1.0, &srq[N], 0, &sRSQrq[N], nx[N], 0);

    // for (int_t j = 0; j < nx[N]; j++)
    //     qp_q[N][j] = cost->W[N][j*(nx[N]+nu[N]+1)]*(w[w_idx+j]-y_ref[N][j]);
}


static void update_variables(const ocp_nlp_in *nlp, ocp_nlp_gn_sqp_memory *mem, real_t *w) {
    const int_t N = nlp->N;
    const int_t *nx = nlp->nx;
    const int_t *nu = nlp->nu;
    sim_solver *sim = nlp->sim;

    #if 0
    for (int_t i = 0; i < N; i++)
        for (int_t j = 0; j < nx[i+1]; j++)
            sim[i].in->S_adj[j] = -mem->qp_solver->qp_out->pi[i][j];

    int_t w_idx = 0;
    for (int_t i = 0; i <= N; i++) {
        for (int_t j = 0; j < nx[i]; j++) {
            w[w_idx+j] += mem->qp_solver->qp_out->x[i][j];
        }
        for (int_t j = 0; j < nu[i]; j++)
            w[w_idx+nx[i]+j] += mem->qp_solver->qp_out->u[i][j];
        w_idx += nx[i]+nu[i];
    }
    #endif
}


static void store_trajectories(const ocp_nlp_in *nlp, ocp_nlp_memory *memory, ocp_nlp_out *out,
    real_t *w) {

    const int_t N = nlp->N;
    const int_t *nx = nlp->nx;
    const int_t *nu = nlp->nu;

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
int_t ocp_nlp_gn_sqp(const ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out, void *nlp_args_,
    void *nlp_mem_, void *nlp_work_) {

    ocp_nlp_gn_sqp_memory *gn_sqp_mem = (ocp_nlp_gn_sqp_memory *) nlp_mem_;
    ocp_nlp_gn_sqp_work *work = (ocp_nlp_gn_sqp_work*) nlp_work_;
    ocp_nlp_gn_sqp_cast_workspace(work, gn_sqp_mem);


    initialize_objective(nlp_in, gn_sqp_mem, work);

    initialize_trajectories(nlp_in, gn_sqp_mem, work);

    // TODO(roversch): Do we need this here?
    int_t **qp_idxb = (int_t **) gn_sqp_mem->qp_solver->qp_in->idxb;
    for (int_t i = 0; i <= nlp_in->N; i++) {
        for (int_t j = 0; j < nlp_in->nb[i]; j++) {
            if (nlp_in->idxb[i][j] < nlp_in->nx[i]) {
                // state bound
                qp_idxb[i][j] = nlp_in->idxb[i][j]+nlp_in->nu[i];
            } else {
                // input bound
                qp_idxb[i][j] = nlp_in->idxb[i][j]-nlp_in->nx[i];
            }
        }
    }

    int_t max_sqp_iterations = ((ocp_nlp_gn_sqp_args *) nlp_args_)->common->maxIter;

    acados_timer timer;
    real_t total_time = 0;
    acados_tic(&timer);
    for (int_t sqp_iter = 0; sqp_iter < max_sqp_iterations; sqp_iter++) {

        multiple_shooting(nlp_in, gn_sqp_mem, work->common->w);

        print_ocp_qp_dims(gn_sqp_mem->qp_solver->qp_in->size);
        print_ocp_qp_in(gn_sqp_mem->qp_solver->qp_in);

        int_t qp_status = gn_sqp_mem->qp_solver->fun(
            gn_sqp_mem->qp_solver->qp_in,
            gn_sqp_mem->qp_solver->qp_out,
            gn_sqp_mem->qp_solver->args,
            gn_sqp_mem->qp_solver->mem);

        printf("qp_status = %d\n", qp_status);
        exit(1);

        #if 0
        if (qp_status != 0) {
            printf("QP solver returned error status %d\n", qp_status);
            return -1;
        }

        update_variables(nlp_in, gn_sqp_mem, work->common->w);

        for (int_t i = 0; i < nlp_in->N; i++) {
            sim_RK_opts *opts = (sim_RK_opts*) nlp_in->sim[i].args;
            nlp_in->sim[i].in->sens_adj = (opts->scheme.type != exact);
            if (nlp_in->freezeSens) {  // freeze inexact sensitivities after first SQP iteration !!
                opts->scheme.freeze = true;
            }
        }
        #endif
    }

    total_time += acados_toc(&timer);
    store_trajectories(nlp_in, gn_sqp_mem->common, nlp_out, work->common->w);

    return 0;
}

void ocp_nlp_gn_sqp_create_memory(const ocp_nlp_in *in, void *args_, void *memory_) {

    ocp_nlp_gn_sqp_args *args = (ocp_nlp_gn_sqp_args *)args_;
    ocp_nlp_gn_sqp_memory *mem = (ocp_nlp_gn_sqp_memory *)memory_;

    // TODO(dimitris): CLEAN THIS UP ONCE NLP SIZE EXISTS!!!!
    int nbx[in->N+1];
    int nbu[in->N+1];
    int ns[in->N+1];
    for (int ii = 0; ii < in->N+1; ii++)
        ns[ii] = 0;

    form_nbx_nbu((int) in->N, nbx, nbu, (int*)in->nb, (int*)in->nx, (int*)in->nu, (int**)in->idxb);

    ocp_qp_dims dims;
    dims.N = (int) in->N;
    dims.nx = (int*)in->nx;
    dims.nu = (int*)in->nu;
    dims.nb = (int*)in->nb;
    dims.ng = (int*)in->nc;
    dims.ns = ns;
    dims.nbu = nbu;
    dims.nbx = nbx;

    ocp_qp_in *dummy_qp = create_ocp_qp_in(&dims);

    int_t **idxb = (int_t **) dummy_qp->idxb;
    for (int_t i = 0; i < in->N; i++)
        for (int_t j = 0; j < in->nb[i]; j++)
            idxb[i][j] = in->idxb[i][j];

    mem->qp_solver = create_ocp_qp_solver(dummy_qp, args->qp_solver_name, NULL);



    ocp_nlp_create_memory(in, mem->common);
}

void ocp_nlp_gn_sqp_free_memory(void *mem_) {
    ocp_nlp_gn_sqp_memory *mem = (ocp_nlp_gn_sqp_memory *)mem_;

    int_t N = mem->qp_solver->qp_in->size->N;
    ocp_nlp_free_memory(N, mem->common);

    mem->qp_solver->destroy(mem->qp_solver->mem, mem->qp_solver->work);

    free(mem->qp_solver->qp_in);
    free(mem->qp_solver->qp_out);
    free(mem->qp_solver->args);
    free(mem->qp_solver);
    // TODO(dimitris): where do we free the integrators?
}
