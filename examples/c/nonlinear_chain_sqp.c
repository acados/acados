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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"

#include "acados/ocp_qp/ocp_qp_condensing_qpoases.h"
#include "acados/sim/sim_erk_integrator.h"
#include "acados/sim/sim_lifted_irk_integrator.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"
#include "examples/c/chain_model/chain_model.h"

#define NN 10
#define UMAX 2
#define PARALLEL 0

extern int vde_chain_nm2(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
extern int vde_chain_nm3(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
extern int vde_chain_nm4(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
extern int vde_chain_nm5(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
extern int vde_chain_nm6(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
extern int vde_chain_nm7(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
extern int vde_chain_nm8(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
extern int vde_chain_nm9(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);

// static void shift_states(real_t *w, real_t *x_end, int_t NN) {
//    for (int_t i = 0; i < NN; i++) {
//        for (int_t j = 0; j < NX; j++) w[i*(NX+NU)+j] = w[(i+1)*(NX+NU)+j];
//    }
//    for (int_t j = 0; j < NX; j++) w[NN*(NX+NU)+j] = x_end[j];
// }
//
// static void shift_controls(real_t *w, real_t *u_end, int_t NN) {
//    for (int_t i = 0; i < NN-1; i++) {
//        for (int_t j = 0; j < NU; j++) w[i*(NX+NU)+NX+j] = w[(i+1)*(NX+NU)+NX+j];
//    }
//    for (int_t j = 0; j < NU; j++) w[(NN-1)*(NX+NU)+NX+j] = u_end[j];
// }

// Simple SQP example for acados
int main() {
    int_t nil;

    for (int_t implicit = 0; implicit < 4; implicit++) {
        if (implicit == 0) {
            printf("\n\n--------------------------------------------------------------------\n");
            printf("------------------ Explicit Runge-Kutta of order 4 -----------------\n");
            printf("--------------------------------------------------------------------\n");
        } else {
            printf("\n\n--------------------------------------------------------------------\n");
            printf("-------------- Lifted Implicit Runge-Kutta of order %d --------------\n",
                    2*implicit);
            printf("--------------------------------------------------------------------\n");
        }

    for (int_t NMF = 1; NMF < 9; NMF++) {
        printf("\n------------ NUMBER OF FREE MASSES =  %d ------------\n", NMF);
    int_t NX = 6*NMF;
    int_t NU = 3;
    int_t jj;

    // Problem data
    real_t  *x0;
    real_t  *w;  // States and controls stacked
    real_t  *Q;
    real_t  *R;
    real_t  *xref;
    real_t  *uref;
    int_t   max_sqp_iters       = 20;
    int_t   max_iters           = 1;
    real_t  *x_end;
    real_t  *u_end;
    FILE *initStates, *refStates;
    real_t feas, stepX, stepU;

    d_zeros(&x0, NX, 1);
    d_zeros(&xref, NX, 1);
    d_zeros(&uref, NU, 1);
    d_zeros(&x_end, NX, 1);
    d_zeros(&u_end, NU, 1);
    d_zeros(&w, NN*(NX+NU)+NX, 1);
    d_zeros(&Q, NX, NX);
    d_zeros(&R, NU, NU);

    switch (NMF) {
    case 1:
        initStates = fopen(X0_NM2_FILE, "r");
        break;
    case 2:
        initStates = fopen(X0_NM3_FILE, "r");
        break;
    case 3:
        initStates = fopen(X0_NM4_FILE, "r");
        break;
    case 4:
        initStates = fopen(X0_NM5_FILE, "r");
        break;
    case 5:
        initStates = fopen(X0_NM6_FILE, "r");
        break;
    case 6:
        initStates = fopen(X0_NM7_FILE, "r");
        break;
    case 7:
        initStates = fopen(X0_NM8_FILE, "r");
        break;
    default:
        initStates = fopen(X0_NM9_FILE, "r");
        break;
    }
    for (int_t i = 0; i < NX; i++) {
        nil = fscanf(initStates, "%lf", &x0[i]);
    }
    fclose(initStates);

    switch (NMF) {
    case 1:
        refStates = fopen(XN_NM2_FILE, "r");
        break;
    case 2:
        refStates = fopen(XN_NM3_FILE, "r");
        break;
    case 3:
        refStates = fopen(XN_NM4_FILE, "r");
        break;
    case 4:
        refStates = fopen(XN_NM5_FILE, "r");
        break;
    case 5:
        refStates = fopen(XN_NM6_FILE, "r");
        break;
    case 6:
        refStates = fopen(XN_NM7_FILE, "r");
        break;
    case 7:
        refStates = fopen(XN_NM8_FILE, "r");
        break;
    default:
        refStates = fopen(XN_NM9_FILE, "r");
        break;
    }
    for (int_t i = 0; i < NX; i++) {
        nil = fscanf(refStates, "%lf", &xref[i]);
    }
    fclose(refStates);

    for (int_t i = 0; i < NN; i++) {
        for (int_t j = 0; j < NX; j++) w[i*(NX+NU)+j] = xref[j];
    }
    for (int_t j = 0; j < NX; j++) w[NN*(NX+NU)+j] = xref[j];

    for (int_t i = 0; i < NX; i++) Q[i*(NX+1)] = 1.0;
    for (int_t i = 0; i < NU; i++) R[i*(NU+1)] = 0.01;

    // Integrator structs
    real_t T = 0.2;
    sim_in  sim_in[NN];
    sim_out sim_out[NN];
    sim_info info[NN];

    sim_RK_opts rk_opts[NN];
    void *sim_work = NULL;
    sim_lifted_irk_memory irk_mem[NN];

    for (jj = 0; jj < NN; jj++) {
        sim_in[jj].nSteps = 2;
        sim_in[jj].step = T/sim_in[jj].nSteps;
        sim_in[jj].nx = NX;
        sim_in[jj].nu = NU;

        sim_in[jj].sens_forw = true;
        sim_in[jj].sens_adj = false;
        sim_in[jj].sens_hess = false;
        sim_in[jj].nsens_forw = NX+NU;

        switch (NMF) {
        case 1:
            sim_in[jj].vde = &vde_chain_nm2;
            sim_in[jj].VDE_forw = &VDE_fun_nm2;
            sim_in[jj].jac_fun = &jac_fun_nm2;
            break;
        case 2:
            sim_in[jj].vde = &vde_chain_nm3;
            sim_in[jj].VDE_forw = &VDE_fun_nm3;
            sim_in[jj].jac_fun = &jac_fun_nm3;
            break;
        case 3:
            sim_in[jj].vde = &vde_chain_nm4;
            sim_in[jj].VDE_forw = &VDE_fun_nm4;
            sim_in[jj].jac_fun = &jac_fun_nm4;
            break;
        case 4:
            sim_in[jj].vde = &vde_chain_nm5;
            sim_in[jj].VDE_forw = &VDE_fun_nm5;
            sim_in[jj].jac_fun = &jac_fun_nm5;
            break;
        case 5:
            sim_in[jj].vde = &vde_chain_nm6;
            sim_in[jj].VDE_forw = &VDE_fun_nm6;
            sim_in[jj].jac_fun = &jac_fun_nm6;
            break;
        case 6:
            sim_in[jj].vde = &vde_chain_nm7;
            sim_in[jj].VDE_forw = &VDE_fun_nm7;
            sim_in[jj].jac_fun = &jac_fun_nm7;
            break;
        case 7:
            sim_in[jj].vde = &vde_chain_nm8;
            sim_in[jj].VDE_forw = &VDE_fun_nm8;
            sim_in[jj].jac_fun = &jac_fun_nm8;
            break;
        default:
            sim_in[jj].vde = &vde_chain_nm9;
            sim_in[jj].VDE_forw = &VDE_fun_nm9;
            sim_in[jj].jac_fun = &jac_fun_nm9;
            break;
        }

        sim_in[jj].x = malloc(sizeof(*sim_in[jj].x) * (NX));
        sim_in[jj].u = malloc(sizeof(*sim_in[jj].u) * (NU));
        sim_in[jj].S_forw = malloc(sizeof(*sim_in[jj].S_forw) * (NX*(NX+NU)));
        for (int_t i = 0; i < NX*(NX+NU); i++) sim_in[jj].S_forw[i] = 0.0;
        for (int_t i = 0; i < NX; i++) sim_in[jj].S_forw[i*(NX+1)] = 1.0;

        sim_out[jj].xn = malloc(sizeof(*sim_out[jj].xn) * (NX));
        sim_out[jj].S_forw = malloc(sizeof(*sim_out[jj].S_forw) * (NX*(NX+NU)));
        sim_out[jj].info = &info[jj];

        int_t workspace_size;
        if (implicit > 0) {
            sim_irk_create_arguments(&rk_opts[jj], implicit, "Gauss");

            workspace_size = sim_lifted_irk_calculate_workspace_size(&sim_in[jj], &rk_opts[jj]);
            sim_lifted_irk_create_memory(&sim_in[jj], &rk_opts[jj], &irk_mem[jj]);
        } else {
            sim_erk_create_arguments(&rk_opts[jj], 4);
            workspace_size = sim_erk_calculate_workspace_size(&sim_in[jj], &rk_opts[jj]);
        }
        if (jj == 0) sim_work = (void *) malloc(workspace_size);
    }

    int_t nx[NN+1] = {0};
    int_t nu[NN] = {0};
    int_t nb[NN+1] = {0};
    int_t nc[NN+1] = {0};
    for (int_t i = 0; i < NN; i++) {
        nx[i] = NX;
        nu[i] = NU;
    }
    nx[NN] = NX;

    /************************************************
    * box constraints
    ************************************************/

    int *idxb0;
    int_zeros(&idxb0, NX+NU, 1);
    real_t *lb0;
    d_zeros(&lb0, NX+NU, 1);
    real_t *ub0;
    d_zeros(&ub0, NX+NU, 1);
    for (jj = 0; jj < NX; jj++) {
        lb0[jj] = x0[jj];   //   xmin
        ub0[jj] = x0[jj];   //   xmax
        idxb0[jj] = jj;
    }
    for (; jj < NX+NU; jj++) {
        lb0[jj] = -UMAX;  //   umin
        ub0[jj] = UMAX;   //   umax
        idxb0[jj] = jj;
    }
    nb[0] = NX+NU;

    int *idxb1;
    int_zeros(&idxb1, NU, 1);
    double *lb1[NN-1];
    double *ub1[NN-1];
    for (int_t i = 0; i < NN-1; i++) {
    d_zeros(&lb1[i], NU, 1);
    d_zeros(&ub1[i], NU, 1);
//    for (jj = 0; jj < nbx; jj++) {
//        lb1[jj] = -4.0;  //   xmin
//        ub1[jj] = +4.0;   //   xmax
//        idxb1[jj] = jj;
//    }
    for (jj = 0; jj < NU; jj++) {
        lb1[i][jj] = -UMAX;  //   umin
        ub1[i][jj] = UMAX;   //   umax
        idxb1[jj] = NX+jj;
    }
    }

    real_t *pA[NN];
    real_t *pB[NN];
    real_t *pb[NN];
    real_t *pQ[NN+1];
    real_t *pS[NN];
    real_t *pR[NN];
    real_t *pq[NN+1];
    real_t *pr[NN];
    real_t *px[NN+1];
    real_t *pu[NN];
    real_t *ppi[NN];
    for (int_t i = 0; i < NN; i++) {
        d_zeros(&pA[i], nx[i+1], nx[i]);
        d_zeros(&pB[i], nx[i+1], nu[i]);
        d_zeros(&pb[i], nx[i+1], 1);
        d_zeros(&pS[i], nu[i], nx[i]);
        d_zeros(&pq[i], nx[i], 1);
        d_zeros(&pr[i], nu[i], 1);
        d_zeros(&px[i], nx[i], 1);
        d_zeros(&pu[i], nu[i], 1);
        d_zeros(&ppi[i], nx[i+1], 1);
    }
    d_zeros(&pq[NN], nx[NN], 1);
    d_zeros(&px[NN], nx[NN], 1);

    real_t *hlb[NN+1];
    real_t *hub[NN+1];
    int *hidxb[NN+1];

    hlb[0] = lb0;
    hub[0] = ub0;
    hidxb[0] = idxb0;
    for (int_t i = 1; i < NN; i++) {
        hlb[i] = lb1[i-1];
        hub[i] = ub1[i-1];
        hidxb[i] = idxb1;
        nb[i] = NU;
    }

    // Allocate OCP QP variables
    ocp_qp_in qp_in;
    qp_in.N = NN;
    ocp_qp_out qp_out;
    ocp_qp_condensing_qpoases_args args;
    real_t *work = NULL;
    qp_in.nx = nx;
    qp_in.nu = nu;
    qp_in.nb = nb;
    qp_in.nc = nc;
    for (int_t i = 0; i < NN; i++) {
        pQ[i] = Q;
        pR[i] = R;
    }
    pQ[NN] = Q;
    qp_in.Q = (const real_t **) pQ;
    qp_in.S = (const real_t **) pS;
    qp_in.R = (const real_t **) pR;
    qp_in.q = (const real_t **) pq;
    qp_in.r = (const real_t **) pr;
    qp_in.A = (const real_t **) pA;
    qp_in.B = (const real_t **) pB;
    qp_in.b = (const real_t **) pb;
    qp_in.lb = (const real_t **) hlb;
    qp_in.ub = (const real_t **) hub;
    qp_in.idxb = (const int **) hidxb;
    qp_out.x = px;
    qp_out.u = pu;
    qp_out.pi = ppi;

    acados_timer timer;
    real_t timings = 0;
    real_t timings_sim = 0;
    real_t timings_la = 0;
    real_t timings_ad = 0;
//    for (int_t iter = 0; iter < max_iters; iter++) {
//        printf("\n------ TIME STEP %d ------\n", iter);

    acados_tic(&timer);
    for (int_t sqp_iter = 0; sqp_iter < max_sqp_iters; sqp_iter++) {
        feas = -1e10; stepX = -1e10; stepU = -1e10;

#if PARALLEL
#pragma omp parallel for
#endif
        for (int_t i = 0; i < NN; i++) {
            // Pass state and control to integrator
            for (int_t j = 0; j < NX; j++) sim_in[i].x[j] = w[i*(NX+NU)+j];
            for (int_t j = 0; j < NU; j++) sim_in[i].u[j] = w[i*(NX+NU)+NX+j];
            if (implicit > 0) {
                sim_lifted_irk(&sim_in[i], &sim_out[i], &rk_opts[i], &irk_mem[i], sim_work);
            } else {
                sim_erk(&sim_in[i], &sim_out[i], &rk_opts[i], 0, sim_work);
            }

            for (int_t j = 0; j < NX; j++) {
                pb[i][j] = sim_out[i].xn[j] - w[(i+1)*(NX+NU)+j];
                if (fabs(pb[i][j]) > feas) feas = fabs(pb[i][j]);
                for (int_t k = 0; k < NX; k++)
                    pA[i][j*NX+k] = sim_out[i].S_forw[j*NX+k];  // COLUMN MAJOR FROM CASADI
            }
            for (int_t j = 0; j < NU; j++)
                for (int_t k = 0; k < NX; k++)
                    pB[i][j*NX+k] = sim_out[i].S_forw[(NX+j)*NX+k];  // COLUMN MAJOR FROM CASADI

            timings_sim += sim_out[i].info->CPUtime;
            timings_la += sim_out[i].info->LAtime;
            timings_ad += sim_out[i].info->ADtime;
        }
        for (int_t i = 0; i < NN; i++) {
            // Update bounds:
            if ( i == 0 ) {
                for (int_t j = 0; j < NU; j++) {
                    lb0[NX+j] = -UMAX - w[NX+j];
                    ub0[NX+j] = UMAX - w[NX+j];
                }
            } else {
                for (int_t j = 0; j < NU; j++) {
                    lb1[i-1][j] = -UMAX - w[i*(NX+NU)+NX+j];
                    ub1[i-1][j] = UMAX - w[i*(NX+NU)+NX+j];
                }
            }

            // Construct QP matrices
            for (int_t j = 0; j < NX; j++) {
                pq[i][j] = Q[j*(NX+1)]*(w[i*(NX+NU)+j]-xref[j]);
            }
            for (int_t j = 0; j < NU; j++) {
                pr[i][j] = R[j*(NU+1)]*(w[i*(NX+NU)+NX+j]-uref[j]);
            }
        }
        for (int_t j = 0; j < NX; j++) {
            lb0[j] = (x0[j]-w[j]);
        }
        for (int_t j = 0; j < NX; j++) {
            pq[NN][j] = Q[j*(NX+1)]*(w[NN*(NX+NU)+j]-xref[j]);
        }

        // Set updated bounds:
        hlb[0] = lb0;
        hub[0] = ub0;
        for (int_t i = 1; i < NN; i++) {
            hlb[i] = lb1[i-1];
            hub[i] = ub1[i-1];
        }

        int status = 0;
        status = ocp_qp_condensing_qpoases(&qp_in, &qp_out, &args, NULL, work);
        if (status) {
            printf("qpOASES returned error status %d\n", status);
            return -1;
        }
        for (int_t i = 0; i < NN; i++) {
            for (int_t j = 0; j < NX; j++) {
                w[i*(NX+NU)+j] += qp_out.x[i][j];
                if (fabs(qp_out.x[i][j]) > stepX) stepX = fabs(qp_out.x[i][j]);
            }
            for (int_t j = 0; j < NU; j++) {
                w[i*(NX+NU)+NX+j] += qp_out.u[i][j];
                if (fabs(qp_out.u[i][j]) > stepU) stepU = fabs(qp_out.u[i][j]);
            }
        }
        for (int_t j = 0; j < NX; j++) {
            w[NN*(NX+NU)+j] += qp_out.x[NN][j];
            if (fabs(qp_out.x[NN][j]) > stepX) stepX = fabs(qp_out.x[NN][j]);
        }

        if (sqp_iter == max_sqp_iters-1) {
            fprintf(stdout, "--- ITERATION %d, Infeasibility: %+.3e , step X: %+.3e, "
                    "step U: %+.3e \n", sqp_iter, feas, stepX, stepU);
        }
    }
//        for (int_t i = 0; i < NX; i++) x0[i] = w[NX+NU+i];
//        shift_states(w, x_end, NN);
//        shift_controls(w, u_end, NN);
        timings += acados_toc(&timer);
//    }

        printf("\nAverage of %.3f ms in the integrator,\n",
                1e3*timings_sim/(max_sqp_iters*max_iters));
        printf("  of which %.3f ms spent in CasADi and\n",
                1e3*timings_ad/(max_sqp_iters*max_iters));
        printf("  of which %.3f ms spent in BLASFEO.\n",
                1e3*timings_la/(max_sqp_iters*max_iters));
        printf("--Total of %.3f ms per SQP iteration.--\n",
                1e3*timings/(max_sqp_iters*max_iters));

//    #ifdef DEBUG
//    print_matrix_name("stdout", "sol", w, NX+NU, NN);
//    #endif  // DEBUG
    }
    }
    return 0*nil;
}
