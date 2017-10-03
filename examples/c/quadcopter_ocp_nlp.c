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
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"

#include "acados/ocp_nlp/ocp_nlp_gn_sqp.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_hpmpc.h"
#include "acados/sim/sim_common.h"
#include "acados/sim/sim_erk_integrator.h"
#include "acados/sim/sim_lifted_irk_integrator.h"
#include "acados/sim/casadi_wrapper.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"

#include "examples/c/acados_gnuplot/acados_gnuplot.h"
#include "examples/c/quadcopter_model/ls_res_eval.h"

#define NN 50
#define TT 2.0
#define Ns 1
#define NX 11
#define NU 4
#define NSIM 2000
#define UMAX 3.0
#define NR 18
#define NR_END 14

#define MAX_SQP_ITERS 1
// #define N_SQP_HACK 100
#define N_SQP_HACK 1

// #define ALPHA 0.1
#define ALPHA 0.9


#define INITIAL_ANGLE_RAD 0.0
#define ANGLE_STEP 3.14/4.0
#define STEP_PERIOD 10.0
#define SIM_SCENARIO 1  // 0: stabilization, 1: tracking

#define MU_TIGHT 10
#define MM 2
#define LAM_INIT 1.0
#define T_INIT 0.1

#define OMEGA_REF 40.0

// #define PLOT_OL_RESULTS
#define PLOT_CL_RESULTS
// #define FP_EXCEPTIONS
#define PLOT_CONTROLS 0  // plot rates:0 plot controls:1

extern int ls_res_Fun(const real_t **arg, real_t **res, int *iw, real_t *w, int mem);
extern int ls_res_end_Fun(const real_t **arg, real_t **res, int *iw, real_t *w, int mem);
//
extern int vdeFun(const real_t **arg, real_t **res, int *iw, real_t *w, int mem);

int main() {
    const int d = 0;
    const int INEXACT = 0;

    int_t jj;

    // Problem data
    ocp_nlp_ls_cost ls_cost;

    // allocate memory for ls_cost
    ls_cost.ls_res = malloc(sizeof(char *)*(NN+1));
    ls_cost.ls_res_eval = malloc(sizeof(char *)*(NN+1));
    ls_cost.nr = (int_t *)malloc(sizeof(int_t)*(NN+1));

    // assign residual function pointers
    for (int_t i = 0; i < NN; i++) {
        ls_cost.ls_res[i] = &ls_res_Fun;
        ls_cost.ls_res_eval[i] = &ls_res_eval_quadcopter;
        ls_cost.nr[i] = NR;
    }

    ls_cost.ls_res[NN] = &ls_res_end_Fun;
    ls_cost.ls_res_eval[NN] = &ls_res_eval_end_quadcopter;
    ls_cost.nr[NN] = NR_END;

    real_t *W, *WN;
    real_t *uref;
    int_t max_sqp_iters = MAX_SQP_ITERS;
    real_t *x_end;
    real_t *u_end;

    d_zeros(&W, NX + NU, NX + NU);
    d_zeros(&WN, NX, NX);
    d_zeros(&uref, NU, 1);
    d_zeros(&x_end, NX, 1);
    d_zeros(&u_end, NU, 1);

    real_t initial_angle_rad = INITIAL_ANGLE_RAD;

    real_t x0[11] = {
        cos(initial_angle_rad/2),
        sin(initial_angle_rad/2),
        0.0000000000000000e+00,
        0.0000000000000000e+00,
        0.0000000000000000e+00,
        0.0000000000000000e+00,
        0.0000000000000000e+00,
        OMEGA_REF,
        OMEGA_REF,
        OMEGA_REF,
        OMEGA_REF,
    };

    real_t y_ref[18] =  {
        0.0000000000000000e+00,
        0.0000000000000000e+00,
        0.0000000000000000e+00,
        1.0000000000000000e+00,
        0.0000000000000000e+00,
        0.0000000000000000e+00,
        0.0000000000000000e+00,
        0.0000000000000000e+00,
        0.0000000000000000e+00,
        0.0000000000000000e+00,
        OMEGA_REF,
        OMEGA_REF,
        OMEGA_REF,
        OMEGA_REF,
        0.0000000000000000e+00,
        0.0000000000000000e+00,
        0.0000000000000000e+00,
        0.0000000000000000e+00,
    };

    real_t y_ref_end[14] =  {
        0.0000000000000000e+00,
        0.0000000000000000e+00,
        0.0000000000000000e+00,
        1.0000000000000000e+00,
        0.0000000000000000e+00,
        0.0000000000000000e+00,
        0.0000000000000000e+00,
        0.0000000000000000e+00,
        0.0000000000000000e+00,
        0.0000000000000000e+00,
        OMEGA_REF,
        OMEGA_REF,
        OMEGA_REF,
        OMEGA_REF,
    };

    real_t W_diag[15] = {
        1.0e+2,
        1.0e+2,
        1.0e+2,
        1.0e+2,
        1.0e-2,
        1.0e-2,
        1.0e-2,
        1.0e-2,
        1.0e-2,
        1.0e-2,
        1.0e-2,
        1.0e-2,
        1.0e-2,
        1.0e-2,
        1.0e-2,
    };

    for (int_t i = 0; i < NX; i++) W[i * (NX + NU + 1)] = W_diag[i];
    for (int_t i = 0; i < NU; i++) W[(NX + i) * (NX + NU + 1)] = W_diag[NX+i];
    for (int_t i = 0; i < NX; i++) WN[i * (NX + 1)] = W_diag[i];

    ls_cost.W = (real_t **)malloc(sizeof(*ls_cost.W) * (NN + 1));
    for (int_t i = 0; i < NN; i++) ls_cost.W[i] = W;
    ls_cost.W[NN] = WN;
    ls_cost.y_ref = (real_t **)malloc(sizeof(*ls_cost.y_ref) * (NN + 1));
    for (int_t i = 0; i < NN; i++) {
        ls_cost.y_ref[i] =
            (real_t *)malloc(sizeof(*ls_cost.y_ref[i]) * NR);
        for (int_t j = 0; j < NR; j++) ls_cost.y_ref[i][j] = y_ref[j];
    }
    ls_cost.y_ref[NN] = (real_t *)malloc(sizeof(*ls_cost.y_ref[NN]) * (NR_END));
    for (int_t j = 0; j < NR_END; j++) ls_cost.y_ref[NN][j] = y_ref_end[j];

    // Integrator structs
    real_t Ts = TT / NN;
    sim_in sim_in[NN];
    sim_out sim_out[NN];
    sim_info info[NN];
    sim_solver integrators[NN];

    sim_RK_opts rk_opts[NN];
    void *sim_work;
    sim_lifted_irk_memory irk_mem[NN];

    for (jj = 0; jj < NN; jj++) {
        integrators[jj].in = &sim_in[jj];
        integrators[jj].out = &sim_out[jj];
        integrators[jj].args = &rk_opts[jj];
        if (d > 0) {
            integrators[jj].fun = &sim_lifted_irk;
            integrators[jj].mem = &irk_mem[jj];
        } else {
            integrators[jj].fun = &sim_erk;
            integrators[jj].mem = 0;
        }

        sim_in[jj].num_steps = Ns;
        sim_in[jj].step = Ts / sim_in[jj].num_steps;
        sim_in[jj].nx = NX;
        sim_in[jj].nu = NU;

        sim_in[jj].sens_forw = true;
        sim_in[jj].sens_adj = false;
        sim_in[jj].sens_hess = false;
        sim_in[jj].num_forw_sens = NX + NU;

        sim_in[jj].vde = &vdeFun;
        sim_in[jj].VDE_forw = &vde_fun;

        sim_in[jj].x = (real_t *)malloc(sizeof(*sim_in[jj].x) * (NX));
        sim_in[jj].u = (real_t *)malloc(sizeof(*sim_in[jj].u) * (NU));
        sim_in[jj].S_forw =
            (real_t *)malloc(sizeof(*sim_in[jj].S_forw) * (NX * (NX + NU)));
        for (int_t i = 0; i < NX * (NX + NU); i++) sim_in[jj].S_forw[i] = 0.0;
        for (int_t i = 0; i < NX; i++) sim_in[jj].S_forw[i * (NX + 1)] = 1.0;

        sim_in[jj].S_adj =
            (real_t *)malloc(sizeof(*sim_in[jj].S_adj) * (NX + NU));
        for (int_t i = 0; i < NX + NU; i++) sim_in[jj].S_adj[i] = 0.0;

        sim_in[jj].grad_K =
            (real_t *)malloc(sizeof(*sim_in[jj].grad_K) * (d * NX));
        for (int_t i = 0; i < d * NX; i++) sim_in[jj].grad_K[i] = 0.0;

        sim_out[jj].xn = (real_t *)malloc(sizeof(*sim_out[jj].xn) * (NX));
        sim_out[jj].S_forw =
            (real_t *)malloc(sizeof(*sim_out[jj].S_forw) * (NX * (NX + NU)));
        sim_out[jj].info = &info[jj];
        sim_out[jj].grad =
            (real_t *)malloc(sizeof(*sim_out[jj].grad) * (NX + NU));

        int_t workspace_size;
        if (d > 0) {
            sim_irk_create_arguments(&rk_opts[jj], d, "Gauss");
            if (INEXACT == 0) {
                sim_irk_create_Newton_scheme(&rk_opts[jj], d, "Gauss", exact);
            } else if (INEXACT == 1 || INEXACT == 3) {
                sim_irk_create_Newton_scheme(&rk_opts[jj], d, "Gauss",
                                             simplified_in);
            } else if (INEXACT == 2 || INEXACT == 4) {
                sim_irk_create_Newton_scheme(&rk_opts[jj], d, "Gauss",
                                             simplified_inis);
            }

            workspace_size = sim_lifted_irk_calculate_workspace_size(
                &sim_in[jj], &rk_opts[jj]);
            sim_lifted_irk_create_memory(&sim_in[jj], &rk_opts[jj],
                                         &irk_mem[jj]);
        } else {
            sim_erk_create_arguments(&rk_opts[jj], 4);
            workspace_size =
                sim_erk_calculate_workspace_size(&sim_in[jj], &rk_opts[jj]);
        }
        // TODO(roversch): Next line is leaking memory!
        sim_work = (void *)malloc(workspace_size);
        integrators[jj].work = sim_work;
    }

    int_t nx[NN + 1] = {0};
    int_t nu[NN + 1] = {0};
    int_t nb[NN + 1] = {0};
    int_t nc[NN + 1] = {0};
    int_t ng[NN + 1] = {0};
    for (int_t i = 0; i < NN; i++) {
        nx[i] = NX;
        nu[i] = NU;
    }
    nx[NN] = NX;

    /************************************************
     * box constraints
     ************************************************/

    int *idxb0;
    int_zeros(&idxb0, NX + NU, 1);
    real_t *lb0;
    d_zeros(&lb0, NX + NU, 1);
    real_t *ub0;
    d_zeros(&ub0, NX + NU, 1);
#ifdef FLIP_BOUNDS
    for (jj = 0; jj < NU; jj++) {
        lb0[jj] = -UMAX;  // umin
        ub0[jj] = UMAX;   // umax
        idxb0[jj] = jj;
    }
    for (jj = 0; jj < NX; jj++) {
        lb0[NU+jj] = x0[jj];  // xmin
        ub0[NU+jj] = x0[jj];  // xmax
        idxb0[NU+jj] = NU+jj;
    }
#else
    for (jj = 0; jj < NX; jj++) {
        lb0[jj] = x0[jj];  // xmin
        ub0[jj] = x0[jj];  // xmax
        idxb0[jj] = jj;
    }
    for (jj = 0; jj < NU; jj++) {
        lb0[NX+jj] = -UMAX;  // umin
        ub0[NX+jj] = UMAX;   // umax
        idxb0[NX+jj] = NX+jj;
    }
#endif

    nb[0] = NX + NU;

    int *idxb1;
    int_zeros(&idxb1, NX + NU, 1);
    double *lb1[NN - 1];
    double *ub1[NN - 1];
    for (int_t i = 0; i < NN - 1; i++) {
        d_zeros(&lb1[i], NX + NU, 1);
        d_zeros(&ub1[i], NX + NU, 1);
#ifdef FLIP_BOUNDS
        for (jj = 0; jj < NU; jj++) {
            lb1[i][jj] = -UMAX;  // umin
            ub1[i][jj] = UMAX;   // umax
            idxb1[jj] = jj;
        }
        for (jj = 0; jj < NX; jj++) {
            lb1[i][NU+jj] = -1e12;
            ub1[i][NU+jj] = 1e12;
            idxb1[NU+jj] = NU + jj;
        }
#else
        for (jj = 0; jj < NX; jj++) {
            lb1[i][jj] = -1e12;
            ub1[i][jj] = 1e12;
            idxb1[jj] = jj;
        }
        for (jj = 0; jj < NU; jj++) {
            lb1[i][NX+jj] = -UMAX;  // umin
            ub1[i][NX+jj] = UMAX;   // umax
            idxb1[NX+jj] = NX+jj;
        }
#endif
        nb[i + 1] = NX + NU;
    }

    int *idxbN;
    int_zeros(&idxbN, NX, 1);
    real_t *lbN;
    d_zeros(&lbN, NX, 1);
    real_t *ubN;
    d_zeros(&ubN, NX, 1);
    for (jj = 0; jj < NX; jj++) {
        lbN[jj] = -1e12;  // xmin
        ubN[jj] = 1e12;  // xmax
        idxbN[jj] = jj;
    }
    nb[NN] = NX;

    real_t *hlb[NN + 1];
    real_t *hub[NN + 1];
    int *hidxb[NN + 1];

    hlb[0] = lb0;
    hub[0] = ub0;
    hidxb[0] = idxb0;
    for (int_t i = 1; i < NN; i++) {
        hlb[i] = lb1[i - 1];
        hub[i] = ub1[i - 1];
        hidxb[i] = idxb1;
    }
    hlb[NN] = lbN;
    hub[NN] = ubN;
    hidxb[NN] = idxbN;

    ocp_nlp_in nlp_in;
    nlp_in.N = NN;
    nlp_in.nx = nx;
    nlp_in.nu = nu;
    nlp_in.nb = nb;
    nlp_in.nc = nc;
    nlp_in.ng = ng;
    nlp_in.idxb = (const int_t **)hidxb;
    nlp_in.lb = (const real_t **)hlb;
    nlp_in.ub = (const real_t **)hub;
    nlp_in.sim = integrators;
    nlp_in.cost = &ls_cost;
    nlp_in.freezeSens = false;
    if (INEXACT > 2) nlp_in.freezeSens = true;

    ocp_nlp_out nlp_out;
    nlp_out.x = (real_t **)malloc(sizeof(*nlp_out.x) * (NN + 1));
    nlp_out.u = (real_t **)malloc(sizeof(*nlp_out.u) * NN);
    nlp_out.lam = (real_t **)malloc(sizeof(*nlp_out.lam) * NN);
    for (int_t i = 0; i < NN; i++) {
        nlp_out.x[i] = (real_t *)malloc(sizeof(*nlp_out.x[i]) * (NX));
        nlp_out.u[i] = (real_t *)malloc(sizeof(*nlp_out.u[i]) * (NU));
        nlp_out.lam[i] = (real_t *)malloc(sizeof(*nlp_out.lam[i]) * (NX));
    }
    nlp_out.x[NN] = (real_t *)malloc(sizeof(*nlp_out.x[NN]) * (NX));

    ocp_nlp_gn_sqp_args nlp_args;
    ocp_nlp_args nlp_common_args;
    nlp_args.lin_res = 0;
    nlp_args.common = &nlp_common_args;
    nlp_args.common->maxIter = max_sqp_iters;

    snprintf(nlp_args.qp_solver_name, sizeof(nlp_args.qp_solver_name), "%s",
             "hpmpc");  // supported: "condensing_qpoases", "ooqp", "qpdunes"

    ocp_nlp_gn_sqp_memory nlp_mem;
    ocp_nlp_memory nlp_mem_common;
    nlp_mem.common = &nlp_mem_common;
    nlp_mem.common->step_size = ALPHA;
    ocp_nlp_gn_sqp_create_memory(&nlp_in, &nlp_args, &nlp_mem);
    ocp_qp_hpmpc_args *qp_args = (ocp_qp_hpmpc_args *)nlp_mem.qp_solver->args;

    // TODO(ANDREA) UGLY HACK: how to move this into ocp_nlp_gn_sqp_create_memory?
    double *lam_in[NN+1];
    double *t_in[NN+1];
    double *pt[NN+1];

    for (int_t ii = 0; ii < NN; ii++) {
        d_zeros(&lam_in[ii], 2*nb[ii]+2*ng[ii], 1);
        d_zeros(&t_in[ii], 2*nb[ii]+2*ng[ii], 1);
        d_zeros(&pt[ii], 2*nb[ii]+2*ng[ii], 1);
    }

    d_zeros(&pt[NN], 2*nb[NN]+2*ng[NN], 1);
    d_zeros(&lam_in[NN], 2*nb[NN]+2*ng[NN], 1);
    d_zeros(&t_in[NN], 2*nb[NN]+2*ng[NN], 1);

    qp_args->lam0 = lam_in;
    qp_args->t0 = t_in;
    qp_args->sigma_mu = MU_TIGHT;
    qp_args->M = MM;

    ocp_qp_solver *qp_solver = (ocp_qp_solver *)nlp_mem.qp_solver;
    qp_solver->qp_out->t = pt;

    // initialize nlp dual variables
    for (int_t i = MM; i < NN; i++) {
        for (int_t j  = 0; j < 2*nb[i]; j++) {
            lam_in[i][j] = LAM_INIT;
            t_in[i][j] = T_INIT;
        }
    }

    for (int_t j  = 0; j < 2*NX; j++) {
        lam_in[NN][j] = LAM_INIT;
        t_in[NN][j] = T_INIT;
    }

    int_t work_space_size =
        ocp_nlp_gn_sqp_calculate_workspace_size(&nlp_in, &nlp_args, &nlp_mem);
    void *nlp_work = (void *)malloc(work_space_size);

    for (int_t i = 0; i < NN; i++) {
        for (int_t j = 0; j < NX; j++)
            nlp_mem.common->x[i][j] = x0[j];  // resX(j,i)
        for (int_t j = 0; j < NU; j++)
            nlp_mem.common->u[i][j] = 0.1;  // resU(j, i)
    }
    for (int_t j = 0; j < NX; j++)
        nlp_mem.common->x[NN][j] = x0[j];  // resX(j, NN)

    int_t status;

#ifdef PLOT_CL_RESULTS
    real_t w_cl[(NX+NU)*NSIM] = {0.0};
    for (int_t i = 0; i < NX; i++) w_cl[i] = x0[i];
    real_t T_SIM = TT/NN*NSIM;
#endif
    for (int_t sim_iter = 0; sim_iter < NSIM; sim_iter++) {
        if (SIM_SCENARIO == 1) {
            int_t step_samples = (int_t)(STEP_PERIOD*NN/TT);
            int_t ref_phase = (int_t)((sim_iter/step_samples)%2);
            switch (ref_phase) {
                case 0:
                    y_ref[0] = -ANGLE_STEP;
                    y_ref_end[0] = -ANGLE_STEP;
                    for (int_t i = 0; i < NN; i++) {
                        for (int_t j = 0; j < NR; j++) ls_cost.y_ref[i][j] = y_ref[j];
                    }
                    for (int_t j = 0; j < NR_END; j++) ls_cost.y_ref[NN][j] = y_ref_end[j];
                    break;
                case 1:
                    y_ref[0] = ANGLE_STEP;
                    y_ref_end[0] = ANGLE_STEP;
                    for (int_t i = 0; i < NN; i++) {
                        for (int_t j = 0; j < NR; j++) ls_cost.y_ref[i][j] = y_ref[j];
                    }
                    for (int_t j = 0; j < NR_END; j++) ls_cost.y_ref[NN][j] = y_ref_end[j];
                    break;
            }
        }
        for (int_t sqp_iter_hack = 0; sqp_iter_hack < N_SQP_HACK; sqp_iter_hack++) {
            status = ocp_nlp_gn_sqp(&nlp_in, &nlp_out, &nlp_args, &nlp_mem, nlp_work);

            // TODO(Andrea): UGLY HACK udpate of t and lam should take place inside
            // ocp_nlp_gn_sqp
            for (int_t i = MM; i < NN; i++) {
                for (int_t j  = 0; j < 2*nb[i]; j++) {
                    lam_in[i][j] = lam_in[i][j] +
                    ALPHA*(qp_solver->qp_out->lam[i][j] - lam_in[i][j]);
                    t_in[i][j] = qp_solver->qp_out->t[i][j] +
                    ALPHA*(qp_solver->qp_out->t[i][j] - t_in[i][j]);
                }
            }

            for (int_t j  = 0; j < 2*NX; j++) {
                lam_in[NN][j] = lam_in[NN][j] +
                ALPHA*(qp_solver->qp_out->lam[NN][j] - lam_in[NN][j]);
                t_in[NN][j] = qp_solver->qp_out->t[NN][j] +
                ALPHA*(qp_solver->qp_out->t[NN][j] - t_in[NN][j]);
            }
        }
        // forward simulation
        for (int_t j = 0; j < nx[0]; j++) integrators[0].in->x[j] = x0[j];
        for (int_t j = 0; j < nu[0]; j++) integrators[0].in->u[j] = nlp_out.u[0][j];
        integrators[0].fun(integrators[0].in, integrators[0].out,
            integrators[0].args, integrators[0].mem, integrators[0].work);

        // TODO(rien): transition functions for changing dimensions not yet implemented!
#ifdef PLOT_CL_RESULTS
        for (int_t j = 0; j < nx[0]; j++)
            w_cl[sim_iter*(nx[0]+nu[0]) + j] = x0[j];
        for (int_t j = 0; j < nu[0]; j++)
            w_cl[sim_iter*(nx[0]+nu[0]) + nx[0] + j] = nlp_out.u[0][j];

        for (int_t j = 0; j < nx[0]; j++)
            x0[j] = integrators[0].out->xn[j];
#endif

#ifdef FLIP_BOUNDS
        for (jj = 0; jj < nx[0]; jj++) {
            lb0[NU+jj] = integrators[0].out->xn[jj];  // xmin
            ub0[NU+jj] = integrators[0].out->xn[jj];  // xmax
        }
#else
        for (jj = 0; jj < nx[0]; jj++) {
            lb0[jj] = integrators[0].out->xn[jj];  // xmin
            ub0[jj] = integrators[0].out->xn[jj];  // xmax
        }
#endif

        printf("status = %i\n", status);
    }

#if defined(PLOT_OL_RESULTS) || defined(PLOT_CL_RESULTS)
#if defined (PLOT_OL_RESULTS)
        real_t gnu_plot_w[(NX+NU)*(NN+1)];
        for (int_t i = 0; i < NN; i++) {
            for (int_t j = 0; j < NX; j++)
                gnu_plot_w[i*(NX+NU)+j] = nlp_out.x[i][j];
            for (int_t j = 0; j < NU; j++)
                gnu_plot_w[i*(NX+NU)+NX+j] = nlp_out.u[i][j];
        }
        int_t N_plt = NN;
        real_t T_plt = TT;
#elif defined(PLOT_CL_RESULTS)
        real_t *gnu_plot_w = w_cl;
        int_t N_plt = NSIM;
        real_t T_plt = T_SIM;
#endif

        real_t *gnuplot_data[8];

        real_t q_1[N_plt];
        real_t q_2[N_plt];
        real_t q_3[N_plt];
        real_t q_4[N_plt];

#if PLOT_CONTROLS
        real_t w_1[N_plt];
        real_t w_2[N_plt];
        real_t w_3[N_plt];
        real_t w_4[N_plt];
#else
        real_t rw_1[N_plt];
        real_t rw_2[N_plt];
        real_t rw_3[N_plt];
        real_t rw_4[N_plt];
#endif

        for (int_t i = 0; i < N_plt; i++) {
            q_1[i] = gnu_plot_w[i*(NX+NU) + 0];
            q_2[i] = gnu_plot_w[i*(NX+NU) + 1];
            q_3[i] = gnu_plot_w[i*(NX+NU) + 2];
            q_4[i] = gnu_plot_w[i*(NX+NU) + 3];

#if PLOT_CONTROLS
            w_1[i] = gnu_plot_w[i*(NX+NU) + 7];
            w_2[i] = gnu_plot_w[i*(NX+NU) + 8];
            w_3[i] = gnu_plot_w[i*(NX+NU) + 9];
            w_4[i] = gnu_plot_w[i*(NX+NU) + 10];
#else
            rw_1[i] = gnu_plot_w[i*(NX+NU) + 11];
            rw_2[i] = gnu_plot_w[i*(NX+NU) + 12];
            rw_3[i] = gnu_plot_w[i*(NX+NU) + 13];
            rw_4[i] = gnu_plot_w[i*(NX+NU) + 14];
#endif
        }

        gnuplot_data[0] = q_1;
        gnuplot_data[1] = q_2;
        gnuplot_data[2] = q_3;
        gnuplot_data[3] = q_4;

#if PLOT_CONTROLS
        gnuplot_data[4] = w_1;
        gnuplot_data[5] = w_2;
        gnuplot_data[6] = w_3;
        gnuplot_data[7] = w_4;
#else
        gnuplot_data[4] = rw_1;
        gnuplot_data[5] = rw_2;
        gnuplot_data[6] = rw_3;
        gnuplot_data[7] = rw_4;
#endif

        char plot_1[10] = "q_1";
        char plot_2[10] = "q_2";
        char plot_3[10] = "q_3";
        char plot_4[10] = "q_4";
#if PLOT_CONTROLS
        char plot_5[10] = "w_1";
        char plot_6[10] = "w_2";
        char plot_7[10] = "w_3";
        char plot_8[10] = "w_4";
#else
        char plot_5[10]  = "rw_1";
        char plot_6[10] = "rw_2";
        char plot_7[10] = "rw_3";
        char plot_8[10] = "rw_4";
#endif
        char *labels[8];
        labels[0] = plot_1;
        labels[1] = plot_2;
        labels[2] = plot_3;
        labels[3] = plot_4;

        labels[4] = plot_5;
        labels[5] = plot_6;
        labels[6] = plot_7;
        labels[7] = plot_8;

        acados_gnuplot(gnuplot_data, 8, T_plt, N_plt, labels, 4, 2);
#endif  // PLOT_OL_RESULTS

    ocp_nlp_gn_sqp_free_memory(&nlp_mem);

#if 0
    for (int_t i = 0; i < NN; i++) {
        for (int_t j = 0; j < NX; j++)
            printf("%.3f  ", nlp_out.x[i][j]);
        for (int_t j = 0; j < NU; j++)
            printf("%.3f  ", nlp_out.u[i][j]);
        printf("\n");

    }
    for (int_t j = 0; j < NX; j++)
        printf("%.3f  ", nlp_out.x[NN][j]);
    printf("\n");
#endif

    d_free(W);
    d_free(WN);
    d_free(uref);
    d_free(x_end);
    d_free(u_end);

    int_free(idxb0);
    d_free(lb0);
    d_free(ub0);
    int_free(idxb1);
    for (jj = 0; jj < NN - 1; jj++) {
        d_free(lb1[jj]);
        d_free(ub1[jj]);
    }
    int_free(idxbN);
    d_free(lbN);
    d_free(ubN);

    for (jj = 0; jj < NN; jj++) {
        free(sim_in[jj].x);
        free(sim_in[jj].u);
        free(sim_in[jj].S_forw);
        free(sim_in[jj].S_adj);
        free(sim_in[jj].grad_K);
        free(sim_out[jj].xn);
        free(sim_out[jj].S_forw);
        free(sim_out[jj].grad);
        free(ls_cost.y_ref[jj]);
    }
    free(ls_cost.y_ref[NN]);
    free(ls_cost.y_ref);
    free(ls_cost.W);

    for (jj = 0; jj < NN; jj++) {
        free(nlp_out.x[jj]);
        free(nlp_out.u[jj]);
        free(nlp_out.lam[jj]);
    }
    free(nlp_out.x[NN]);
    free(nlp_out.x);
    free(nlp_out.u);
    free(nlp_out.lam);

    free(nlp_work);
}
