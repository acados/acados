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

#include <iostream>
#include <string>

#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"

#include "catch/include/catch.hpp"

#include "acados/ocp_nlp/ocp_nlp_sm_gn.h"
#include "acados/ocp_nlp/ocp_nlp_sqp.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_qpdunes.h"
#include "acados/ocp_qp/ocp_qp_condensing_qpoases.h"
#include "acados/sim/sim_casadi_wrapper.h"
#include "acados/sim/sim_common.h"
#include "acados/sim/sim_erk_integrator.h"
#include "acados/sim/sim_lifted_irk_integrator.h"
#include "acados/utils/casadi_wrapper.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"

#include "test/ocp_nlp/chain/chain_model.h"
#include "test/ocp_nlp/chain/chain_ocp.h"
#include "test/test_utils/eigen.h"
#include "test/test_utils/read_matrix.h"

real_t COMPARISON_TOLERANCE_IPOPT = 1e-6;

#define NN 15
#define TT 3.0
#define Ns 2

using Eigen::MatrixXd;
using Eigen::VectorXd;

TEST_CASE("GN-SQP for nonlinear optimal control of chain of masses",
          "[nonlinear optimization]") {
    // TODO(nielsvd): re-implement (Frozen) IN/INIS
    for (int INEXACT = 0; INEXACT < 5; INEXACT++) {
        int d_start = 0;
        if (INEXACT > 0) d_start = 2;

        for (int d = d_start; d < 4; d++) {  // RK4 in case d == 0
            for (int NMF = 1; NMF < 4; NMF++) {
                if (INEXACT == 0) {
                    printf(
                        "\n----- NUMBER OF FREE MASSES = %d, d = %d (Exact "
                        "Newton) -----\n",
                        NMF, d);
                } else if (INEXACT == 1) {
                    printf(
                        "\n----- NUMBER OF FREE MASSES = %d, d = %d (IN "
                        "Scheme) -----\n",
                        NMF, d);
                } else if (INEXACT == 2) {
                    printf(
                        "\n----- NUMBER OF FREE MASSES = %d, d = %d (INIS "
                        "Scheme) -----\n",
                        NMF, d);
                } else if (INEXACT == 3) {
                    printf(
                        "\n----- NUMBER OF FREE MASSES = %d, d = %d (FROZEN IN "
                        "Scheme) -----\n",
                        NMF, d);
                } else if (INEXACT == 4) {
                    printf(
                        "\n----- NUMBER OF FREE MASSES = %d, d = %d (FROZEN "
                        "INIS Scheme) -----\n",
                        NMF, d);
                }
                int_t NX = 6 * NMF;
                int_t NU = 3;
                int_t jj;

                real_t wall_pos = -0.01;
                int_t UMAX = 10;

                // Problem data
                int_t N = NN;
                ocp_nlp_ls_cost ls_cost;
                real_t *W, *WN;
                real_t *uref;
                int_t max_sqp_iters = 20;
                real_t *x_end;
                real_t *u_end;

                /************************************************
                 * cost
                 ************************************************/
                d_zeros(&W, NX + NU, NX + NU);
                d_zeros(&WN, NX, NX);
                d_zeros(&uref, NU, 1);
                d_zeros(&x_end, NX, 1);
                d_zeros(&u_end, NU, 1);

                std::string NMFdat =
                    std::to_string(NMF + 1) + "_d" + std::to_string(d) + ".dat";
                VectorXd x0 =
                    readMatrix(std::string("ocp_nlp/chain/x0_nm") + NMFdat);
                VectorXd xref = readMatrix("ocp_nlp/chain/xN_nm" + NMFdat);

                MatrixXd resX = readMatrix("ocp_nlp/chain/resX_nm" + NMFdat);
                MatrixXd resU = readMatrix("ocp_nlp/chain/resU_nm" + NMFdat);

                for (int_t i = 0; i < NX; i++) W[i * (NX + NU + 1)] = 1e-2;
                for (int_t i = 0; i < NU; i++)
                    W[(NX + i) * (NX + NU + 1)] = 1.0;
                for (int_t i = 0; i < NX; i++) WN[i * (NX + 1)] = 1e-2;

                ls_cost.N = N;
                ls_cost.W = (real_t **)malloc(sizeof(*ls_cost.W) * (N + 1));
                for (int_t i = 0; i < N; i++) ls_cost.W[i] = W;
                ls_cost.W[N] = WN;
                ls_cost.y_ref =
                    (real_t **)malloc(sizeof(*ls_cost.y_ref) * (N + 1));
                ls_cost.fun =
                    (ocp_nlp_function **)malloc(sizeof(*ls_cost.fun) * (N + 1));
                for (int_t i = 0; i < N; i++) {
                    ls_cost.fun[i] =
                        (ocp_nlp_function *)malloc(sizeof(ocp_nlp_function));
                    // Initialize LS cost
                    ls_cost.fun[i]->nx = NX;
                    ls_cost.fun[i]->nu = NU;
                    ls_cost.fun[i]->np = 0;
                    ls_cost.fun[i]->ny = (NX + NU);
                    ls_cost.fun[i]->in =
                        (casadi_wrapper_in *)malloc(sizeof(casadi_wrapper_in));
                    ls_cost.fun[i]->in->compute_jac = true;
                    ls_cost.fun[i]->in->compute_hess = false;
                    ls_cost.fun[i]->out = (casadi_wrapper_out *)malloc(
                        sizeof(casadi_wrapper_out));
                    ls_cost.fun[i]->args = casadi_wrapper_create_arguments();
                    switch (NMF) {
                        case 1:
                            ls_cost.fun[i]->args->fun = &ls_cost_nm2;
                            ls_cost.fun[i]->args->dims = &ls_cost_nm2_work;
                            ls_cost.fun[i]->args->sparsity =
                                &ls_cost_nm2_sparsity_out;
                            break;
                        case 2:
                            ls_cost.fun[i]->args->fun = &ls_cost_nm3;
                            ls_cost.fun[i]->args->dims = &ls_cost_nm3_work;
                            ls_cost.fun[i]->args->sparsity =
                                &ls_cost_nm3_sparsity_out;
                            break;
                        case 3:
                            ls_cost.fun[i]->args->fun = &ls_cost_nm4;
                            ls_cost.fun[i]->args->dims = &ls_cost_nm4_work;
                            ls_cost.fun[i]->args->sparsity =
                                &ls_cost_nm4_sparsity_out;
                            break;
                        default:
                            REQUIRE(1 == 0);
                            break;
                    }
                    casadi_wrapper_initialize(ls_cost.fun[i]->in,
                                              ls_cost.fun[i]->args,
                                              &ls_cost.fun[i]->work);

                    ls_cost.y_ref[i] =
                        (real_t *)malloc(sizeof(*ls_cost.y_ref[i]) * (NX + NU));
                    for (int_t j = 0; j < NX; j++)
                        ls_cost.y_ref[i][j] = xref[j];
                    for (int_t j = 0; j < NU; j++)
                        ls_cost.y_ref[i][NX + j] = 0.0;
                }
                ls_cost.fun[N] =
                    (ocp_nlp_function *)malloc(sizeof(ocp_nlp_function));
                ls_cost.fun[N]->nx = NX;
                ls_cost.fun[N]->nu = 0;
                ls_cost.fun[N]->np = 0;
                ls_cost.fun[N]->ny = NX;
                ls_cost.fun[N]->in =
                    (casadi_wrapper_in *)malloc(sizeof(casadi_wrapper_in));
                ls_cost.fun[N]->in->compute_jac = true;
                ls_cost.fun[N]->in->compute_hess = false;
                ls_cost.fun[N]->out =
                    (casadi_wrapper_out *)malloc(sizeof(casadi_wrapper_out));
                ls_cost.fun[N]->args = casadi_wrapper_create_arguments();
                switch (NMF) {
                    case 1:
                        ls_cost.fun[N]->args->fun = &ls_costN_nm2;
                        ls_cost.fun[N]->args->dims = &ls_costN_nm2_work;
                        ls_cost.fun[N]->args->sparsity =
                            &ls_costN_nm2_sparsity_out;
                        break;
                    case 2:
                        ls_cost.fun[N]->args->fun = &ls_costN_nm3;
                        ls_cost.fun[N]->args->dims = &ls_costN_nm3_work;
                        ls_cost.fun[N]->args->sparsity =
                            &ls_costN_nm3_sparsity_out;
                        break;
                    case 3:
                        ls_cost.fun[N]->args->fun = &ls_costN_nm4;
                        ls_cost.fun[N]->args->dims = &ls_costN_nm4_work;
                        ls_cost.fun[N]->args->sparsity =
                            &ls_costN_nm4_sparsity_out;
                        break;
                    default:
                        REQUIRE(1 == 0);
                        break;
                }
                casadi_wrapper_initialize(ls_cost.fun[N]->in,
                                          ls_cost.fun[N]->args,
                                          &ls_cost.fun[N]->work);

                ls_cost.y_ref[N] =
                    (real_t *)malloc(sizeof(*ls_cost.y_ref[N]) * (NX));
                for (int_t j = 0; j < NX; j++) ls_cost.y_ref[N][j] = xref(j);

                /************************************************
                 * simulators
                 ************************************************/
                real_t Ts = TT / N;
                sim_in sim_in[N];
                sim_out sim_out[N];
                sim_info info[N];
                sim_solver *integrators[N];

                sim_rk_opts rk_opts[N];
                void *sim_work = NULL;
                sim_lifted_irk_memory irk_mem[N];

                // TODO(rien): can I move this somewhere inside the integrator?
                struct d_strmat str_mat[N];
                struct d_strmat str_sol[N];

                for (jj = 0; jj < N; jj++) {
                    integrators[jj] = (sim_solver *)malloc(sizeof(sim_solver));
                    integrators[jj]->in = &sim_in[jj];
                    integrators[jj]->out = &sim_out[jj];
                    integrators[jj]->args = &rk_opts[jj];
                    if (d > 0) {
                        integrators[jj]->fun = &sim_lifted_irk;
                        integrators[jj]->mem = &irk_mem[jj];
                    } else {
                        integrators[jj]->fun = &sim_erk;
                        integrators[jj]->mem = 0;
                    }

                    sim_in[jj].num_steps = Ns;
                    sim_in[jj].step = Ts / sim_in[jj].num_steps;
                    sim_in[jj].nx = NX;
                    sim_in[jj].nu = NU;

                    sim_in[jj].sens_forw = true;
                    sim_in[jj].sens_adj = false;
                    sim_in[jj].sens_hess = false;
                    sim_in[jj].num_forw_sens = NX + NU;

                    switch (NMF) {
                        case 1:
                            sim_in[jj].vde = &vde_chain_nm2;
                            sim_in[jj].forward_vde_wrapper = &vde_fun;
                            sim_in[jj].jac = &jac_chain_nm2;
                            sim_in[jj].jacobian_wrapper = &jac_fun;
                            break;
                        case 2:
                            sim_in[jj].vde = &vde_chain_nm3;
                            sim_in[jj].forward_vde_wrapper = &vde_fun;
                            sim_in[jj].jac = &jac_chain_nm3;
                            sim_in[jj].jacobian_wrapper = &jac_fun;
                            break;
                        case 3:
                            sim_in[jj].vde = &vde_chain_nm4;
                            sim_in[jj].forward_vde_wrapper = &vde_fun;
                            sim_in[jj].jac = &jac_chain_nm4;
                            sim_in[jj].jacobian_wrapper = &jac_fun;
                            break;
                        default:
                            REQUIRE(1 == 0);
                            break;
                    }

                    sim_in[jj].x =
                        (real_t *)malloc(sizeof(*sim_in[jj].x) * (NX));
                    sim_in[jj].u =
                        (real_t *)malloc(sizeof(*sim_in[jj].u) * (NU));
                    sim_in[jj].S_forw = (real_t *)malloc(
                        sizeof(*sim_in[jj].S_forw) * (NX * (NX + NU)));
                    for (int_t i = 0; i < NX * (NX + NU); i++)
                        sim_in[jj].S_forw[i] = 0.0;
                    for (int_t i = 0; i < NX; i++)
                        sim_in[jj].S_forw[i * (NX + 1)] = 1.0;

                    sim_in[jj].S_adj =
                        (real_t *)malloc(sizeof(*sim_in[jj].S_adj) * (NX + NU));
                    for (int_t i = 0; i < NX + NU; i++)
                        sim_in[jj].S_adj[i] = 0.0;

                    sim_in[jj].grad_K =
                        (real_t *)malloc(sizeof(*sim_in[jj].grad_K) * (d * NX));
                    for (int_t i = 0; i < d * NX; i++)
                        sim_in[jj].grad_K[i] = 0.0;

                    sim_out[jj].xn =
                        (real_t *)malloc(sizeof(*sim_out[jj].xn) * (NX));
                    sim_out[jj].S_forw = (real_t *)malloc(
                        sizeof(*sim_out[jj].S_forw) * (NX * (NX + NU)));
                    sim_out[jj].info = &info[jj];
                    sim_out[jj].grad =
                        (real_t *)malloc(sizeof(*sim_out[jj].grad) * (NX + NU));

                    int_t workspace_size;
                    if (d > 0) {
                        sim_irk_create_arguments(&rk_opts[jj], d, "Gauss");
                        if (INEXACT == 0) {
                            sim_irk_create_Newton_scheme(&rk_opts[jj], d,
                                                         "Gauss", exact);
                        } else if (INEXACT == 1 || INEXACT == 3) {
                            sim_irk_create_Newton_scheme(
                                &rk_opts[jj], d, "Gauss", simplified_in);
                        } else if (INEXACT == 2 || INEXACT == 4) {
                            sim_irk_create_Newton_scheme(
                                &rk_opts[jj], d, "Gauss", simplified_inis);
                        }

                        workspace_size =
                            sim_lifted_irk_calculate_workspace_size(
                                &sim_in[jj], &rk_opts[jj]);
                        sim_lifted_irk_create_memory(&sim_in[jj], &rk_opts[jj],
                                                     &irk_mem[jj]);
                    } else {
                        sim_erk_create_arguments(&rk_opts[jj], 4);
                        workspace_size = sim_erk_calculate_workspace_size(
                            &sim_in[jj], &rk_opts[jj]);
                    }
                    if (jj == 0) sim_work = (void *)malloc(workspace_size);
                    integrators[jj]->work = sim_work;
                }

                int_t nx[NN + 1] = {0};
                int_t nu[NN + 1] = {0};
                int_t nb[NN + 1] = {0};
                int_t nc[NN + 1] = {0};
                int_t ng[NN + 1] = {0};
                for (int_t i = 0; i < N; i++) {
                    nx[i] = NX;
                    nu[i] = NU;
                }
                nx[N] = NX;
                nu[N] = 0;

                /************************************************
                 * box constraints
                 ************************************************/

                int *idxb0;
                int_zeros(&idxb0, NX + NU, 1);
                real_t *lb0;
                d_zeros(&lb0, NX + NU, 1);
                real_t *ub0;
                d_zeros(&ub0, NX + NU, 1);
                for (jj = 0; jj < NX; jj++) {
                    lb0[jj] = x0(jj);  // xmin
                    ub0[jj] = x0(jj);  // xmax
                    idxb0[jj] = jj;
                }
                for (; jj < NX + NU; jj++) {
                    lb0[jj] = -UMAX;  // umin
                    ub0[jj] = UMAX;   // umax
                    idxb0[jj] = jj;
                }
                nb[0] = NX + NU;

                int *idxb1;
                int_zeros(&idxb1, NMF + NU, 1);
                double *lb1[N - 1];
                double *ub1[N - 1];
                for (int_t i = 0; i < N - 1; i++) {
                    d_zeros(&lb1[i], NMF + NU, 1);
                    d_zeros(&ub1[i], NMF + NU, 1);
                    for (jj = 0; jj < NMF; jj++) {
                        lb1[i][jj] = wall_pos;  // wall position
                        ub1[i][jj] = 1e12;
                        idxb1[jj] = 6 * jj + 1;
                    }
                    for (jj = 0; jj < NU; jj++) {
                        lb1[i][NMF + jj] = -UMAX;  // umin
                        ub1[i][NMF + jj] = UMAX;   // umax
                        idxb1[NMF + jj] = NX + jj;
                    }
                    nb[i + 1] = NMF + NU;
                }

                int *idxbN;
                int_zeros(&idxbN, NX, 1);
                real_t *lbN;
                d_zeros(&lbN, NX, 1);
                real_t *ubN;
                d_zeros(&ubN, NX, 1);
                for (jj = 0; jj < NX; jj++) {
                    lbN[jj] = xref(jj);  // xmin
                    ubN[jj] = xref(jj);  // xmax
                    idxbN[jj] = jj;
                }
                nb[N] = NX;

                real_t *hlb[N + 1];
                real_t *hub[N + 1];
                int *hidxb[N + 1];

                hlb[0] = lb0;
                hub[0] = ub0;
                hidxb[0] = idxb0;
                for (int_t i = 1; i < N; i++) {
                    hlb[i] = lb1[i - 1];
                    hub[i] = ub1[i - 1];
                    hidxb[i] = idxb1;
                }
                hlb[N] = lbN;
                hub[N] = ubN;
                hidxb[N] = idxbN;

                /************************************************
                 * nonlinear path constraints
                 ************************************************/
                ocp_nlp_function **path_constraints =
                    (ocp_nlp_function **)malloc(sizeof(ocp_nlp_function *) *
                                                (N + 1));
                for (int_t i = 0; i < N; i++) {
                    // Initialize path constraints
                    path_constraints[i] =
                        (ocp_nlp_function *)malloc(sizeof(ocp_nlp_function));
                    path_constraints[i]->nx = NX;
                    path_constraints[i]->nu = NU;
                    path_constraints[i]->np = 0;
                    path_constraints[i]->ny = (NX + NU);
                    path_constraints[i]->in =
                        (casadi_wrapper_in *)malloc(sizeof(casadi_wrapper_in));
                    path_constraints[i]->in->compute_jac = true;
                    path_constraints[i]->in->compute_hess = false;
                    path_constraints[i]->out = (casadi_wrapper_out *)malloc(
                        sizeof(casadi_wrapper_out));
                    path_constraints[i]->args =
                        casadi_wrapper_create_arguments();
                    switch (NMF) {
                        case 1:
                            path_constraints[i]->args->fun = &pathcon_nm2;
                            path_constraints[i]->args->dims = &pathcon_nm2_work;
                            path_constraints[i]->args->sparsity =
                                &pathcon_nm2_sparsity_out;
                            break;
                        case 2:
                            path_constraints[i]->args->fun = &pathcon_nm3;
                            path_constraints[i]->args->dims = &pathcon_nm3_work;
                            path_constraints[i]->args->sparsity =
                                &pathcon_nm3_sparsity_out;
                            break;
                        case 3:
                            path_constraints[i]->args->fun = &pathcon_nm4;
                            path_constraints[i]->args->dims = &pathcon_nm4_work;
                            path_constraints[i]->args->sparsity =
                                &pathcon_nm4_sparsity_out;
                            break;
                        default:
                            REQUIRE(1 == 0);
                            break;
                    }
                    casadi_wrapper_initialize(path_constraints[i]->in,
                                              path_constraints[i]->args,
                                              &path_constraints[i]->work);
                }
                path_constraints[N] =
                    (ocp_nlp_function *)malloc(sizeof(ocp_nlp_function));
                path_constraints[N]->nx = NX;
                path_constraints[N]->nu = 0;
                path_constraints[N]->np = 0;
                path_constraints[N]->ny = NX;
                path_constraints[N]->in =
                    (casadi_wrapper_in *)malloc(sizeof(casadi_wrapper_in));
                path_constraints[N]->in->compute_jac = true;
                path_constraints[N]->in->compute_hess = false;
                path_constraints[N]->out =
                    (casadi_wrapper_out *)malloc(sizeof(casadi_wrapper_out));
                path_constraints[N]->args = casadi_wrapper_create_arguments();
                switch (NMF) {
                    case 1:
                        path_constraints[N]->args->fun = &pathconN_nm2;
                        path_constraints[N]->args->dims = &pathconN_nm2_work;
                        path_constraints[N]->args->sparsity =
                            &pathconN_nm2_sparsity_out;
                        break;
                    case 2:
                        path_constraints[N]->args->fun = &pathconN_nm3;
                        path_constraints[N]->args->dims = &pathconN_nm3_work;
                        path_constraints[N]->args->sparsity =
                            &pathconN_nm3_sparsity_out;
                        break;
                    case 3:
                        path_constraints[N]->args->fun = &pathconN_nm4;
                        path_constraints[N]->args->dims = &pathconN_nm4_work;
                        path_constraints[N]->args->sparsity =
                            &pathconN_nm4_sparsity_out;
                        break;
                    default:
                        REQUIRE(1 == 0);
                        break;
                }
                casadi_wrapper_initialize(path_constraints[N]->in,
                                          path_constraints[N]->args,
                                          &path_constraints[N]->work);

                /************************************************
                 * sensitivity method
                 ************************************************/
                ocp_nlp_sm sensitivity_method;
                sensitivity_method.fun = &ocp_nlp_sm_gn;
                sensitivity_method.initialize = &ocp_nlp_sm_gn_initialize;
                sensitivity_method.destroy = &ocp_nlp_sm_gn_destroy;
                sensitivity_method.args = ocp_nlp_sm_gn_create_arguments();
                if (INEXACT > 2) {
                    ((ocp_nlp_sm_gn_args *)sensitivity_method.args)
                        ->freezeSens = true;
                }

                /************************************************
                 * QP solver
                 ************************************************/
                ocp_qp_solver qp_solver;
                qp_solver.fun = &ocp_qp_qpdunes;
                qp_solver.initialize = &ocp_qp_qpdunes_initialize;
                qp_solver.destroy = &ocp_qp_qpdunes_destroy;
                qp_solver.qp_in = create_ocp_qp_in(N, nx, nu, nb, ng);
                qp_solver.qp_out = create_ocp_qp_out(N, nx, nu, nb, ng);
                // TODO(nielsvd): lines below should go
                int_t **idxb = (int_t **) qp_solver.qp_in->idxb;
                for (int_t i = 0; i <= N; i++)
                    for (int_t j = 0; j < nb[i]; j++) idxb[i][j] = hidxb[i][j];
                qp_solver.args = (void *)ocp_qp_qpdunes_create_arguments(
                    QPDUNES_NONLINEAR_MPC);  // qp_solver.qp_in); //

                /************************************************
                 * SQP method
                 ************************************************/

                ocp_nlp_in nlp_in;
                nlp_in.N = N;
                nlp_in.nx = nx;
                nlp_in.nu = nu;
                nlp_in.nb = nb;
                nlp_in.ng = ng;
                nlp_in.idxb = (const int_t **)hidxb;
                nlp_in.lb = (const real_t **)hlb;
                nlp_in.ub = (const real_t **)hub;
                nlp_in.lg = NULL;
                nlp_in.ug = NULL;
                nlp_in.sim = (void **)&integrators;
                nlp_in.cost = (void *)&ls_cost;
                nlp_in.path_constraints = (void **)path_constraints;

                ocp_nlp_out nlp_out;
                nlp_out.x = (real_t **)malloc(sizeof(*nlp_out.x) * (N + 1));
                nlp_out.u = (real_t **)malloc(sizeof(*nlp_out.u) * (N + 1));
                nlp_out.pi = (real_t **)malloc(sizeof(*nlp_out.pi) * (N + 1));
                nlp_out.lam = (real_t **)malloc(sizeof(*nlp_out.lam) * (N + 1));
                // Allocate output variables
                for (int_t i = 0; i < N; i++) {
                    nlp_out.x[i] =
                        (real_t *)malloc(sizeof(*nlp_out.x[i]) * (NX));
                    nlp_out.u[i] =
                        (real_t *)malloc(sizeof(*nlp_out.u[i]) * (NU));
                    nlp_out.pi[i] =
                        (real_t *)malloc(sizeof(*nlp_out.pi[i]) * (NX));
                    nlp_out.lam[i] = (real_t *)malloc(
                        sizeof(*nlp_out.lam[i]) * 2 * nb[i] + 2 * ng[i]);
                }
                nlp_out.x[N] = (real_t *)malloc(sizeof(*nlp_out.x[N]) * (NX));
                nlp_out.u[N] = (real_t *)malloc(sizeof(*nlp_out.u[N]) * 0);
                nlp_out.pi[N] = (real_t *)malloc(sizeof(*nlp_out.pi[N]) * 0);
                nlp_out.lam[N] = (real_t *)malloc(
                    sizeof(*nlp_out.lam[N]) * 2 * nb[N] + 2 * ng[N]);

                ocp_nlp_sqp_args *nlp_args = ocp_nlp_sqp_create_arguments();
                nlp_args->maxIter = max_sqp_iters;
                nlp_args->sensitivity_method = &sensitivity_method;
                nlp_args->qp_solver = &qp_solver;

                ocp_nlp_sqp_memory *nlp_mem;
                ocp_nlp_sqp_workspace *nlp_work;
                ocp_nlp_sqp_initialize(&nlp_in, nlp_args, (void **)&nlp_mem,
                                       (void **)&nlp_work);

                // TOOD(nielsvd): should go, old interface


                // TODO(nielsvd): set memory to zero during allocation
                real_t **nlp_x_mem = (real_t **)nlp_mem->common->x;
                real_t **nlp_u_mem = (real_t **)nlp_mem->common->u;
                for (int_t i = 0; i < N; i++) {
                    for (int_t j = 0; j < NX; j++)
                        nlp_x_mem[i][j] = xref[j];  // resX(j,i)
                    for (int_t j = 0; j < NU; j++)
                        nlp_u_mem[i][j] = 0.0;  // resU(j, i)
                }
                for (int_t j = 0; j < NX; j++)
                    nlp_x_mem[N][j] = xref[j];  // resX(j, N)

                int_t status;

                status =
                    ocp_nlp_sqp(&nlp_in, &nlp_out, nlp_args, nlp_mem, nlp_work);
                REQUIRE(status == 0);

                real_t out_x[NX * (N + 1)], err_x[NX * (N + 1)];
                real_t out_u[NU * N], err_u[NU * N];
                for (int_t i = 0; i < N; i++) {
                    for (int_t j = 0; j < NX; j++)
                        out_x[i * NX + j] = nlp_out.x[i][j];
                    for (int_t j = 0; j < NU; j++)
                        out_u[i * NU + j] = nlp_out.u[i][j];
                }
                for (int_t j = 0; j < NX; j++)
                    out_x[N * NX + j] = nlp_out.x[N][j];

                for (int_t i = 0; i < N; i++) {
                    for (int_t j = 0; j < NX; j++)
                        err_x[i * NX + j] =
                            fabs(out_x[i * NX + j] - resX(j, i));
                    for (int_t j = 0; j < NU; j++)
                        err_u[i * NU + j] =
                            fabs(out_u[i * NU + j] - resU(j, i));
                }
                for (int_t j = 0; j < NX; j++)
                    err_x[N * NX + j] = fabs(out_x[N * NX + j] - resX(j, N));

                // print_matrix_name((char*)"stdout", (char*)"out_x", out_x, NX,
                // N+1);
                // print_matrix_name((char*)"stdout", (char*)"out_u", out_u, NU,
                // N);

                print_matrix_name((char *)"stdout", (char *)"err_x", err_x, NX,
                                  N + 1);
                print_matrix_name((char *)"stdout", (char *)"err_u", err_u, NU,
                                  N);

                std::cout << resX << std::endl;
                std::cout << resU << std::endl;

                MatrixXd SQP_x = Eigen::Map<MatrixXd>(&out_x[0], NX, N + 1);
                MatrixXd SQP_u = Eigen::Map<MatrixXd>(&out_u[0], NU, N);

                std::cout << "SQP_x:" << std::endl;
                std::cout << SQP_x << std::endl;

                REQUIRE(SQP_x.isApprox(resX, COMPARISON_TOLERANCE_IPOPT));
                REQUIRE(SQP_u.isApprox(resU, COMPARISON_TOLERANCE_IPOPT));

                d_free(W);
                d_free(WN);
                d_free(uref);
                d_free(x_end);
                d_free(u_end);

                int_free(idxb0);
                d_free(lb0);
                d_free(ub0);
                int_free(idxb1);
                for (jj = 0; jj < N - 1; jj++) {
                    d_free(lb1[jj]);
                    d_free(ub1[jj]);
                }
                int_free(idxbN);
                d_free(lbN);
                d_free(ubN);

                // LS cost and path constraints
                for (int_t i = 0; i <= N; i++) {
                    // Least-squares cost
                    free(ls_cost.fun[i]->in);
                    free(ls_cost.fun[i]->out);
                    free(ls_cost.fun[i]->args);
                    casadi_wrapper_destroy(ls_cost.fun[i]->work);
                    free(ls_cost.y_ref[i]);
                    free(ls_cost.fun[i]);
                    // Path constraints
                    free(path_constraints[i]->in);
                    free(path_constraints[i]->out);
                    free(path_constraints[i]->args);
                    casadi_wrapper_destroy(path_constraints[i]->work);
                    free(path_constraints[i]);
                }
                free(path_constraints);
                free(ls_cost.W);
                free(ls_cost.y_ref);

                // Integrators
                for (jj = 0; jj < N; jj++) {
                    free(sim_in[jj].x);
                    free(sim_in[jj].u);
                    free(sim_in[jj].S_forw);
                    free(sim_in[jj].S_adj);
                    free(sim_in[jj].grad_K);
                    free(sim_out[jj].xn);
                    free(sim_out[jj].S_forw);
                    free(sim_out[jj].grad);
                }

                // NLP memory and workspace
                ocp_nlp_sqp_destroy(nlp_mem, nlp_work);
                // NLP arguments
                free(nlp_args);
                // NLP output
                for (int_t i = 0; i <= N; i++) {
                    free(nlp_out.x[i]);
                    free(nlp_out.u[i]);
                    free(nlp_out.pi[i]);
                    free(nlp_out.lam[i]);
                }
                free(nlp_out.x);
                free(nlp_out.u);
                free(nlp_out.lam);
                free(nlp_out.pi);
            }
        }
    }
}
