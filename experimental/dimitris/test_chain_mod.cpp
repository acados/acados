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
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_qpdunes.h"
#include "acados/ocp_nlp/ocp_nlp_gn_sqp.h"
#include "acados/sim/sim_common.h"
#include "acados/sim/sim_erk_integrator.h"
#include "acados/sim/sim_lifted_irk_integrator.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"
#include "test/ocp_nlp/chain/Chain_model.h"
#include "test/test_utils/eigen.h"
#include "test/test_utils/read_matrix.h"
#include "test/test_utils/zeros.h"
#include "test/test_utils/read_ocp_qp_in.h"

#define NN 15
#define TT 3.0
#define Ns 2

using Eigen::MatrixXd;
using Eigen::VectorXd;

int main() {
    for (int INEXACT = 0; INEXACT < 5; INEXACT++) {
    int d_start = 0;
    if (INEXACT > 0) d_start = 2;

    for (int d = d_start; d < 4; d++) {  // RK4 in case d == 0
    for (int NMF = 1; NMF < 4; NMF++) {
        if (INEXACT == 0) {
            printf("\n----- NUMBER OF FREE MASSES = %d, d = %d (Exact Newton) -----\n", NMF, d);
        } else if (INEXACT == 1) {
            printf("\n----- NUMBER OF FREE MASSES = %d, d = %d (IN Scheme) -----\n", NMF, d);
        } else if (INEXACT == 2) {
            printf("\n----- NUMBER OF FREE MASSES = %d, d = %d (INIS Scheme) -----\n", NMF, d);
        } else if (INEXACT == 3) {
            printf("\n----- NUMBER OF FREE MASSES = %d, d = %d (FROZEN IN Scheme) -----\n", NMF, d);
        } else if (INEXACT == 4) {
            printf("\n----- NUMBER OF FREE MASSES = %d, d = %d (FROZEN INIS Scheme) -----\n",
                    NMF, d);
        }
        int_t NX = 6*NMF;
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

        d_zeros(&W, NX+NU, NX+NU);
        d_zeros(&WN, NX, NX);
        d_zeros(&uref, NU, 1);
        d_zeros(&x_end, NX, 1);
        d_zeros(&u_end, NU, 1);

        std::string NMFdat = std::to_string(NMF+1) + "_d" + std::to_string(d) + ".dat";
        std::string x0_str = "../../build/test/chain/x0_nm" + NMFdat;
        std::string xref_str = "../../build/test/chain/xN_nm" + NMFdat;
        std::string resX_str = "../../build/test/chain/resX_nm" + NMFdat;
        std::string resU_str = "../../build/test/chain/resU_nm" + NMFdat;
        real_t *x0_tmp = (real_t*)malloc(sizeof(real_t)*NX);
        real_t *xref_tmp = (real_t*)malloc(sizeof(real_t)*NX);
        real_t *resX_tmp = (real_t*)malloc(sizeof(real_t)*NX*(NN+1));
        real_t *resU_tmp = (real_t*)malloc(sizeof(real_t)*NU*(NN+1));
        read_double_vector_from_txt(x0_tmp, NX, x0_str.c_str());
        read_double_vector_from_txt(xref_tmp, NX, xref_str.c_str());
        read_double_matrix_from_txt(resX_tmp, NX, NN+1, resX_str.c_str());
        read_double_matrix_from_txt(resU_tmp, NU, NN, resU_str.c_str());
        VectorXd x0 = Eigen::Map<VectorXd>(x0_tmp, NX);
        VectorXd xref = Eigen::Map<VectorXd>(xref_tmp, NX);
        MatrixXd resX = Eigen::Map<MatrixXd>(resX_tmp, NX, N+1);
        MatrixXd resU = Eigen::Map<MatrixXd>(resU_tmp, NU, N);

        // NOTE: code above avoids readMatrix for valgrind
        // VectorXd x0 = readMatrix("../../build/test/chain/x0_nm" + NMFdat);
        // VectorXd xref = readMatrix("../../build/test/chain/xN_nm" + NMFdat);
        // MatrixXd resX = readMatrix("../../build/test/chain/resX_nm" + NMFdat);
        // MatrixXd resU = readMatrix("../../build/test/chain/resU_nm" + NMFdat);

        for (int_t i = 0; i < NX; i++) W[i*(NX+NU+1)] = 1e-2;
        for (int_t i = 0; i < NU; i++) W[(NX+i)*(NX+NU+1)] = 1.0;
        for (int_t i = 0; i < NX; i++) WN[i*(NX+1)] = 1e-2;

        ls_cost.W = (real_t **) malloc(sizeof(*ls_cost.W) * (N+1));
        for (int_t i = 0; i < NN; i++) ls_cost.W[i] = W;
        ls_cost.W[NN] = WN;
        ls_cost.y_ref = (real_t **) malloc(sizeof(*ls_cost.y_ref) * (N+1));
        for (int_t i = 0; i < NN; i++) {
            ls_cost.y_ref[i] = (real_t *) malloc(sizeof(*ls_cost.y_ref[i]) * (NX+NU));
            for (int_t j = 0; j < NX; j++) ls_cost.y_ref[i][j] = xref(j);
            for (int_t j = 0; j < NU; j++) ls_cost.y_ref[i][NX+j] = 0.0;
        }
        ls_cost.y_ref[N] = (real_t *) malloc(sizeof(*ls_cost.y_ref[N]) * (NX));
        for (int_t j = 0; j < NX; j++) ls_cost.y_ref[N][j] = xref(j);

        // Integrator structs
        real_t Ts = TT/NN;
        sim_in  sim_in[NN];
        sim_out sim_out[NN];
        sim_info info[NN];
        sim_solver integrators[NN];

        sim_RK_opts rk_opts[NN];
        sim_erk_workspace erk_work[NN];
        sim_lifted_irk_workspace irk_work[NN];
        sim_lifted_irk_memory irk_mem[NN];

        // TODO(rien): can I move this somewhere inside the integrator?
        struct d_strmat str_mat[NN];
        struct d_strmat str_sol[NN];

        for (jj = 0; jj < NN; jj++) {
            integrators[jj].in = &sim_in[jj];
            integrators[jj].out = &sim_out[jj];
            if (d > 0) {
                integrators[jj].fun = &sim_lifted_irk;
                integrators[jj].mem = &irk_mem[jj];
                integrators[jj].work = &irk_work[jj];
            } else {
                integrators[jj].fun = &sim_erk;
                integrators[jj].mem = 0;
                integrators[jj].work = &erk_work[jj];
            }

            sim_in[jj].nSteps = Ns;
            sim_in[jj].step = Ts/sim_in[jj].nSteps;
            sim_in[jj].nx = NX;
            sim_in[jj].nu = NU;

            sim_in[jj].opts = &rk_opts[jj];

            sim_in[jj].sens_forw = true;
            sim_in[jj].sens_adj = false;
            sim_in[jj].sens_hess = false;
            sim_in[jj].nsens_forw = NX+NU;

            switch (NMF) {
            case 1:
                sim_in[jj].VDE_forw = &VDE_fun_nm2;
                sim_in[jj].jac_fun = &jac_fun_nm2;
                break;
            case 2:
                sim_in[jj].VDE_forw = &VDE_fun_nm3;
                sim_in[jj].jac_fun = &jac_fun_nm3;
                break;
            case 3:
                sim_in[jj].VDE_forw = &VDE_fun_nm4;
                sim_in[jj].jac_fun = &jac_fun_nm4;
                break;
            default:
                break;
            }

            sim_in[jj].x = (real_t *) malloc(sizeof(*sim_in[jj].x) * (NX));
            sim_in[jj].u = (real_t *) malloc(sizeof(*sim_in[jj].u) * (NU));
            sim_in[jj].S_forw = (real_t *) malloc(sizeof(*sim_in[jj].S_forw) * (NX*(NX+NU)));
            for (int_t i = 0; i < NX*(NX+NU); i++) sim_in[jj].S_forw[i] = 0.0;
            for (int_t i = 0; i < NX; i++) sim_in[jj].S_forw[i*(NX+1)] = 1.0;

            sim_in[jj].S_adj = (real_t *) malloc(sizeof(*sim_in[jj].S_adj) * (NX+NU));
            for (int_t i = 0; i < NX+NU; i++) sim_in[jj].S_adj[i] = 0.0;

            sim_in[jj].grad_K = (real_t *) malloc(sizeof(*sim_in[jj].grad_K) * (d*NX));
            for (int_t i = 0; i < d*NX; i++) sim_in[jj].grad_K[i] = 0.0;

            sim_out[jj].xn = (real_t *) malloc(sizeof(*sim_out[jj].xn) * (NX));
            sim_out[jj].S_forw = (real_t *) malloc(sizeof(*sim_out[jj].S_forw) * (NX*(NX+NU)));
            sim_out[jj].info = &info[jj];
            sim_out[jj].grad = (real_t *) malloc(sizeof(*sim_out[jj].grad ) * (NX+NU));

            irk_work[jj].str_mat = &str_mat[jj];
            irk_work[jj].str_sol = &str_sol[jj];
            if (d > 0) {
                sim_irk_create_opts(d, "Gauss", &rk_opts[jj]);
                if (INEXACT == 0) {
                    sim_irk_create_Newton_scheme(d, "Gauss", &rk_opts[jj], exact);
                } else if (INEXACT == 1 || INEXACT == 3) {
                    sim_irk_create_Newton_scheme(d, "Gauss", &rk_opts[jj], simplified_in);
                } else if (INEXACT == 2 || INEXACT == 4) {
                    sim_irk_create_Newton_scheme(d, "Gauss", &rk_opts[jj], simplified_inis);
                }

                sim_lifted_irk_create_workspace(&sim_in[jj], &irk_work[jj]);
                sim_lifted_irk_create_memory(&sim_in[jj], &irk_mem[jj]);
            } else {
                sim_erk_create_arguments(4, &rk_opts[jj]);
                sim_erk_create_workspace(&sim_in[jj], &erk_work[jj]);
            }
        }

        int_t nx[NN+1] = {0};
        int_t nu[NN] = {0};
        int_t nb[NN+1] = {0};
        int_t nc[NN+1] = {0};
        int_t ng[NN+1] = {0};
        for (int_t i = 0; i < N; i++) {
            nx[i] = NX;
            nu[i] = NU;
        }
        nx[N] = NX;

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
            lb0[jj] = x0(jj);   //   xmin
            ub0[jj] = x0(jj);   //   xmax
            idxb0[jj] = jj;
        }
        for (; jj < NX+NU; jj++) {
            lb0[jj] = -UMAX;  //   umin
            ub0[jj] = UMAX;   //   umax
            idxb0[jj] = jj;
        }
        nb[0] = NX+NU;

        int *idxb1;
        int_zeros(&idxb1, NMF+NU, 1);
        double *lb1[N-1];
        double *ub1[N-1];
        for (int_t i = 0; i < N-1; i++) {
            d_zeros(&lb1[i], NMF+NU, 1);
            d_zeros(&ub1[i], NMF+NU, 1);
            for (jj = 0; jj < NMF; jj++) {
                lb1[i][jj] = wall_pos;      // wall position
                ub1[i][jj] = 1e8;  // NOTE: CAN'T BE 1e12 for OOQP ATM!!!!
                idxb1[jj] = 6*jj+1;
            }
            for (jj = 0; jj < NU; jj++) {
                lb1[i][NMF+jj] = -UMAX;  //   umin
                ub1[i][NMF+jj] = UMAX;   //   umax
                idxb1[NMF+jj] = NX+jj;
            }
            nb[i+1] = NMF+NU;
        }

        int *idxbN;
        int_zeros(&idxbN, NX, 1);
        real_t *lbN;
        d_zeros(&lbN, NX, 1);
        real_t *ubN;
        d_zeros(&ubN, NX, 1);
        for (jj = 0; jj < NX; jj++) {
            lbN[jj] = xref(jj);   //   xmin
            ubN[jj] = xref(jj);   //   xmax
            idxbN[jj] = jj;
        }
        nb[NN] = NX;

        real_t *hlb[N+1];
        real_t *hub[N+1];
        int *hidxb[N+1];

        hlb[0] = lb0;
        hub[0] = ub0;
        hidxb[0] = idxb0;
        for (int_t i = 1; i < N; i++) {
            hlb[i] = lb1[i-1];
            hub[i] = ub1[i-1];
            hidxb[i] = idxb1;
        }
        hlb[N] = lbN;
        hub[N] = ubN;
        hidxb[N] = idxbN;

        ocp_nlp_in nlp_in;
        nlp_in.N = NN;
        nlp_in.nx = nx;
        nlp_in.nu = nu;
        nlp_in.nb = nb;
        nlp_in.nc = nc;
        nlp_in.ng = ng;
        nlp_in.idxb = (const int_t **) hidxb;
        nlp_in.lb = (const real_t **) hlb;
        nlp_in.ub = (const real_t **) hub;
        nlp_in.sim = integrators;
        nlp_in.cost = &ls_cost;
        nlp_in.freezeSens = false;
        if (INEXACT > 2) nlp_in.freezeSens = true;

        ocp_nlp_out nlp_out;
        nlp_out.x = (real_t **) malloc(sizeof(*nlp_out.x) * (N+1));
        nlp_out.u = (real_t **) malloc(sizeof(*nlp_out.u) * N);
        nlp_out.lam = (real_t **) malloc(sizeof(*nlp_out.lam) * N);
        for (int_t i = 0; i < NN; i++) {
            nlp_out.x[i] = (real_t *) malloc(sizeof(*nlp_out.x[i]) * (NX));
            nlp_out.u[i] = (real_t *) malloc(sizeof(*nlp_out.u[i]) * (NU));
            nlp_out.lam[i] = (real_t *) malloc(sizeof(*nlp_out.lam[i]) * (NX));
        }
        nlp_out.x[N] = (real_t *) malloc(sizeof(*nlp_out.x[N]) * (NX));

        ocp_nlp_gn_sqp_opts nlp_args;
        ocp_nlp_args nlp_common_args;
        nlp_args.common = &nlp_common_args;
        nlp_args.common->maxIter = max_sqp_iters;

        sprintf(nlp_args.qp_solver_name, "ooqp");

        ocp_nlp_gn_sqp_memory nlp_mem;
        ocp_nlp_mem nlp_mem_common;
        nlp_mem.common = &nlp_mem_common;
        ocp_nlp_gn_sqp_create_memory(&nlp_in, &nlp_args, &nlp_mem);

        int_t work_space_size = ocp_nlp_gn_sqp_calculate_workspace_size(&nlp_in, &nlp_args);
        void *nlp_work = (void*)malloc(work_space_size);

        if (0) {
            for (int_t i = 0; i < NN; i++) {
                for (int_t j = 0; j < NX; j++) nlp_mem.common->x[i][j] = resX(j, i);// xref(j);  // resX(j, i);
                for (int_t j = 0; j < NU; j++) nlp_mem.common->u[i][j] = resU(j, i); //0.0;  // resU(j, i);
            }
        for (int_t j = 0; j < NX; j++) nlp_mem.common->x[NN][j] = resX(j, N);//xref(j);  // resX(j, N);
        } else {
            for (int_t i = 0; i < NN; i++) {
                for (int_t j = 0; j < NX; j++) nlp_mem.common->x[i][j] = xref(j);
                for (int_t j = 0; j < NU; j++) nlp_mem.common->u[i][j] = 0.0;
            }
        for (int_t j = 0; j < NX; j++) nlp_mem.common->x[NN][j] = xref(j);
        }

        int_t status;
        status = ocp_nlp_gn_sqp(&nlp_in, &nlp_out, &nlp_args, &nlp_mem, nlp_work);
        ocp_nlp_gn_sqp_free_memory(&nlp_mem);

        real_t out_x[NX*(N+1)], err_x[NX*(N+1)];
        real_t out_u[NU*N], err_u[NU*N];
        for (int_t i = 0; i < NN; i++) {
            for (int_t j = 0; j < NX; j++) out_x[i*NX+j] = nlp_out.x[i][j];
            for (int_t j = 0; j < NU; j++) out_u[i*NU+j] = nlp_out.u[i][j];
        }
        for (int_t j = 0; j < NX; j++) out_x[N*NX+j] = nlp_out.x[N][j];

        for (int_t i = 0; i < NN; i++) {
            for (int_t j = 0; j < NX; j++) err_x[i*NX+j] = fabs(out_x[i*NX+j] - resX(j, i));
            for (int_t j = 0; j < NU; j++) err_u[i*NU+j] = fabs(out_u[i*NU+j] - resU(j, i));
        }
        for (int_t j = 0; j < NX; j++) err_x[N*NX+j] = fabs(out_x[N*NX+j] - resX(j, N));

//        print_matrix_name((char*)"stdout", (char*)"out_x", out_x, NX, N+1);
//        print_matrix_name((char*)"stdout", (char*)"out_u", out_u, NU, N);

        print_matrix_name((char*)"stdout", (char*)"err_x", err_x, NX, N+1);
        print_matrix_name((char*)"stdout", (char*)"err_u", err_u, NU, N);

        MatrixXd SQP_x = Eigen::Map<MatrixXd>(&out_x[0], NX, N+1);
        MatrixXd SQP_u = Eigen::Map<MatrixXd>(&out_u[0], NU, N);
        printf("correct u:\n");
        std::cout << resU << std::endl;
        printf("sqp u:\n");
        std::cout << SQP_u << std::endl;
        // printf("correct x:\n");
        // std::cout << resX << std::endl;
        // printf("sqp x:\n");
        // std::cout << SQP_x << std::endl;

        // REQUIRE(SQP_x.isApprox(resX, COMPARISON_TOLERANCE_IPOPT));
        // REQUIRE(SQP_u.isApprox(resU, COMPARISON_TOLERANCE_IPOPT));

        d_free(W);
        d_free(WN);
        d_free(uref);
        d_free(x_end);
        d_free(u_end);

        int_free(idxb0);
        d_free(lb0);
        d_free(ub0);
        int_free(idxb1);
        for (jj = 0; jj < N-1; jj++) {
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
        free(ls_cost.y_ref[N]);
        free(ls_cost.y_ref);
        free(ls_cost.W);

        for (jj = 0; jj < N; jj++) {
            free(nlp_out.x[jj]);
            free(nlp_out.u[jj]);
            free(nlp_out.lam[jj]);
        }
        free(nlp_out.x[N]);
        free(nlp_out.x);
        free(nlp_out.u);
        free(nlp_out.lam);

        free(nlp_work);

        if (1) {
            free(x0_tmp);
            free(xref_tmp);
            free(resX_tmp);
            free(resU_tmp);
        }
    }
    }
    }
    return 0;
}
