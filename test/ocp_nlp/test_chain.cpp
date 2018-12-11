/*
 *    This file is part of acados.
 *
 *    acados is free software; you can redistribute it and/or modify it under
 *    the terms of the GNU Lesser General Public License as published by the
 *    Free Software Foundation; either version 3 of the License, or (at your
 *    option) any later version.
 *
 *    acados is distributed in the hope that it will be useful, but WITHOUT ANY
 *    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 *    FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
 *    more details.
 *
 *    You should have received a copy of the GNU Lesser General Public License
 *    along with acados; if not, write to the Free Software Foundation, Inc., 51
 *    Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
 *
 */

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>

#include "test/test_utils/eigen.h"
#include "catch/include/catch.hpp"

#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"

#include "acados_c/external_function_interface.h"
#include "acados_c/ocp_qp_interface.h"
#include "acados_c/ocp_nlp_interface.h"

// TODO(dimitris): use only the strictly necessary includes here

#include "acados/utils/mem.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"

#include "acados/ocp_nlp/ocp_nlp_sqp.h"
#include "acados/ocp_nlp/ocp_nlp_cost_common.h"
#include "acados/ocp_nlp/ocp_nlp_cost_ls.h"
#include "acados/ocp_nlp/ocp_nlp_cost_nls.h"
#include "acados/ocp_nlp/ocp_nlp_cost_external.h"
#include "acados/ocp_nlp/ocp_nlp_dynamics_cont.h"
#include "acados/ocp_nlp/ocp_nlp_dynamics_disc.h"
#include "acados/ocp_nlp/ocp_nlp_constraints_bgh.h"

#include "examples/c/chain_model/chain_model.h"
#include "examples/c/implicit_chain_model/chain_model_impl.h"

// x0
#include "examples/c/chain_model/x0_nm2.c"
#include "examples/c/chain_model/x0_nm3.c"
#include "examples/c/chain_model/x0_nm4.c"
#include "examples/c/chain_model/x0_nm5.c"
#include "examples/c/chain_model/x0_nm6.c"

// xN
#include "examples/c/chain_model/xN_nm2.c"
#include "examples/c/chain_model/xN_nm3.c"
#include "examples/c/chain_model/xN_nm4.c"
#include "examples/c/chain_model/xN_nm5.c"
#include "examples/c/chain_model/xN_nm6.c"

#define TF 3.75
#define MAX_SQP_ITERS 20
#define NREP 10
#define TOL 1e-6

typedef enum {
    BOX = 0,
    GENERAL,
    GENERAL_NONLINEAR
} constraints_t;

ocp_qp_solver_t qp_solver_enum(std::string const& inString)
{
    if (inString == "SPARSE_HPIPM") return PARTIAL_CONDENSING_HPIPM;
    if (inString == "SPARSE_HPMPC") return PARTIAL_CONDENSING_HPMPC;
    if (inString == "SPARSE_QPDUNES") return PARTIAL_CONDENSING_QPDUNES;

    if (inString == "DENSE_HPIPM") return FULL_CONDENSING_HPIPM;
    if (inString == "DENSE_QPOASES") return FULL_CONDENSING_QPOASES;
#ifdef ACADOS_WITH_QORE
    if (inString == "DENSE_QORE") return FULL_CONDENSING_QORE;
#endif
#ifdef ACADOS_WITH_OOQP
    if (inString == "DENSE_OOQP") return FULL_CONDENSING_OOQP;
    if (inString == "SPARSE_OOQP") return PARTIAL_CONDENSING_OOQP;
#endif

    return (ocp_qp_solver_t) -1;
}

constraints_t constraints_enum(std::string const& inString)
{
    if (inString == "BOX") return BOX;
    if (inString == "GENERAL") return GENERAL;
    if (inString == "NONLINEAR+GENERAL") return GENERAL_NONLINEAR;

    return (constraints_t) -1;
}

ocp_nlp_dynamics_t nlp_dynamics_enum(std::string const& inString)
{
    if (inString == "CONTINUOUS") return CONTINUOUS_MODEL;
    if (inString == "DISCRETE") return DISCRETE_MODEL;

    return (ocp_nlp_dynamics_t) -1;
}

sim_solver_t integrator_enum(std::string const& inString)
{
    if (inString == "ERK") return ERK;
    if (inString == "IRK") return IRK;
    if (inString == "LIFTED_IRK") return LIFTED_IRK;

    return (sim_solver_t) -1;
}

ocp_nlp_cost_t cost_enum(std::string const& inString)
{
    if (inString == "LINEAR_LS") return LINEAR_LS;
    if (inString == "NONLINEAR_LS") return NONLINEAR_LS;
    if (inString == "EXTERNAL") return EXTERNALLY_PROVIDED;

    return (ocp_nlp_cost_t) -1;
}


static void select_dynamics_casadi(int N, int num_free_masses,
    external_function_casadi *forw_vde,
    external_function_casadi *impl_ode_fun,
    external_function_casadi *impl_ode_fun_jac_x_xdot,
    external_function_casadi *impl_ode_fun_jac_x_xdot_u,
    external_function_casadi *impl_ode_jac_x_xdot_u,
    external_function_casadi *erk4_casadi)
{
    switch (num_free_masses)
    {
        case 1:
            for (int ii = 0; ii < N; ii++)
            {
                forw_vde[ii].casadi_fun = &vde_chain_nm2;
                forw_vde[ii].casadi_work = &vde_chain_nm2_work;
                forw_vde[ii].casadi_sparsity_in = &vde_chain_nm2_sparsity_in;
                forw_vde[ii].casadi_sparsity_out = &vde_chain_nm2_sparsity_out;
                forw_vde[ii].casadi_n_in = &vde_chain_nm2_n_in;
                forw_vde[ii].casadi_n_out = &vde_chain_nm2_n_out;

                impl_ode_fun[ii].casadi_fun = &casadi_impl_ode_fun_chain_nm2;
                impl_ode_fun[ii].casadi_work = &casadi_impl_ode_fun_chain_nm2_work;
                impl_ode_fun[ii].casadi_sparsity_in = &casadi_impl_ode_fun_chain_nm2_sparsity_in;
                impl_ode_fun[ii].casadi_sparsity_out = &casadi_impl_ode_fun_chain_nm2_sparsity_out;
                impl_ode_fun[ii].casadi_n_in = &casadi_impl_ode_fun_chain_nm2_n_in;
                impl_ode_fun[ii].casadi_n_out = &casadi_impl_ode_fun_chain_nm2_n_out;

                impl_ode_fun_jac_x_xdot[ii].casadi_fun =
                &casadi_impl_ode_fun_jac_x_xdot_chain_nm2;
                impl_ode_fun_jac_x_xdot[ii].casadi_work =
                &casadi_impl_ode_fun_jac_x_xdot_chain_nm2_work;
                impl_ode_fun_jac_x_xdot[ii].casadi_sparsity_in =
                &casadi_impl_ode_fun_jac_x_xdot_chain_nm2_sparsity_in;
                impl_ode_fun_jac_x_xdot[ii].casadi_sparsity_out =
                &casadi_impl_ode_fun_jac_x_xdot_chain_nm2_sparsity_out;
                impl_ode_fun_jac_x_xdot[ii].casadi_n_in =
                &casadi_impl_ode_fun_jac_x_xdot_chain_nm2_n_in;
                impl_ode_fun_jac_x_xdot[ii].casadi_n_out =
                &casadi_impl_ode_fun_jac_x_xdot_chain_nm2_n_out;

                impl_ode_fun_jac_x_xdot_u[ii].casadi_fun =
                &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm2;
                impl_ode_fun_jac_x_xdot_u[ii].casadi_work =
                &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm2_work;
                impl_ode_fun_jac_x_xdot_u[ii].casadi_sparsity_in =
                &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm2_sparsity_in;
                impl_ode_fun_jac_x_xdot_u[ii].casadi_sparsity_out =
                &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm2_sparsity_out;
                impl_ode_fun_jac_x_xdot_u[ii].casadi_n_in =
                &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm2_n_in;
                impl_ode_fun_jac_x_xdot_u[ii].casadi_n_out =
                &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm2_n_out;

                impl_ode_jac_x_xdot_u[ii].casadi_fun =
                &casadi_impl_ode_jac_x_xdot_u_chain_nm2;
                impl_ode_jac_x_xdot_u[ii].casadi_work =
                &casadi_impl_ode_jac_x_xdot_u_chain_nm2_work;
                impl_ode_jac_x_xdot_u[ii].casadi_sparsity_in =
                &casadi_impl_ode_jac_x_xdot_u_chain_nm2_sparsity_in;
                impl_ode_jac_x_xdot_u[ii].casadi_sparsity_out =
                &casadi_impl_ode_jac_x_xdot_u_chain_nm2_sparsity_out;
                impl_ode_jac_x_xdot_u[ii].casadi_n_in =
                &casadi_impl_ode_jac_x_xdot_u_chain_nm2_n_in;
                impl_ode_jac_x_xdot_u[ii].casadi_n_out =
                &casadi_impl_ode_jac_x_xdot_u_chain_nm2_n_out;

                erk4_casadi[ii].casadi_fun = &casadi_erk4_chain_nm2;
                erk4_casadi[ii].casadi_work = &casadi_erk4_chain_nm2_work;
                erk4_casadi[ii].casadi_sparsity_in = &casadi_erk4_chain_nm2_sparsity_in;
                erk4_casadi[ii].casadi_sparsity_out = &casadi_erk4_chain_nm2_sparsity_out;
                erk4_casadi[ii].casadi_n_in = &casadi_erk4_chain_nm2_n_in;
                erk4_casadi[ii].casadi_n_out = &casadi_erk4_chain_nm2_n_out;
            }
            break;
        case 2:
            for (int ii = 0; ii < N; ii++)
            {
                forw_vde[ii].casadi_fun = &vde_chain_nm3;
                forw_vde[ii].casadi_work = &vde_chain_nm3_work;
                forw_vde[ii].casadi_sparsity_in = &vde_chain_nm3_sparsity_in;
                forw_vde[ii].casadi_sparsity_out = &vde_chain_nm3_sparsity_out;
                forw_vde[ii].casadi_n_in = &vde_chain_nm3_n_in;
                forw_vde[ii].casadi_n_out = &vde_chain_nm3_n_out;

                impl_ode_fun[ii].casadi_fun = &casadi_impl_ode_fun_chain_nm3;
                impl_ode_fun[ii].casadi_work = &casadi_impl_ode_fun_chain_nm3_work;
                impl_ode_fun[ii].casadi_sparsity_in = &casadi_impl_ode_fun_chain_nm3_sparsity_in;
                impl_ode_fun[ii].casadi_sparsity_out = &casadi_impl_ode_fun_chain_nm3_sparsity_out;
                impl_ode_fun[ii].casadi_n_in = &casadi_impl_ode_fun_chain_nm3_n_in;
                impl_ode_fun[ii].casadi_n_out = &casadi_impl_ode_fun_chain_nm3_n_out;

                impl_ode_fun_jac_x_xdot[ii].casadi_fun =
                &casadi_impl_ode_fun_jac_x_xdot_chain_nm3;
                impl_ode_fun_jac_x_xdot[ii].casadi_work =
                &casadi_impl_ode_fun_jac_x_xdot_chain_nm3_work;
                impl_ode_fun_jac_x_xdot[ii].casadi_sparsity_in =
                &casadi_impl_ode_fun_jac_x_xdot_chain_nm3_sparsity_in;
                impl_ode_fun_jac_x_xdot[ii].casadi_sparsity_out =
                &casadi_impl_ode_fun_jac_x_xdot_chain_nm3_sparsity_out;
                impl_ode_fun_jac_x_xdot[ii].casadi_n_in =
                &casadi_impl_ode_fun_jac_x_xdot_chain_nm3_n_in;
                impl_ode_fun_jac_x_xdot[ii].casadi_n_out =
                &casadi_impl_ode_fun_jac_x_xdot_chain_nm3_n_out;

                impl_ode_fun_jac_x_xdot_u[ii].casadi_fun =
                &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm3;
                impl_ode_fun_jac_x_xdot_u[ii].casadi_work =
                &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm3_work;
                impl_ode_fun_jac_x_xdot_u[ii].casadi_sparsity_in =
                &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm3_sparsity_in;
                impl_ode_fun_jac_x_xdot_u[ii].casadi_sparsity_out =
                &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm3_sparsity_out;
                impl_ode_fun_jac_x_xdot_u[ii].casadi_n_in =
                &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm3_n_in;
                impl_ode_fun_jac_x_xdot_u[ii].casadi_n_out =
                &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm3_n_out;

                impl_ode_jac_x_xdot_u[ii].casadi_fun =
                &casadi_impl_ode_jac_x_xdot_u_chain_nm3;
                impl_ode_jac_x_xdot_u[ii].casadi_work =
                &casadi_impl_ode_jac_x_xdot_u_chain_nm3_work;
                impl_ode_jac_x_xdot_u[ii].casadi_sparsity_in =
                &casadi_impl_ode_jac_x_xdot_u_chain_nm3_sparsity_in;
                impl_ode_jac_x_xdot_u[ii].casadi_sparsity_out =
                &casadi_impl_ode_jac_x_xdot_u_chain_nm3_sparsity_out;
                impl_ode_jac_x_xdot_u[ii].casadi_n_in =
                &casadi_impl_ode_jac_x_xdot_u_chain_nm3_n_in;
                impl_ode_jac_x_xdot_u[ii].casadi_n_out =
                &casadi_impl_ode_jac_x_xdot_u_chain_nm3_n_out;

                erk4_casadi[ii].casadi_fun = &casadi_erk4_chain_nm3;
                erk4_casadi[ii].casadi_work = &casadi_erk4_chain_nm3_work;
                erk4_casadi[ii].casadi_sparsity_in = &casadi_erk4_chain_nm3_sparsity_in;
                erk4_casadi[ii].casadi_sparsity_out = &casadi_erk4_chain_nm3_sparsity_out;
                erk4_casadi[ii].casadi_n_in = &casadi_erk4_chain_nm3_n_in;
                erk4_casadi[ii].casadi_n_out = &casadi_erk4_chain_nm3_n_out;
            }
            break;
        case 3:
            for (int ii = 0; ii < N; ii++)
            {
                forw_vde[ii].casadi_fun = &vde_chain_nm4;
                forw_vde[ii].casadi_work = &vde_chain_nm4_work;
                forw_vde[ii].casadi_sparsity_in = &vde_chain_nm4_sparsity_in;
                forw_vde[ii].casadi_sparsity_out = &vde_chain_nm4_sparsity_out;
                forw_vde[ii].casadi_n_in = &vde_chain_nm4_n_in;
                forw_vde[ii].casadi_n_out = &vde_chain_nm4_n_out;

                impl_ode_fun[ii].casadi_fun = &casadi_impl_ode_fun_chain_nm4;
                impl_ode_fun[ii].casadi_work = &casadi_impl_ode_fun_chain_nm4_work;
                impl_ode_fun[ii].casadi_sparsity_in = &casadi_impl_ode_fun_chain_nm4_sparsity_in;
                impl_ode_fun[ii].casadi_sparsity_out = &casadi_impl_ode_fun_chain_nm4_sparsity_out;
                impl_ode_fun[ii].casadi_n_in = &casadi_impl_ode_fun_chain_nm4_n_in;
                impl_ode_fun[ii].casadi_n_out = &casadi_impl_ode_fun_chain_nm4_n_out;

                impl_ode_fun_jac_x_xdot[ii].casadi_fun = &casadi_impl_ode_fun_jac_x_xdot_chain_nm4;
                impl_ode_fun_jac_x_xdot[ii].casadi_work =
                &casadi_impl_ode_fun_jac_x_xdot_chain_nm4_work;
                impl_ode_fun_jac_x_xdot[ii].casadi_sparsity_in =
                &casadi_impl_ode_fun_jac_x_xdot_chain_nm4_sparsity_in;
                impl_ode_fun_jac_x_xdot[ii].casadi_sparsity_out =
                &casadi_impl_ode_fun_jac_x_xdot_chain_nm4_sparsity_out;
                impl_ode_fun_jac_x_xdot[ii].casadi_n_in =
                &casadi_impl_ode_fun_jac_x_xdot_chain_nm4_n_in;
                impl_ode_fun_jac_x_xdot[ii].casadi_n_out =
                &casadi_impl_ode_fun_jac_x_xdot_chain_nm4_n_out;

                impl_ode_fun_jac_x_xdot_u[ii].casadi_fun =
                &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm4;
                impl_ode_fun_jac_x_xdot_u[ii].casadi_work =
                &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm4_work;
                impl_ode_fun_jac_x_xdot_u[ii].casadi_sparsity_in =
                &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm4_sparsity_in;
                impl_ode_fun_jac_x_xdot_u[ii].casadi_sparsity_out =
                &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm4_sparsity_out;
                impl_ode_fun_jac_x_xdot_u[ii].casadi_n_in =
                &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm4_n_in;
                impl_ode_fun_jac_x_xdot_u[ii].casadi_n_out =
                &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm4_n_out;

                impl_ode_jac_x_xdot_u[ii].casadi_fun =
                &casadi_impl_ode_jac_x_xdot_u_chain_nm4;
                impl_ode_jac_x_xdot_u[ii].casadi_work =
                &casadi_impl_ode_jac_x_xdot_u_chain_nm4_work;
                impl_ode_jac_x_xdot_u[ii].casadi_sparsity_in =
                &casadi_impl_ode_jac_x_xdot_u_chain_nm4_sparsity_in;
                impl_ode_jac_x_xdot_u[ii].casadi_sparsity_out =
                &casadi_impl_ode_jac_x_xdot_u_chain_nm4_sparsity_out;
                impl_ode_jac_x_xdot_u[ii].casadi_n_in =
                &casadi_impl_ode_jac_x_xdot_u_chain_nm4_n_in;
                impl_ode_jac_x_xdot_u[ii].casadi_n_out =
                &casadi_impl_ode_jac_x_xdot_u_chain_nm4_n_out;

                erk4_casadi[ii].casadi_fun = &casadi_erk4_chain_nm4;
                erk4_casadi[ii].casadi_work = &casadi_erk4_chain_nm4_work;
                erk4_casadi[ii].casadi_sparsity_in = &casadi_erk4_chain_nm4_sparsity_in;
                erk4_casadi[ii].casadi_sparsity_out = &casadi_erk4_chain_nm4_sparsity_out;
                erk4_casadi[ii].casadi_n_in = &casadi_erk4_chain_nm4_n_in;
                erk4_casadi[ii].casadi_n_out = &casadi_erk4_chain_nm4_n_out;
            }
            break;
        case 4:
            for (int ii = 0; ii < N; ii++)
            {
                forw_vde[ii].casadi_fun = &vde_chain_nm5;
                forw_vde[ii].casadi_work = &vde_chain_nm5_work;
                forw_vde[ii].casadi_sparsity_in = &vde_chain_nm5_sparsity_in;
                forw_vde[ii].casadi_sparsity_out = &vde_chain_nm5_sparsity_out;
                forw_vde[ii].casadi_n_in = &vde_chain_nm5_n_in;
                forw_vde[ii].casadi_n_out = &vde_chain_nm5_n_out;

                impl_ode_fun[ii].casadi_fun = &casadi_impl_ode_fun_chain_nm5;
                impl_ode_fun[ii].casadi_work = &casadi_impl_ode_fun_chain_nm5_work;
                impl_ode_fun[ii].casadi_sparsity_in = &casadi_impl_ode_fun_chain_nm5_sparsity_in;
                impl_ode_fun[ii].casadi_sparsity_out = &casadi_impl_ode_fun_chain_nm5_sparsity_out;
                impl_ode_fun[ii].casadi_n_in = &casadi_impl_ode_fun_chain_nm5_n_in;
                impl_ode_fun[ii].casadi_n_out = &casadi_impl_ode_fun_chain_nm5_n_out;

                impl_ode_fun_jac_x_xdot[ii].casadi_fun = &casadi_impl_ode_fun_jac_x_xdot_chain_nm5;
                impl_ode_fun_jac_x_xdot[ii].casadi_work =
                &casadi_impl_ode_fun_jac_x_xdot_chain_nm5_work;
                impl_ode_fun_jac_x_xdot[ii].casadi_sparsity_in =
                &casadi_impl_ode_fun_jac_x_xdot_chain_nm5_sparsity_in;
                impl_ode_fun_jac_x_xdot[ii].casadi_sparsity_out =
                &casadi_impl_ode_fun_jac_x_xdot_chain_nm5_sparsity_out;
                impl_ode_fun_jac_x_xdot[ii].casadi_n_in =
                &casadi_impl_ode_fun_jac_x_xdot_chain_nm5_n_in;
                impl_ode_fun_jac_x_xdot[ii].casadi_n_out =
                &casadi_impl_ode_fun_jac_x_xdot_chain_nm5_n_out;

                impl_ode_fun_jac_x_xdot_u[ii].casadi_fun =
                &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm5;
                impl_ode_fun_jac_x_xdot_u[ii].casadi_work =
                &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm5_work;
                impl_ode_fun_jac_x_xdot_u[ii].casadi_sparsity_in =
                &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm5_sparsity_in;
                impl_ode_fun_jac_x_xdot_u[ii].casadi_sparsity_out =
                &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm5_sparsity_out;
                impl_ode_fun_jac_x_xdot_u[ii].casadi_n_in =
                &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm5_n_in;
                impl_ode_fun_jac_x_xdot_u[ii].casadi_n_out =
                &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm5_n_out;

                impl_ode_jac_x_xdot_u[ii].casadi_fun = &casadi_impl_ode_jac_x_xdot_u_chain_nm5;
                impl_ode_jac_x_xdot_u[ii].casadi_work =
                &casadi_impl_ode_jac_x_xdot_u_chain_nm5_work;
                impl_ode_jac_x_xdot_u[ii].casadi_sparsity_in =
                &casadi_impl_ode_jac_x_xdot_u_chain_nm5_sparsity_in;
                impl_ode_jac_x_xdot_u[ii].casadi_sparsity_out =
                &casadi_impl_ode_jac_x_xdot_u_chain_nm5_sparsity_out;
                impl_ode_jac_x_xdot_u[ii].casadi_n_in =
                &casadi_impl_ode_jac_x_xdot_u_chain_nm5_n_in;
                impl_ode_jac_x_xdot_u[ii].casadi_n_out =
                &casadi_impl_ode_jac_x_xdot_u_chain_nm5_n_out;

                // erk4_casadi[ii].casadi_fun = &casadi_erk4_chain_nm5;
                // erk4_casadi[ii].casadi_work = &casadi_erk4_chain_nm5_work;
                // erk4_casadi[ii].casadi_sparsity_in = &casadi_erk4_chain_nm5_sparsity_in;
                // erk4_casadi[ii].casadi_sparsity_out = &casadi_erk4_chain_nm5_sparsity_out;
                // erk4_casadi[ii].casadi_n_in = &casadi_erk4_chain_nm5_n_in;
                // erk4_casadi[ii].casadi_n_out = &casadi_erk4_chain_nm5_n_out;
            }
            break;
        case 5:
            for (int ii = 0; ii < N; ii++)
            {
                forw_vde[ii].casadi_fun = &vde_chain_nm6;
                forw_vde[ii].casadi_work = &vde_chain_nm6_work;
                forw_vde[ii].casadi_sparsity_in = &vde_chain_nm6_sparsity_in;
                forw_vde[ii].casadi_sparsity_out = &vde_chain_nm6_sparsity_out;
                forw_vde[ii].casadi_n_in = &vde_chain_nm6_n_in;
                forw_vde[ii].casadi_n_out = &vde_chain_nm6_n_out;

                impl_ode_fun[ii].casadi_fun = &casadi_impl_ode_fun_chain_nm6;
                impl_ode_fun[ii].casadi_work = &casadi_impl_ode_fun_chain_nm6_work;
                impl_ode_fun[ii].casadi_sparsity_in = &casadi_impl_ode_fun_chain_nm6_sparsity_in;
                impl_ode_fun[ii].casadi_sparsity_out = &casadi_impl_ode_fun_chain_nm6_sparsity_out;
                impl_ode_fun[ii].casadi_n_in = &casadi_impl_ode_fun_chain_nm6_n_in;
                impl_ode_fun[ii].casadi_n_out = &casadi_impl_ode_fun_chain_nm6_n_out;

                impl_ode_fun_jac_x_xdot[ii].casadi_fun = &casadi_impl_ode_fun_jac_x_xdot_chain_nm6;
                impl_ode_fun_jac_x_xdot[ii].casadi_work =
                &casadi_impl_ode_fun_jac_x_xdot_chain_nm6_work;
                impl_ode_fun_jac_x_xdot[ii].casadi_sparsity_in =
                &casadi_impl_ode_fun_jac_x_xdot_chain_nm6_sparsity_in;
                impl_ode_fun_jac_x_xdot[ii].casadi_sparsity_out =
                &casadi_impl_ode_fun_jac_x_xdot_chain_nm6_sparsity_out;
                impl_ode_fun_jac_x_xdot[ii].casadi_n_in =
                &casadi_impl_ode_fun_jac_x_xdot_chain_nm6_n_in;
                impl_ode_fun_jac_x_xdot[ii].casadi_n_out =
                &casadi_impl_ode_fun_jac_x_xdot_chain_nm6_n_out;

                impl_ode_fun_jac_x_xdot_u[ii].casadi_fun =
                &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm6;
                impl_ode_fun_jac_x_xdot_u[ii].casadi_work =
                &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm6_work;
                impl_ode_fun_jac_x_xdot_u[ii].casadi_sparsity_in =
                &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm6_sparsity_in;
                impl_ode_fun_jac_x_xdot_u[ii].casadi_sparsity_out =
                &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm6_sparsity_out;
                impl_ode_fun_jac_x_xdot_u[ii].casadi_n_in =
                &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm6_n_in;
                impl_ode_fun_jac_x_xdot_u[ii].casadi_n_out =
                &casadi_impl_ode_fun_jac_x_xdot_u_chain_nm6_n_out;

                impl_ode_jac_x_xdot_u[ii].casadi_fun = &casadi_impl_ode_jac_x_xdot_u_chain_nm6;
                impl_ode_jac_x_xdot_u[ii].casadi_work =
                &casadi_impl_ode_jac_x_xdot_u_chain_nm6_work;
                impl_ode_jac_x_xdot_u[ii].casadi_sparsity_in =
                &casadi_impl_ode_jac_x_xdot_u_chain_nm6_sparsity_in;
                impl_ode_jac_x_xdot_u[ii].casadi_sparsity_out =
                &casadi_impl_ode_jac_x_xdot_u_chain_nm6_sparsity_out;
                impl_ode_jac_x_xdot_u[ii].casadi_n_in =
                &casadi_impl_ode_jac_x_xdot_u_chain_nm6_n_in;
                impl_ode_jac_x_xdot_u[ii].casadi_n_out =
                &casadi_impl_ode_jac_x_xdot_u_chain_nm6_n_out;

                // erk4_casadi[ii].casadi_fun = &casadi_erk4_chain_nm6;
                // erk4_casadi[ii].casadi_work = &casadi_erk4_chain_nm6_work;
                // erk4_casadi[ii].casadi_sparsity_in = &casadi_erk4_chain_nm6_sparsity_in;
                // erk4_casadi[ii].casadi_sparsity_out = &casadi_erk4_chain_nm6_sparsity_out;
                // erk4_casadi[ii].casadi_n_in = &casadi_erk4_chain_nm6_n_in;
                // erk4_casadi[ii].casadi_n_out = &casadi_erk4_chain_nm6_n_out;
            }
            break;
        default:
            printf("Problem size not available\n");
            exit(1);
            break;
    }
    return;
}



static void select_ls_stage_cost_jac_casadi(int indx,
    int N,
    int num_free_masses,
    external_function_casadi *ls_cost_jac)
{
    switch (num_free_masses)
    {
        case 1:
            if (indx < N)
            {
                ls_cost_jac->casadi_fun = &ls_cost_nm2;
                ls_cost_jac->casadi_work = &ls_cost_nm2_work;
                ls_cost_jac->casadi_sparsity_in = &ls_cost_nm2_sparsity_in;
                ls_cost_jac->casadi_sparsity_out = &ls_cost_nm2_sparsity_out;
                ls_cost_jac->casadi_n_in = &ls_cost_nm2_n_in;
                ls_cost_jac->casadi_n_out = &ls_cost_nm2_n_out;
            }
            else
            {
                ls_cost_jac->casadi_fun = &ls_costN_nm2;
                ls_cost_jac->casadi_work = &ls_costN_nm2_work;
                ls_cost_jac->casadi_sparsity_in = &ls_costN_nm2_sparsity_in;
                ls_cost_jac->casadi_sparsity_out = &ls_costN_nm2_sparsity_out;
                ls_cost_jac->casadi_n_in = &ls_costN_nm2_n_in;
                ls_cost_jac->casadi_n_out = &ls_costN_nm2_n_out;
            }
            break;
        case 2:
            if (indx < N)
            {
                ls_cost_jac->casadi_fun = &ls_cost_nm3;
                ls_cost_jac->casadi_work = &ls_cost_nm3_work;
                ls_cost_jac->casadi_sparsity_in = &ls_cost_nm3_sparsity_in;
                ls_cost_jac->casadi_sparsity_out = &ls_cost_nm3_sparsity_out;
                ls_cost_jac->casadi_n_in = &ls_cost_nm3_n_in;
                ls_cost_jac->casadi_n_out = &ls_cost_nm3_n_out;
            }
            else
            {
                ls_cost_jac->casadi_fun = &ls_costN_nm3;
                ls_cost_jac->casadi_work = &ls_costN_nm3_work;
                ls_cost_jac->casadi_sparsity_in = &ls_costN_nm3_sparsity_in;
                ls_cost_jac->casadi_sparsity_out = &ls_costN_nm3_sparsity_out;
                ls_cost_jac->casadi_n_in = &ls_costN_nm3_n_in;
                ls_cost_jac->casadi_n_out = &ls_costN_nm3_n_out;
            }
            break;
        case 3:
            if (indx < N)
            {
                ls_cost_jac->casadi_fun = &ls_cost_nm4;
                ls_cost_jac->casadi_work = &ls_cost_nm4_work;
                ls_cost_jac->casadi_sparsity_in = &ls_cost_nm4_sparsity_in;
                ls_cost_jac->casadi_sparsity_out = &ls_cost_nm4_sparsity_out;
                ls_cost_jac->casadi_n_in = &ls_cost_nm4_n_in;
                ls_cost_jac->casadi_n_out = &ls_cost_nm4_n_out;
            }
            else
            {
                ls_cost_jac->casadi_fun = &ls_costN_nm4;
                ls_cost_jac->casadi_work = &ls_costN_nm4_work;
                ls_cost_jac->casadi_sparsity_in = &ls_costN_nm4_sparsity_in;
                ls_cost_jac->casadi_sparsity_out = &ls_costN_nm4_sparsity_out;
                ls_cost_jac->casadi_n_in = &ls_costN_nm4_n_in;
                ls_cost_jac->casadi_n_out = &ls_costN_nm4_n_out;
            }
            break;
        case 4:
            if (indx < N)
            {
                ls_cost_jac->casadi_fun = &ls_cost_nm5;
                ls_cost_jac->casadi_work = &ls_cost_nm5_work;
                ls_cost_jac->casadi_sparsity_in = &ls_cost_nm5_sparsity_in;
                ls_cost_jac->casadi_sparsity_out = &ls_cost_nm5_sparsity_out;
                ls_cost_jac->casadi_n_in = &ls_cost_nm5_n_in;
                ls_cost_jac->casadi_n_out = &ls_cost_nm5_n_out;
            }
            else
            {
                ls_cost_jac->casadi_fun = &ls_costN_nm5;
                ls_cost_jac->casadi_work = &ls_costN_nm5_work;
                ls_cost_jac->casadi_sparsity_in = &ls_costN_nm5_sparsity_in;
                ls_cost_jac->casadi_sparsity_out = &ls_costN_nm5_sparsity_out;
                ls_cost_jac->casadi_n_in = &ls_costN_nm5_n_in;
                ls_cost_jac->casadi_n_out = &ls_costN_nm5_n_out;
            }
            break;
        case 5:
            if (indx < N)
            {
                ls_cost_jac->casadi_fun = &ls_cost_nm6;
                ls_cost_jac->casadi_work = &ls_cost_nm6_work;
                ls_cost_jac->casadi_sparsity_in = &ls_cost_nm6_sparsity_in;
                ls_cost_jac->casadi_sparsity_out = &ls_cost_nm6_sparsity_out;
                ls_cost_jac->casadi_n_in = &ls_cost_nm6_n_in;
                ls_cost_jac->casadi_n_out = &ls_cost_nm6_n_out;
            }
            else
            {
                ls_cost_jac->casadi_fun = &ls_costN_nm6;
                ls_cost_jac->casadi_work = &ls_costN_nm6_work;
                ls_cost_jac->casadi_sparsity_in = &ls_costN_nm6_sparsity_in;
                ls_cost_jac->casadi_sparsity_out = &ls_costN_nm6_sparsity_out;
                ls_cost_jac->casadi_n_in = &ls_costN_nm6_n_in;
                ls_cost_jac->casadi_n_out = &ls_costN_nm6_n_out;
            }
            break;
        default:
            printf("Problem size not available\n");
            exit(1);
            break;
    }

    return;
}



void read_initial_state(const int nx, const int num_free_masses, double *x0)
{
    double *ptr;
    switch (num_free_masses)
    {
        case 1:
            ptr = x0_nm2;
            break;
        case 2:
            ptr = x0_nm3;
            break;
        case 3:
            ptr = x0_nm4;
            break;
        case 4:
            ptr = x0_nm5;
            break;
        case 5:
            ptr = x0_nm6;
            break;
        default:
            printf("\nwrong number of free masses\n");
            exit(1);
            break;
    }
    for (int i = 0; i < nx; i++)
        x0[i] = ptr[i];
}



void read_final_state(const int nx, const int num_free_masses, double *xN)
{
    double *ptr;
    switch (num_free_masses)
    {
        case 1:
            ptr = xN_nm2;
            break;
        case 2:
            ptr = xN_nm3;
            break;
        case 3:
            ptr = xN_nm4;
            break;
        case 4:
            ptr = xN_nm5;
            break;
        case 5:
            ptr = xN_nm6;
            break;
        default:
            printf("\nwrong number of free masses\n");
            exit(1);
            break;
    }
    for (int i = 0; i < nx; i++)
        xN[i] = ptr[i];
}



// hand-generated external function for externally provided hessian and gradient
void ext_cost_nm2(void *fun,
    ext_fun_arg_t *type_in,
    void **in,
    ext_fun_arg_t *type_out,
    void **out)
{

    int ii;

    int nu = 3;
    int nx = 6;

    int nv = nu+nx;

    // ref
    double *ref = (double *)calloc(nx+nu, sizeof(double));
    for (ii = 0; ii < nu; ii++)
        ref[ii] = 0.0;
    for (ii = 0; ii < nx; ii++)
        ref[nu+ii] = xN_nm2[ii];

    // Hessian
    double *hess = (double *)out[1];
    for (ii = 0; ii < nv*nv; ii++)
        hess[ii] = 0.0;
    for (ii = 0; ii < nu; ii++)
        hess[ii*(nv+1)] = 1.0;
    for (; ii < nu+nx; ii++)
        hess[ii*(nv+1)] = 1e-2;

    // gradient
    double *ux = (double *)in[0];
    double *grad = (double *)out[0];
    for (ii = 0; ii < nv; ii++)
        grad[ii] = 0.0;
    for (ii = 0; ii < nv; ii++)
        grad[ii] = hess[ii*(nv+1)] * (ux[ii] - ref[ii]);

    free(ref);
    return;

}

void ext_cost_nm3(void *fun,
    ext_fun_arg_t *type_in,
    void **in,
    ext_fun_arg_t *type_out,
    void **out)
{

    int ii;

    int nu = 3;
    int nx = 12;

    int nv = nu+nx;

    // ref
    double *ref = (double *)calloc(nx+nu, sizeof(double));
    for (ii = 0; ii < nu; ii++)
        ref[ii] = 0.0;
    for (ii = 0; ii < nx; ii++)
        ref[nu+ii] = xN_nm3[ii];

    // Hessian
    double *hess = (double *)out[1];
    for (ii = 0; ii < nv*nv; ii++)
        hess[ii] = 0.0;
    for (ii = 0; ii < nu; ii++)
        hess[ii*(nv+1)] = 1.0;
    for (; ii < nu+nx; ii++)
        hess[ii*(nv+1)] = 1e-2;

    // gradient
    double *ux = (double *)in[0];
    double *grad = (double *)out[0];
    for (ii = 0; ii < nv; ii++)
        grad[ii] = 0.0;
    for (ii = 0; ii < nv; ii++)
        grad[ii] = hess[ii*(nv+1)] * (ux[ii] - ref[ii]);

    free(ref);
    return;

}

void ext_cost_nm4(void *fun,
    ext_fun_arg_t *type_in,
    void **in,
    ext_fun_arg_t *type_out,
    void **out)
{

    int ii;

    int nu = 3;
    int nx = 18;

    int nv = nu+nx;

    // ref
    double *ref = (double *)calloc(nx+nu, sizeof(double));
    for (ii = 0; ii < nu; ii++)
        ref[ii] = 0.0;
    for (ii = 0; ii < nx; ii++)
        ref[nu+ii] = xN_nm4[ii];

    // Hessian
    double *hess = (double *)out[1];
    for (ii = 0; ii < nv*nv; ii++)
        hess[ii] = 0.0;
    for (ii = 0; ii < nu; ii++)
        hess[ii*(nv+1)] = 1.0;
    for (; ii < nu+nx; ii++)
        hess[ii*(nv+1)] = 1e-2;

    // gradient
    double *ux = (double *)in[0];
    double *grad = (double *)out[0];
    for (ii = 0; ii < nv; ii++)
        grad[ii] = 0.0;
    for (ii = 0; ii < nv; ii++)
        grad[ii] = hess[ii*(nv+1)] * (ux[ii] - ref[ii]);

    free(ref);
    return;

}

void ext_cost_nm5(void *fun,
    ext_fun_arg_t *type_in,
    void **in,
    ext_fun_arg_t *type_out,
    void **out)
{

    int ii;

    int nu = 3;
    int nx = 24;

    int nv = nu+nx;

    // ref
    double *ref = (double *)calloc(nx+nu, sizeof(double));
    for (ii = 0; ii < nu; ii++)
        ref[ii] = 0.0;
    for (ii = 0; ii < nx; ii++)
        ref[nu+ii] = xN_nm5[ii];

    // Hessian
    double *hess = (double *)out[1];
    for (ii = 0; ii < nv*nv; ii++)
        hess[ii] = 0.0;
    for (ii = 0; ii < nu; ii++)
        hess[ii*(nv+1)] = 1.0;
    for (; ii < nu+nx; ii++)
        hess[ii*(nv+1)] = 1e-2;

    // gradient
    double *ux = (double *)in[0];
    double *grad = (double *)out[0];
    for (ii = 0; ii < nv; ii++)
        grad[ii] = 0.0;
    for (ii = 0; ii < nv; ii++)
        grad[ii] = hess[ii*(nv+1)] * (ux[ii] - ref[ii]);

    free(ref);
    return;

}

void ext_cost_nm6(void *fun,
    ext_fun_arg_t *type_in,
    void **in,
    ext_fun_arg_t *type_out,
    void **out)
{

    int ii;

    int nu = 3;
    int nx = 30;

    int nv = nu+nx;

    // ref
    double *ref = (double *)calloc(nx+nu, sizeof(double));
    for (ii = 0; ii < nu; ii++)
        ref[ii] = 0.0;
    for (ii = 0; ii < nx; ii++)
        ref[nu+ii] = xN_nm6[ii];

    // Hessian
    double *hess = (double *)out[1];
    for (ii = 0; ii < nv*nv; ii++)
        hess[ii] = 0.0;
    for (ii = 0; ii < nu; ii++)
        hess[ii*(nv+1)] = 1.0;
    for (; ii < nu+nx; ii++)
        hess[ii*(nv+1)] = 1e-2;

    // gradient
    double *ux = (double *)in[0];
    double *grad = (double *)out[0];
    for (ii = 0; ii < nv; ii++)
        grad[ii] = 0.0;
    for (ii = 0; ii < nv; ii++)
        grad[ii] = hess[ii*(nv+1)] * (ux[ii] - ref[ii]);

    free(ref);
    return;

}



// hand-wirtten box constraints on states as nonlinear constraints
void nonlin_constr_nm2(void *evaluate,
    ext_fun_arg_t *type_in,
    void **in,
    ext_fun_arg_t *type_out,
    void **out)
{

    int ii;

    int nu = 3;
    int nx = 6;

    int nh = nx;

    // fun
    struct blasfeo_dvec_args *fun_args = (struct blasfeo_dvec_args *)out[0];
    struct blasfeo_dvec *fun = fun_args->x;
    int xi = fun_args->xi;
    struct blasfeo_dvec *ux = (struct blasfeo_dvec *)in[0];
    blasfeo_dveccp(nx, ux, nu, fun, xi);

    // jacobian
    struct blasfeo_dmat_args *jac_args = (struct blasfeo_dmat_args *)out[1];
    struct blasfeo_dmat *jac = jac_args->A;
    int ai = jac_args->ai;
    int aj = jac_args->aj;
    blasfeo_dgese(nu+nx, nh, 0.0, jac, ai, aj);
    for (ii = 0; ii < nh; ii++)
        BLASFEO_DMATEL(jac, ai+nu+ii, aj+ii) = 1.0;

    return;

}

void nonlin_constr_nm3(void *evaluate,
    ext_fun_arg_t *type_in,
    void **in,
    ext_fun_arg_t *type_out,
    void **out)
{

    int ii;

    int nu = 3;
    int nx = 12;

    int nh = nx;

    // fun
    struct blasfeo_dvec_args *fun_args = (struct blasfeo_dvec_args *)out[0];
    struct blasfeo_dvec *fun = fun_args->x;
    int xi = fun_args->xi;
    struct blasfeo_dvec *ux = (struct blasfeo_dvec *)in[0];
    blasfeo_dveccp(nx, ux, nu, fun, xi);

    // jacobian
    struct blasfeo_dmat_args *jac_args = (struct blasfeo_dmat_args *)out[1];
    struct blasfeo_dmat *jac = jac_args->A;
    int ai = jac_args->ai;
    int aj = jac_args->aj;
    blasfeo_dgese(nu+nx, nh, 0.0, jac, ai, aj);
    for (ii = 0; ii < nh; ii++)
        BLASFEO_DMATEL(jac, ai+nu+ii, aj+ii) = 1.0;

    return;

}

void nonlin_constr_nm4(void *evaluate,
    ext_fun_arg_t *type_in,
    void **in,
    ext_fun_arg_t *type_out,
    void **out)
{

    int ii;

    int nu = 3;
    int nx = 18;

    int nh = nx;

    // fun
    struct blasfeo_dvec_args *fun_args = (struct blasfeo_dvec_args *)out[0];
    struct blasfeo_dvec *fun = fun_args->x;
    int xi = fun_args->xi;
    struct blasfeo_dvec *ux = (struct blasfeo_dvec *)in[0];
    blasfeo_dveccp(nx, ux, nu, fun, xi);

    // jacobian
    struct blasfeo_dmat_args *jac_args = (struct blasfeo_dmat_args *)out[1];
    struct blasfeo_dmat *jac = jac_args->A;
    int ai = jac_args->ai;
    int aj = jac_args->aj;
    blasfeo_dgese(nu+nx, nh, 0.0, jac, ai, aj);
    for (ii = 0; ii < nh; ii++)
        BLASFEO_DMATEL(jac, ai+nu+ii, aj+ii) = 1.0;

    return;

}

void nonlin_constr_nm5(void *evaluate,
    ext_fun_arg_t *type_in,
    void **in,
    ext_fun_arg_t *type_out,
    void **out)
{

    int ii;

    int nu = 3;
    int nx = 24;

    int nh = nx;

    // fun
    struct blasfeo_dvec_args *fun_args = (struct blasfeo_dvec_args *)out[0];
    struct blasfeo_dvec *fun = fun_args->x;
    int xi = fun_args->xi;
    struct blasfeo_dvec *ux = (struct blasfeo_dvec *)in[0];
    blasfeo_dveccp(nx, ux, nu, fun, xi);

    // jacobian
    struct blasfeo_dmat_args *jac_args = (struct blasfeo_dmat_args *)out[1];
    struct blasfeo_dmat *jac = jac_args->A;
    int ai = jac_args->ai;
    int aj = jac_args->aj;
    blasfeo_dgese(nu+nx, nh, 0.0, jac, ai, aj);
    for (ii = 0; ii < nh; ii++)
        BLASFEO_DMATEL(jac, ai+nu+ii, aj+ii) = 1.0;

    return;

}

void nonlin_constr_nm6(void *evaluate,
    ext_fun_arg_t *type_in,
    void **in,
    ext_fun_arg_t *type_out,
    void **out)
{

    int ii;

    int nu = 3;
    int nx = 30;

    int nh = nx;

    // fun
    struct blasfeo_dvec_args *fun_args = (struct blasfeo_dvec_args *)out[0];
    struct blasfeo_dvec *fun = fun_args->x;
    int xi = fun_args->xi;
    struct blasfeo_dvec *ux = (struct blasfeo_dvec *)in[0];
    blasfeo_dveccp(nx, ux, nu, fun, xi);

    // jacobian
    struct blasfeo_dmat_args *jac_args = (struct blasfeo_dmat_args *)out[1];
    struct blasfeo_dmat *jac = jac_args->A;
    int ai = jac_args->ai;
    int aj = jac_args->aj;
    blasfeo_dgese(nu+nx, nh, 0.0, jac, ai, aj);
    for (ii = 0; ii < nh; ii++)
        BLASFEO_DMATEL(jac, ai+nu+ii, aj+ii) = 1.0;

    return;

}



void setup_and_solve_nlp(int NN,
    int NMF,
    std::string const& con_str,
    std::string const& cost_str,
    std::string const& qp_solver_str,
    std::string const& model_str,
    std::string const& integrator_str
    )
{
    /************************************************
    * problem dimensions
    ************************************************/

    int  *nx = (int *)calloc(NN+1, sizeof(int));
    int  *nu = (int *)calloc(NN+1, sizeof(int));
    int *nbx = (int *)calloc(NN+1, sizeof(int));
    int *nbu = (int *)calloc(NN+1, sizeof(int));
    int  *nb = (int *)calloc(NN+1, sizeof(int));
    int  *ng = (int *)calloc(NN+1, sizeof(int));
    int  *nh = (int *)calloc(NN+1, sizeof(int));
    int  *nq = (int *)calloc(NN+1, sizeof(int));
    int  *ns = (int *)calloc(NN+1, sizeof(int));
    int  *ny = (int *)calloc(NN+1, sizeof(int));
    int  *nz = (int *)calloc(NN+1, sizeof(int));


    int NX = 6 * NMF;
    int NU = 3;

    nx[0] = NX;
    nu[0] = NU;
    ny[0] = nx[0]+nu[0];
    nz[0] = 0;

    constraints_t con_type = constraints_enum(con_str);
    switch (con_type)
    {
        case BOX:
            nbx[0] = nx[0];
            nbu[0] = nu[0];
            nb[0] = nbu[0]+nbx[0];
            ng[0] = 0;
            nh[0] = 0;
            break;
        case GENERAL:
            nbx[0] = 0;
            nbu[0] = 0;
            nb[0] = 0;
            ng[0] = nu[0]+nx[0];
            nh[0] = 0;
            break;
        case GENERAL_NONLINEAR:
        default:
            nbx[0] = 0;
            nbu[0] = 0;
            nb[0] = 0;
            ng[0] = nu[0];
            nh[0] = nx[0];
            break;
    }

    for (int i = 1; i < NN; i++)
    {
        nx[i] = NX;
        nu[i] = NU;
        nbx[i] = NMF;
        nbu[i] = NU;
        nb[i] = nbu[i]+nbx[i];
        ng[i] = 0;
        nh[i] = 0;
        ny[i] = nx[i]+nu[i];
        nz[i] = 0;
    }

    nx[NN] = NX;
    nu[NN] = 0;
    nbx[NN] = NX;
    nbu[NN] = 0;
    nb[NN] = nbu[NN]+nbx[NN];
    ng[NN] = 0;
    nh[NN] = 0;
    ny[NN] = nx[NN]+nu[NN];
    nz[NN] = 0;

    /************************************************
    * problem data
    ************************************************/

    double wall_pos = -0.01;
    double UMAX = 10;

    double x_pos_inf = +1e4;
    double x_neg_inf = -1e4;

    double *xref = (double *)malloc(NX*sizeof(double));
    read_final_state(NX, NMF, xref);

    double uref[3] = {0.0, 0.0, 0.0};

    double *diag_cost_x = (double *)malloc(NX*sizeof(double));

    for (int i = 0; i < NX; i++)
        diag_cost_x[i] = 1e-2;

    double diag_cost_u[3] = {1.0, 1.0, 1.0};


    // idxb0
    int *idxb0 = (int *)malloc(nb[0]*sizeof(int));

    for (int i = 0; i < nb[0]; i++) idxb0[i] = i;

    // idxb1
    int *idxb1 = (int *)malloc(nb[1]*sizeof(int));
    for (int i = 0; i < NU; i++) idxb1[i] = i;

    for (int i = 0; i < NMF; i++) idxb1[NU+i] = NU + 6*i + 1;

    // idxbN
    int *idxbN = (int *)malloc(nb[NN]*sizeof(int));
    for (int i = 0; i < nb[NN]; i++)
        idxbN[i] = i;

    // lb0, ub0
    double *lb0 = (double *)malloc((NX+NU)*sizeof(double));
    double *ub0 = (double *)malloc((NX+NU)*sizeof(double));

    for (int i = 0; i < NU; i++)
    {
        lb0[i] = -UMAX;
        ub0[i] = +UMAX;
    }
    read_initial_state(NX, NMF, lb0+NU);
    read_initial_state(NX, NMF, ub0+NU);

    // lb1, ub1
    double *lb1 = (double *)malloc((NMF+NU)*sizeof(double));
    double *ub1 = (double *)malloc((NMF+NU)*sizeof(double));

    for (int j = 0; j < NU; j++)
    {
        lb1[j] = -UMAX;  // umin
        ub1[j] = +UMAX;  // umax
    }
    for (int j = 0; j < NMF; j++)
    {
        lb1[NU+j] = wall_pos;  // wall position
        ub1[NU+j] = x_pos_inf;
    }

    // lbN, ubN
    double *lbN = (double *)malloc(NX*sizeof(double));
    double *ubN = (double *)malloc(NX*sizeof(double));

    for (int i = 0; i < NX; i++)
    {
        lbN[i] = x_neg_inf;
        ubN[i] = x_pos_inf;
    }

    /************************************************
    * plan + config
    ************************************************/

    ocp_nlp_solver_plan *plan = ocp_nlp_plan_create(NN);

    // TODO(dimitris): not necessarily GN, depends on cost module
    plan->nlp_solver = SQP;

    ocp_nlp_cost_t cost_type = cost_enum(cost_str);
    switch (cost_type)
    {
        case LINEAR_LS:
            for (int i = 0; i <= NN; i++)
            {
                plan->nlp_cost[i] = LINEAR_LS;
            }
            break;
        case NONLINEAR_LS:
            for (int i = 0; i <= NN; i++)
            {
                plan->nlp_cost[i] = NONLINEAR_LS;
            }
            break;
        case EXTERNALLY_PROVIDED:
            for (int i = 0; i < NN; i++)
            {
                plan->nlp_cost[i] = EXTERNALLY_PROVIDED;
            }
            plan->nlp_cost[NN] = LINEAR_LS;
            break;
        default:
            for (int i = 0; i <= NN; i++)
            {
                if (i%3 == 0)
                    plan->nlp_cost[i] = EXTERNALLY_PROVIDED;
                else if (i%3 == 1)
                    plan->nlp_cost[i] = LINEAR_LS;
                else if (i%3 == 2)
                    plan->nlp_cost[i] = NONLINEAR_LS;
            }
            plan->nlp_cost[NN] = LINEAR_LS;
            break;
    }

    ocp_qp_solver_t qp_solver_type = qp_solver_enum(qp_solver_str);
    plan->ocp_qp_solver_plan.qp_solver = qp_solver_type;

    ocp_nlp_dynamics_t model_type = nlp_dynamics_enum(model_str);
    switch (model_type)
    {
        case CONTINUOUS_MODEL:
            for (int i = 0; i < NN; i++)
            {
                plan->nlp_dynamics[i] = CONTINUOUS_MODEL;
            }
            break;
        case DISCRETE_MODEL:
            for (int i = 0; i < NN; i++)
            {
                plan->nlp_dynamics[i] = DISCRETE_MODEL;
            }
            break;
        default:
            for (int i = 0; i < NN; i++)
            {
                if (i < NN/2)
                {
                    plan->nlp_dynamics[i] = CONTINUOUS_MODEL;
                }
                else
                {
                    plan->nlp_dynamics[i] = DISCRETE_MODEL;
                }
            }
            break;
    }

    sim_solver_t integrator_type = integrator_enum(integrator_str);
    if (model_type != DISCRETE_MODEL)
    {
        switch (integrator_type)
        {
            case IRK:
                for (int i = 0; i < NN; i++)
                {
                    if (plan->nlp_dynamics[i] == CONTINUOUS_MODEL)
                        plan->sim_solver_plan[i].sim_solver = IRK;
                }
                break;
            case ERK:
                for (int i = 0; i < NN; i++)
                {
                    if (plan->nlp_dynamics[i] == CONTINUOUS_MODEL)
                        plan->sim_solver_plan[i].sim_solver = ERK;
                }
                break;
            case LIFTED_IRK:
                for (int i = 0; i < NN; i++)
                {
                    if (plan->nlp_dynamics[i] == CONTINUOUS_MODEL)
                        plan->sim_solver_plan[i].sim_solver = LIFTED_IRK;
                }
                break;
            default:
                for (int i = 0; i < NN; i++)
                {
                    if (plan->nlp_dynamics[i] == CONTINUOUS_MODEL)
                    {
                        if (i%4 == 0)
                            plan->sim_solver_plan[i].sim_solver = IRK;
                        else if (i%4 == 1)
                            plan->sim_solver_plan[i].sim_solver = ERK;
                        else if (i%4 == 2)
                            plan->sim_solver_plan[i].sim_solver = IRK;
                        else if (i%4 == 3)
                            plan->sim_solver_plan[i].sim_solver = LIFTED_IRK;
                    }
                }
                break;
        }
    }

    for (int i = 0; i <= NN; i++)
    {
        plan->nlp_constraints[i] = BGH;
    }

    if (NMF > 3 && model_type != CONTINUOUS_MODEL)
        return;

    // TODO(dimitris): fix minor memory leak
    // here
    ocp_nlp_solver_config *config = ocp_nlp_config_create(*plan);

    /************************************************
    * ocp_nlp_dims
    ************************************************/

    ocp_nlp_dims *dims = ocp_nlp_dims_create(config);

    ocp_nlp_dims_set_opt_vars(config, dims, "nx", nx);
    ocp_nlp_dims_set_opt_vars(config, dims, "nu", nu);
    ocp_nlp_dims_set_opt_vars(config, dims, "nz", nz);
    ocp_nlp_dims_set_opt_vars(config, dims, "ns", ns);

    for (int i = 0; i <= NN; i++)
    {
        if (plan->nlp_cost[i] != EXTERNALLY_PROVIDED)
        {
            ocp_nlp_dims_set_cost(config, dims, i, "ny", &ny[i]);
        }
        ocp_nlp_dims_set_constraints(config, dims, i, "nbx", &nbx[i]);
        ocp_nlp_dims_set_constraints(config, dims, i, "nbu", &nbu[i]);
        ocp_nlp_dims_set_constraints(config, dims, i, "ng", &ng[i]);
        ocp_nlp_dims_set_constraints(config, dims, i, "nh", &nh[i]);
        // ocp_nlp_dims_set_constraints(config, dims, i, "np", &nq[i]);
    }

    /************************************************
    * dynamics
    ************************************************/

    // explicit
    external_function_casadi *expl_vde_for = (external_function_casadi *)
                                              malloc(NN*sizeof(external_function_casadi));

    // implicit
    external_function_casadi *impl_ode_fun = (external_function_casadi *)
                                              malloc(NN*sizeof(external_function_casadi));
    external_function_casadi *impl_ode_fun_jac_x_xdot = (external_function_casadi *)
                                                        malloc(NN*sizeof(external_function_casadi));
    external_function_casadi *impl_ode_fun_jac_x_xdot_u = (external_function_casadi *)
                                                        malloc(NN*sizeof(external_function_casadi));
    external_function_casadi *impl_ode_jac_x_xdot_u = (external_function_casadi *)
                                                        malloc(NN*sizeof(external_function_casadi));

    // discrete model
    external_function_casadi *erk4_casadi = NULL;

    if (NMF < 4)
        erk4_casadi = (external_function_casadi *)malloc(NN*sizeof(external_function_casadi));

    select_dynamics_casadi(NN, NMF, expl_vde_for, impl_ode_fun,
            impl_ode_fun_jac_x_xdot, impl_ode_fun_jac_x_xdot_u, impl_ode_jac_x_xdot_u, erk4_casadi);

    // forw_vde
    external_function_casadi_create_array(NN, expl_vde_for);
    // impl_ode
    external_function_casadi_create_array(NN, impl_ode_fun);
    //
    external_function_casadi_create_array(NN, impl_ode_fun_jac_x_xdot);
    //
    external_function_casadi_create_array(NN, impl_ode_fun_jac_x_xdot_u);
    //
    external_function_casadi_create_array(NN, impl_ode_jac_x_xdot_u);

    if (erk4_casadi != NULL)
        external_function_casadi_create_array(NN, erk4_casadi);

    /************************************************
    * nonlinear least squares
    ************************************************/

    external_function_casadi *ls_cost_jac_casadi = (external_function_casadi *)
                                                    malloc((NN+1)*sizeof(external_function_casadi));
    external_function_generic *ext_cost_generic = (external_function_generic *)
                                                   malloc(NN*sizeof(external_function_casadi));

    for (int i = 0; i <= NN; i++)
    {
        switch (plan->nlp_cost[i])
        {
            case LINEAR_LS:
                // do nothing
                break;

            case NONLINEAR_LS:
                select_ls_stage_cost_jac_casadi(i, NN, NMF, &ls_cost_jac_casadi[i]);
                external_function_casadi_create(&ls_cost_jac_casadi[i]);
                break;

            case EXTERNALLY_PROVIDED:
                // TODO(dimitris): move
                // inside
                // select_ls_stage_cost_jac_casadi?
                switch (NMF)
                {
                    case 1:
                        ext_cost_generic[i].evaluate = &ext_cost_nm2;
                        break;
                    case 2:
                        ext_cost_generic[i].evaluate = &ext_cost_nm3;
                        break;
                    case 3:
                        ext_cost_generic[i].evaluate = &ext_cost_nm4;
                        break;
                    case 4:
                        ext_cost_generic[i].evaluate = &ext_cost_nm5;
                        break;
                    case 5:
                        ext_cost_generic[i].evaluate = &ext_cost_nm6;
                        break;
                    default:
                        printf("\next cost not implemented for this numer of masses\n\n");
                        exit(1);
                }
                break;

            default:
                printf("\ncost not correctly specified\n\n");
                exit(1);
        }
    }

    /************************************************
    * nonlinear constraints
    ************************************************/
    external_function_generic nonlin_constr_generic;

    if (con_type == GENERAL_NONLINEAR)
    {
        switch (NMF)
        {
            case 1:
                nonlin_constr_generic.evaluate = &nonlin_constr_nm2;
                break;
            case 2:
                nonlin_constr_generic.evaluate = &nonlin_constr_nm3;
                break;
            case 3:
                nonlin_constr_generic.evaluate = &nonlin_constr_nm4;
                break;
            case 4:
                nonlin_constr_generic.evaluate = &nonlin_constr_nm5;
                break;
            case 5:
                nonlin_constr_generic.evaluate = &nonlin_constr_nm6;
                break;
            default:
                printf("\nnonlin constr not implemented for this number of masses\n\n");
                exit(1);
        }
    }


    /************************************************
    * nlp_in
    ************************************************/

    ocp_nlp_in *nlp_in = ocp_nlp_in_create(config, dims);

    // sampling times
    for (int ii = 0; ii < NN; ii++)
        nlp_in->Ts[ii] = TF/NN;

    // output definition: y = [x; u]
    /* cost */
    ocp_nlp_cost_ls_model *stage_cost_ls;
    ocp_nlp_cost_nls_model *stage_cost_nls;
    ocp_nlp_cost_external_model *stage_cost_external;

    for (int i = 0; i < NN; i++)
    {
        switch (plan->nlp_cost[i])
        {
            case LINEAR_LS:

                stage_cost_ls = (ocp_nlp_cost_ls_model *) nlp_in->cost[i];

                // Cyt
                blasfeo_dgese(nu[i]+nx[i], ny[i], 0.0, &stage_cost_ls->Cyt, 0, 0);
                for (int j = 0; j < nu[i]; j++)
                    BLASFEO_DMATEL(&stage_cost_ls->Cyt, j, nx[i]+j) = 1.0;
                for (int j = 0; j < nx[i]; j++)
                    BLASFEO_DMATEL(&stage_cost_ls->Cyt, nu[i]+j, j) = 1.0;

                // W
                blasfeo_dgese(ny[i], ny[i], 0.0, &stage_cost_ls->W, 0, 0);
                for (int j = 0; j < nx[i]; j++)
                    BLASFEO_DMATEL(&stage_cost_ls->W, j, j) = diag_cost_x[j];
                for (int j = 0; j < nu[i]; j++)
                    BLASFEO_DMATEL(&stage_cost_ls->W, nx[i]+j, nx[i]+j) = diag_cost_u[j];

                // y_ref
                blasfeo_pack_dvec(nx[i], xref, &stage_cost_ls->y_ref, 0);
                blasfeo_pack_dvec(nu[i], uref, &stage_cost_ls->y_ref, nx[i]);
                break;

            case NONLINEAR_LS:

                stage_cost_nls = (ocp_nlp_cost_nls_model *) nlp_in->cost[i];

                // nls_jac
                stage_cost_nls->nls_jac = (external_function_generic *) &ls_cost_jac_casadi[i];

                // W
                blasfeo_dgese(ny[i], ny[i], 0.0, &stage_cost_nls->W, 0, 0);
                for (int j = 0; j < nx[i]; j++)
                    BLASFEO_DMATEL(&stage_cost_nls->W, j, j) = diag_cost_x[j];
                for (int j = 0; j < nu[i]; j++)
                    BLASFEO_DMATEL(&stage_cost_nls->W, nx[i]+j, nx[i]+j) = diag_cost_u[j];

                // y_ref
                blasfeo_pack_dvec(nx[i], xref, &stage_cost_nls->y_ref, 0);
                blasfeo_pack_dvec(nu[i], uref, &stage_cost_nls->y_ref, nx[i]);
                break;

            case EXTERNALLY_PROVIDED:

                stage_cost_external = (ocp_nlp_cost_external_model *) nlp_in->cost[i];
                stage_cost_external->ext_cost = &ext_cost_generic[i];
                assert(i < NN && "externally provided cost not implemented for last stage!");
                break;

            default:
                printf("\ncost not correctly specified\n\n");
                exit(1);
        }
    }



    /* dynamics */
    int set_fun_status;

    // TODO(dimitris): remove after setting
    // function via nlp interface
    ocp_nlp_dynamics_disc_model *dynamics;

    for (int i = 0; i < NN; i++)
    {
        switch (plan->nlp_dynamics[i])
        {
            case CONTINUOUS_MODEL:

                if (plan->sim_solver_plan[i].sim_solver == ERK)
                {
                    set_fun_status = ocp_nlp_dynamics_model_set(config, nlp_in, i,
                                                            "expl_vde_for", &expl_vde_for[i]);
                    if (set_fun_status != 0) exit(1);
                }
                else if (plan->sim_solver_plan[i].sim_solver == IRK)
                {
                    set_fun_status = ocp_nlp_dynamics_model_set(config, nlp_in, i,
                                                            "impl_ode_fun", &impl_ode_fun[i]);
                    if (set_fun_status != 0) exit(1);
                    set_fun_status = ocp_nlp_dynamics_model_set(config, nlp_in, i,
                                            "impl_ode_fun_jac_x_xdot", &impl_ode_fun_jac_x_xdot[i]);
                    if (set_fun_status != 0) exit(1);
                    set_fun_status = ocp_nlp_dynamics_model_set(config, nlp_in, i,
                                                "impl_ode_jac_x_xdot_u", &impl_ode_jac_x_xdot_u[i]);
                    if (set_fun_status != 0) exit(1);
                }
                else if (plan->sim_solver_plan[i].sim_solver == LIFTED_IRK)
                {
                    set_fun_status = ocp_nlp_dynamics_model_set(config, nlp_in, i,
                                                            "impl_ode_fun", &impl_ode_fun[i]);
                    if (set_fun_status != 0) exit(1);
                    set_fun_status = ocp_nlp_dynamics_model_set(config, nlp_in, i,
                                        "impl_ode_fun_jac_x_xdot_u", &impl_ode_fun_jac_x_xdot_u[i]);
                    if (set_fun_status != 0) exit(1);
                }
                break;
            case DISCRETE_MODEL:
                // TODO(dimitris): do this
                // through the interface and
                // remove header
                if (NMF < 4)
                {
                    dynamics = (ocp_nlp_dynamics_disc_model *)nlp_in->dynamics[i];
                    dynamics->discrete_model = (external_function_generic *) &erk4_casadi[i];
                }
                break;

            default:
                printf("\ndynamics not correctly specified\n\n");
                exit(1);
        }
    }

    /* constraints */
    ocp_nlp_constraints_bgh_model **constraints =
        (ocp_nlp_constraints_bgh_model **) nlp_in->constraints;

    // first stage
    switch (con_type)
    {
        case BOX:
            blasfeo_pack_dvec(nb[0], lb0, &constraints[0]->d, 0);
            blasfeo_pack_dvec(nb[0], ub0, &constraints[0]->d, nb[0]+ng[0]);
            constraints[0]->idxb = idxb0;
            break;
        case GENERAL:
            double *Cu0; d_zeros(&Cu0, ng[0], nu[0]);
            for (int ii = 0; ii < nu[0]; ii++)
                Cu0[ii*(ng[0]+1)] = 1.0;

            double *Cx0; d_zeros(&Cx0, ng[0], nx[0]);
            for (int ii = 0; ii < nx[0]; ii++)
                Cx0[nu[0]+ii*(ng[0]+1)] = 1.0;

            blasfeo_pack_tran_dmat(ng[0], nu[0], Cu0, ng[0], &constraints[0]->DCt, 0, 0);
            blasfeo_pack_tran_dmat(ng[0], nx[0], Cx0, ng[0], &constraints[0]->DCt, nu[0], 0);
            blasfeo_pack_dvec(ng[0], lb0, &constraints[0]->d, nb[0]);
            blasfeo_pack_dvec(ng[0], ub0, &constraints[0]->d, 2*nb[0]+ng[0]);

            d_free(Cu0);
            d_free(Cx0);
            break;
        case GENERAL_NONLINEAR:
        default:
            blasfeo_dgese(nu[0]+nx[0], ng[0], 0.0, &constraints[0]->DCt, 0, 0);
            for (int ii = 0; ii < ng[0]; ii++)
                BLASFEO_DMATEL(&constraints[0]->DCt, ii, ii) = 1.0;

            ocp_nlp_constraints_bgh_model **nl_constr = (ocp_nlp_constraints_bgh_model **)
                                                    nlp_in->constraints;
            nl_constr[0]->h = &nonlin_constr_generic;

            blasfeo_pack_dvec(ng[0]+nh[0], lb0, &constraints[0]->d, nb[0]);
            blasfeo_pack_dvec(ng[0]+nh[0], ub0, &constraints[0]->d, 2*nb[0]+ng[0]+nh[0]);
            break;
    }

    // other stages
    for (int i = 1; i < NN; i++)
    {
        blasfeo_pack_dvec(nb[i], lb1, &constraints[i]->d, 0);
        blasfeo_pack_dvec(nb[i], ub1, &constraints[i]->d, nb[i]+ng[i]);
        constraints[i]->idxb = idxb1;
    }
    blasfeo_pack_dvec(nb[NN], lbN, &constraints[NN]->d, 0);
    blasfeo_pack_dvec(nb[NN], ubN, &constraints[NN]->d, nb[NN]+ng[NN]);
    constraints[NN]->idxb = idxbN;

    /************************************************
    * sqp opts
    ************************************************/

    void *nlp_opts = ocp_nlp_opts_create(config, dims);
    ocp_nlp_sqp_opts *sqp_opts = (ocp_nlp_sqp_opts *) nlp_opts;

    for (int i = 0; i < NN; ++i)
    {
        if (plan->nlp_dynamics[i] == CONTINUOUS_MODEL)
        {
            ocp_nlp_dynamics_cont_opts *dynamics_stage_opts = (ocp_nlp_dynamics_cont_opts *)
                                                              sqp_opts->dynamics[i];
            sim_rk_opts *sim_opts = (sim_rk_opts *)dynamics_stage_opts->sim_solver;

            if (plan->sim_solver_plan[i].sim_solver == ERK)
            {
                sim_opts->ns = 4;
            }
            else if (plan->sim_solver_plan[i].sim_solver == IRK)
            {
                sim_opts->ns = 2;
                sim_opts->jac_reuse = true;
            }
        }
    }
    int maxIter = MAX_SQP_ITERS;
    double min_res_g = 1e-6;
    double min_res_b = 1e-6;
    double min_res_d = 1e-6;
    double min_res_m = 1e-6;

    ocp_nlp_opts_set(config, nlp_opts, "maxIter", &maxIter);
    ocp_nlp_opts_set(config, nlp_opts, "min_res_g", &min_res_g);
    ocp_nlp_opts_set(config, nlp_opts, "min_res_b", &min_res_b);
    ocp_nlp_opts_set(config, nlp_opts, "min_res_d", &min_res_d);
    ocp_nlp_opts_set(config, nlp_opts, "min_res_m", &min_res_m);

    /************************************************
    * ocp_nlp out
    ************************************************/

    ocp_nlp_out *nlp_out = ocp_nlp_out_create(config, dims);

    ocp_nlp_solver *solver = ocp_nlp_create(config, dims, nlp_opts);

    /************************************************
    * sqp solve
    ************************************************/

    int status;

    // warm start output initial guess of
    // solution
    for (int i=0; i <= NN; i++)
    {
        blasfeo_pack_dvec(nu[i], uref, nlp_out->ux+i, 0);
        blasfeo_pack_dvec(nx[i], xref, nlp_out->ux+i, nu[i]);
    }

    // call nlp solver
    status = ocp_nlp_solve(solver, nlp_in, nlp_out);

    double max_res = 0.0;
    double inf_norm_res_g = ((ocp_nlp_sqp_memory *)solver->mem)->nlp_res->inf_norm_res_g;
    double inf_norm_res_b = ((ocp_nlp_sqp_memory *)solver->mem)->nlp_res->inf_norm_res_b;
    double inf_norm_res_d = ((ocp_nlp_sqp_memory *)solver->mem)->nlp_res->inf_norm_res_d;
    double inf_norm_res_m = ((ocp_nlp_sqp_memory *)solver->mem)->nlp_res->inf_norm_res_m;
    max_res = (inf_norm_res_g > max_res) ? inf_norm_res_g : max_res;
    max_res = (inf_norm_res_b > max_res) ? inf_norm_res_b : max_res;
    max_res = (inf_norm_res_d > max_res) ? inf_norm_res_d : max_res;
    max_res = (inf_norm_res_m > max_res) ? inf_norm_res_m : max_res;

    std::cout << "max residuals: " << max_res << std::endl;
    REQUIRE(status == 0);
    REQUIRE(max_res <= TOL);

    /************************************************
    * free memory
    ************************************************/

    // TODO(dimitris): VALGRIND!
    external_function_casadi_free(expl_vde_for);
    free(expl_vde_for);

    external_function_casadi_free(impl_ode_fun);
    external_function_casadi_free(impl_ode_fun_jac_x_xdot);
    external_function_casadi_free(impl_ode_fun_jac_x_xdot_u);
    external_function_casadi_free(impl_ode_jac_x_xdot_u);

    free(impl_ode_fun);
    free(impl_ode_fun_jac_x_xdot);
    free(impl_ode_fun_jac_x_xdot_u);
    free(impl_ode_jac_x_xdot_u);

    if (erk4_casadi != NULL)
    {
        external_function_casadi_free(erk4_casadi);
        free(erk4_casadi);
    }

    ocp_nlp_opts_free(nlp_opts);
    ocp_nlp_in_free(nlp_in);
    ocp_nlp_out_free(nlp_out);
    ocp_nlp_free(solver);
    ocp_nlp_dims_free(dims);
    ocp_nlp_config_free(plan, config);

    free(xref);
    free(diag_cost_x);
    free(lb0);
    free(ub0);
    free(lb1);
    free(ub1);
    free(lbN);
    free(ubN);
    free(idxb0);
    free(idxb1);
    free(idxbN);

    for (int i = 0; i <= NN; i++)
    {
        switch (plan->nlp_cost[i])
        {
            case NONLINEAR_LS:
                external_function_casadi_free(&ls_cost_jac_casadi[i]);
                break;
            default:
                break;
        }
    }

    free(ls_cost_jac_casadi);
    free(ext_cost_generic);

    free(plan);

    free(nx);
    free(nu);
    free(nbx);
    free(nbu);
    free(nb);
    free(ng);
    free(nh);
    free(nq);
    free(ns);
    free(ny);
}



extern bool CHAIN_EXTENSIVE;

/************************************************
* TEST CASE: nonlinear chain
************************************************/

TEST_CASE("chain example", "[NLP solver]")
{
    std::vector<int> horizon_lenghts;
    std::vector<int> num_masses;
    std::vector<std::string> cons;
    std::vector<std::string> models;
    std::vector<std::string> integrators;
    std::vector<std::string> costs;
    std::vector<std::string> qp_solvers = { "SPARSE_HPIPM",
                                            // "SPARSE_HPMPC",
                                            // "SPARSE_QPDUNES",
                                            "DENSE_HPIPM",
                                            "DENSE_QPOASES"
#ifdef ACADOS_WITH_OOQP
                                            // , "DENSE_OOQP"
                                            // , "SPARSE_OOQP"
#endif
#ifdef ACADOS_WITH_QORE
                                            , "DENSE_QORE"
                                            };
#else
    };
#endif

    if (CHAIN_EXTENSIVE)
    {
        horizon_lenghts = {20};
        num_masses = {2, 3, 4};
        cons = {"BOX", "GENERAL", "NONLINEAR+GENERAL"};
        models = {"DISCRETE", "CONTINUOUS", "MIXED"};
        integrators = {"MIXED"};
        costs = {"MIXED"};
    }
    else
    {
        horizon_lenghts = {20, 25};
        num_masses = {2, 3, 4};
        cons = {"NONLINEAR+GENERAL"};
        models = {"DISCRETE", "CONTINUOUS", "MIXED"};
        integrators = {"MIXED"};
        costs = {"MIXED"};
    }

    for (int NN : horizon_lenghts)
    {
        SECTION("Horizon length: " + std::to_string(NN))
        {
            for (int NMF : num_masses)
            {
                SECTION("Number of masses: " + std::to_string(NMF))
                {
                    for (std::string con_str : cons)
                    {
                        SECTION("Type of constraints: " + con_str)
                        {
                            for (std::string cost_str : costs)
                            {
                                SECTION("Stage cost type: " + cost_str)
                                {
                                    for (std::string qp_solver_str : qp_solvers)
                                    {
                                        SECTION("QP solver: " + qp_solver_str)
                                        {
                                            for (std::string model_str : models)
                                            {
                                                SECTION("Type of model: " + model_str)
                                                {
                                                    for (std::string integrator_str : integrators)
                                                    {
                                                        SECTION("Integrator: " + integrator_str)
                                                        {
                                                            setup_and_solve_nlp(NN,
                                                                                NMF,
                                                                                con_str,
                                                                                cost_str,
                                                                                qp_solver_str,
                                                                                model_str,
                                                                                integrator_str);
                                                        }  // integrator
                                                    }
                                                }  // type of model
                                            }
                                        }  // qp solver
                                    }
                                }  // type of stage cost
                            }
                        }  // type of constraints
                    }

                }  // number of masses
            }

        }  // horizon lenght
    }
}  // TEST_CASE
