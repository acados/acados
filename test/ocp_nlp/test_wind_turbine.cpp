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
#include <vector>
#include <cstdlib>

#include "test/test_utils/eigen.h"
#include "catch/include/catch.hpp"

#include <assert.h>
#include <math.h>

#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"

#include "acados_c/external_function_interface.h"
#include "acados_c/ocp_nlp_interface.h"
#include "acados_c/ocp_qp_interface.h"

// TODO(dimitris): use only the strictly necessary includes here

#include "acados/utils/mem.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"

#include "acados/ocp_qp/ocp_qp_partial_condensing_solver.h"

#include "acados/ocp_nlp/ocp_nlp_sqp.h"
#include "acados/ocp_nlp/ocp_nlp_cost_common.h"
#include "acados/ocp_nlp/ocp_nlp_cost_ls.h"
#include "acados/ocp_nlp/ocp_nlp_cost_nls.h"
#include "acados/ocp_nlp/ocp_nlp_cost_external.h"
#include "acados/ocp_nlp/ocp_nlp_dynamics_cont.h"
#include "acados/ocp_nlp/ocp_nlp_constraints_bgh.h"

#include "acados/sim/sim_gnsf.h"

#include "examples/c/wt_model_nx6/nx6p2/wt_model.h"
#include "examples/c/wt_model_nx6/setup.c"

#define NN 40

#define MAX_SQP_ITERS 15
#define TOL 1e-6



ocp_qp_solver_t qp_solver_en(std::string const& inString)
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
    if (inString == "SPARSE_QORE") return PARTIAL_CONDENSING_OOQP;
#endif

    return (ocp_qp_solver_t) -1;
}



sim_solver_t integrator_en(std::string const& inString)
{
    if (inString == "ERK") return ERK;
    if (inString == "IRK") return IRK;
    if (inString == "LIFTED_IRK") return LIFTED_IRK;
    if (inString == "GNSF") return GNSF;

    return (sim_solver_t) -1;
}



static void shift_states(ocp_nlp_dims *dims, ocp_nlp_out *out, double *x_end)
{
    int N = dims->N;

    for (int i = 0; i < N; i++)
         blasfeo_dveccp(dims->nx[i], &out->ux[i], dims->nu[i], &out->ux[i+1], dims->nu[i+1]);
     blasfeo_pack_dvec(dims->nx[N], x_end, &out->ux[N], dims->nu[N]);
}



static void shift_controls(ocp_nlp_dims *dims, ocp_nlp_out *out, double *u_end)
{
    int N = dims->N;

    for (int i = 0; i < N-1; i++)
         blasfeo_dveccp(dims->nu[i], &out->ux[i], 0, &out->ux[i+1], 0);
     blasfeo_pack_dvec(dims->nu[N-1], u_end, &out->ux[N-1], 0);
}



static void select_dynamics_wt_casadi(int N,
    external_function_param_casadi *expl_vde_for,
    external_function_param_casadi *impl_ode_fun,
    external_function_param_casadi *impl_ode_fun_jac_x_xdot,
    external_function_param_casadi *impl_ode_jac_x_xdot_u,
    external_function_param_casadi *impl_ode_fun_jac_x_xdot_u,
    external_function_param_casadi *phi_fun,
    external_function_param_casadi *phi_fun_jac_y,
    external_function_param_casadi *phi_jac_y_uhat,
    external_function_param_casadi *f_lo_jac_x1_x1dot_u_z)
{
    for (int ii = 0; ii < N; ii++)
    {
        expl_vde_for[ii].casadi_fun = &wt_nx6p2_expl_vde_for;
        expl_vde_for[ii].casadi_work = &wt_nx6p2_expl_vde_for_work;
        expl_vde_for[ii].casadi_sparsity_in = &wt_nx6p2_expl_vde_for_sparsity_in;
        expl_vde_for[ii].casadi_sparsity_out = &wt_nx6p2_expl_vde_for_sparsity_out;
        expl_vde_for[ii].casadi_n_in = &wt_nx6p2_expl_vde_for_n_in;
        expl_vde_for[ii].casadi_n_out = &wt_nx6p2_expl_vde_for_n_out;

        impl_ode_fun[ii].casadi_fun = &wt_nx6p2_impl_ode_fun;
        impl_ode_fun[ii].casadi_work = &wt_nx6p2_impl_ode_fun_work;
        impl_ode_fun[ii].casadi_sparsity_in = &wt_nx6p2_impl_ode_fun_sparsity_in;
        impl_ode_fun[ii].casadi_sparsity_out = &wt_nx6p2_impl_ode_fun_sparsity_out;
        impl_ode_fun[ii].casadi_n_in = &wt_nx6p2_impl_ode_fun_n_in;
        impl_ode_fun[ii].casadi_n_out = &wt_nx6p2_impl_ode_fun_n_out;

        impl_ode_fun_jac_x_xdot[ii].casadi_fun = &wt_nx6p2_impl_ode_fun_jac_x_xdot;
        impl_ode_fun_jac_x_xdot[ii].casadi_work = &wt_nx6p2_impl_ode_fun_jac_x_xdot_work;
        impl_ode_fun_jac_x_xdot[ii].casadi_sparsity_in =
            &wt_nx6p2_impl_ode_fun_jac_x_xdot_sparsity_in;
        impl_ode_fun_jac_x_xdot[ii].casadi_sparsity_out =
            &wt_nx6p2_impl_ode_fun_jac_x_xdot_sparsity_out;
        impl_ode_fun_jac_x_xdot[ii].casadi_n_in = &wt_nx6p2_impl_ode_fun_jac_x_xdot_n_in;
        impl_ode_fun_jac_x_xdot[ii].casadi_n_out = &wt_nx6p2_impl_ode_fun_jac_x_xdot_n_out;

        impl_ode_jac_x_xdot_u[ii].casadi_fun = &wt_nx6p2_impl_ode_jac_x_xdot_u;
        impl_ode_jac_x_xdot_u[ii].casadi_work = &wt_nx6p2_impl_ode_jac_x_xdot_u_work;
        impl_ode_jac_x_xdot_u[ii].casadi_sparsity_in = &wt_nx6p2_impl_ode_jac_x_xdot_u_sparsity_in;
        impl_ode_jac_x_xdot_u[ii].casadi_sparsity_out =
            &wt_nx6p2_impl_ode_jac_x_xdot_u_sparsity_out;
        impl_ode_jac_x_xdot_u[ii].casadi_n_in = &wt_nx6p2_impl_ode_jac_x_xdot_u_n_in;
        impl_ode_jac_x_xdot_u[ii].casadi_n_out = &wt_nx6p2_impl_ode_jac_x_xdot_u_n_out;

        impl_ode_fun_jac_x_xdot_u[ii].casadi_fun = &wt_nx6p2_impl_ode_fun_jac_x_xdot_u;
        impl_ode_fun_jac_x_xdot_u[ii].casadi_work = &wt_nx6p2_impl_ode_fun_jac_x_xdot_u_work;
        impl_ode_fun_jac_x_xdot_u[ii].casadi_sparsity_in =
            &wt_nx6p2_impl_ode_fun_jac_x_xdot_u_sparsity_in;
        impl_ode_fun_jac_x_xdot_u[ii].casadi_sparsity_out =
            &wt_nx6p2_impl_ode_fun_jac_x_xdot_u_sparsity_out;
        impl_ode_fun_jac_x_xdot_u[ii].casadi_n_in = &wt_nx6p2_impl_ode_fun_jac_x_xdot_u_n_in;
        impl_ode_fun_jac_x_xdot_u[ii].casadi_n_out = &wt_nx6p2_impl_ode_fun_jac_x_xdot_u_n_out;

        // GNSF functions
        // phi_fun
        phi_fun[ii].casadi_fun            = &wt_nx6p2_phi_fun;
        phi_fun[ii].casadi_work           = &wt_nx6p2_phi_fun_work;
        phi_fun[ii].casadi_sparsity_in    = &wt_nx6p2_phi_fun_sparsity_in;
        phi_fun[ii].casadi_sparsity_out   = &wt_nx6p2_phi_fun_sparsity_out;
        phi_fun[ii].casadi_n_in           = &wt_nx6p2_phi_fun_n_in;
        phi_fun[ii].casadi_n_out          = &wt_nx6p2_phi_fun_n_out;

        phi_fun_jac_y[ii].casadi_fun = &wt_nx6p2_phi_fun_jac_y;
        phi_fun_jac_y[ii].casadi_work = &wt_nx6p2_phi_fun_jac_y_work;
        phi_fun_jac_y[ii].casadi_sparsity_in = &wt_nx6p2_phi_fun_jac_y_sparsity_in;
        phi_fun_jac_y[ii].casadi_sparsity_out = &wt_nx6p2_phi_fun_jac_y_sparsity_out;
        phi_fun_jac_y[ii].casadi_n_in = &wt_nx6p2_phi_fun_jac_y_n_in;
        phi_fun_jac_y[ii].casadi_n_out = &wt_nx6p2_phi_fun_jac_y_n_out;

        phi_jac_y_uhat[ii].casadi_fun = &wt_nx6p2_phi_jac_y_uhat;
        phi_jac_y_uhat[ii].casadi_work = &wt_nx6p2_phi_jac_y_uhat_work;
        phi_jac_y_uhat[ii].casadi_sparsity_in = &wt_nx6p2_phi_jac_y_uhat_sparsity_in;
        phi_jac_y_uhat[ii].casadi_sparsity_out = &wt_nx6p2_phi_jac_y_uhat_sparsity_out;
        phi_jac_y_uhat[ii].casadi_n_in = &wt_nx6p2_phi_jac_y_uhat_n_in;
        phi_jac_y_uhat[ii].casadi_n_out = &wt_nx6p2_phi_jac_y_uhat_n_out;

        // f_lo - linear output function
        f_lo_jac_x1_x1dot_u_z[ii].casadi_fun = &wt_nx6p2_f_lo_fun_jac_x1k1uz;
        f_lo_jac_x1_x1dot_u_z[ii].casadi_work = &wt_nx6p2_f_lo_fun_jac_x1k1uz_work;
        f_lo_jac_x1_x1dot_u_z[ii].casadi_sparsity_in = &wt_nx6p2_f_lo_fun_jac_x1k1uz_sparsity_in;
        f_lo_jac_x1_x1dot_u_z[ii].casadi_sparsity_out = &wt_nx6p2_f_lo_fun_jac_x1k1uz_sparsity_out;
        f_lo_jac_x1_x1dot_u_z[ii].casadi_n_in = &wt_nx6p2_f_lo_fun_jac_x1k1uz_n_in;
        f_lo_jac_x1_x1dot_u_z[ii].casadi_n_out = &wt_nx6p2_f_lo_fun_jac_x1k1uz_n_out;
    }
}



/************************************************
* nonlinear constraint
************************************************/

void ext_fun_h1(void *fun, ext_fun_arg_t *type_in, void **in, ext_fun_arg_t *type_out, void **out)
{
    int nu = 2;
    int nx = 8;
    int nh = 1;

    // scaling
    double alpha = 0.944*97/100;

    // ux
    struct blasfeo_dvec *ux = (struct blasfeo_dvec *) in[0];

    // h
    struct blasfeo_dvec_args *h_args = (struct blasfeo_dvec_args *) out[0];
    struct blasfeo_dvec *h = h_args->x;
    int xi = h_args->xi;
    BLASFEO_DVECEL(h, xi) = alpha * BLASFEO_DVECEL(ux, nu+0) * BLASFEO_DVECEL(ux, nu+5);

    // jac
    struct blasfeo_dmat_args *jac_args = (struct blasfeo_dmat_args *) out[1];
    struct blasfeo_dmat *jac = jac_args->A;
    int ai = jac_args->ai;
    int aj = jac_args->aj;
    blasfeo_dgese(nu+nx, nh, 0.0, jac, ai, aj);
    BLASFEO_DMATEL(jac, ai+nu+0, aj) = alpha * BLASFEO_DVECEL(ux, nu+5);
    BLASFEO_DMATEL(jac, ai+nu+5, aj) = alpha * BLASFEO_DVECEL(ux, nu+0);

    return;

}



void setup_and_solve_nlp(std::string const& integrator_str, std::string const& qp_solver_str)
{
    // _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
    int nx_ = 8;
    int nu_ = 2;
    int ny_ = 4;

    // number of local parametrs for each dynamics model function
    int np = 1;

    std::cout << integrator_str + " " + qp_solver_str << std::endl;

    /************************************************
    * problem dimensions
    ************************************************/

    int nx[NN+1] = {};
    int nu[NN+1] = {};
    int nbx[NN+1] = {};
    int nbu[NN+1] = {};
    int nb[NN+1] = {};
    int ng[NN+1] = {};
    int nh[NN+1] = {};
    int nq[NN+1] = {};
    int ns[NN+1] = {};
    int ny[NN+1] = {};
    int nz[NN+1] = {};

    // TODO(dimitris): setup bounds on states and controls based on ACADO controller
    nx[0] = nx_;
    nu[0] = nu_;
    nbx[0] = nx_;
    nbu[0] = nu_;
    nb[0] = nbu[0]+nbx[0];
    ng[0] = 0;

    // TODO(dimitris): add bilinear constraints later
    nh[0] = 0;
    ns[0] = 0;
    ny[0] = 4;  // ny_
    nz[0] = 0;

    for (int i = 1; i < NN; i++)
    {
        nx[i] = nx_;
        nu[i] = nu_;
        nbx[i] = 3;
        nbu[i] = nu_;
        nb[i] = nbu[i]+nbx[i];
        ng[i] = 0;
        nh[i] = 1;
        ns[i] = 1;
        ny[i] = 4;  // ny_
        nz[i] = 0;
    }

    nx[NN] = nx_;
    nu[NN] = 0;
    nbx[NN] = 3;
    nbu[NN] = 0;
    nb[NN] = nbu[NN]+nbx[NN];
    ng[NN] = 0;
    nh[NN] = 0;
    ns[NN] = 0;
    ny[NN] = 2;
    nz[NN] = 0;

    /************************************************
    * problem data
    ************************************************/

    double *x_end = (double *) malloc(sizeof(double)*nx_);
    double *u_end = (double *) malloc(sizeof(double)*nu_);

    // value of last stage when shifting states and controls
    for (int i = 0; i < nx_; i++) x_end[i] = 0.0;
    for (int i = 0; i < nu_; i++) u_end[i] = 0.0;



    /* constraints */

    // pitch angle rate
    double dbeta_min = - 8.0;
    double dbeta_max =   8.0;
    // generator torque
    double dM_gen_min = - 1.0;
    double dM_gen_max =   1.0;
    // generator angular velocity
    double OmegaR_min =  6.0/60*2*3.14159265359;
    double OmegaR_max = 13.0/60*2*3.14159265359;
    // pitch angle
    double beta_min =  0.0;
    double beta_max = 35.0;
    // generator torque
    double M_gen_min = 0.0;
    double M_gen_max = 5.0;
    // electric power
    double Pel_min = 0.0;
    double Pel_max = 5.0;


    /* soft constraints */

    // first stage
    int *idxs0 = (int *) malloc(ns[0]*sizeof(int));
    double *ls0 = (double *) malloc((ns[0])*sizeof(double));
    double *us0 = (double *) malloc((ns[0])*sizeof(double));

    // middle stage
    int *idxs1 = (int *) malloc(ns[1]*sizeof(int));
    double *ls1 = (double *) malloc((ns[1])*sizeof(double));
    double *us1 = (double *) malloc((ns[1])*sizeof(double));

    // last stage
    int *idxsN = (int *) malloc(ns[NN]*sizeof(int));
    double *lsN = (double *) malloc((ns[NN])*sizeof(double));
    double *usN = (double *) malloc((ns[NN])*sizeof(double));


    /* box constraints */

    // acados inf
    double acados_inf = 1e8;

    // first stage
    int *idxb0 = (int *) malloc(nb[0]*sizeof(int));
    double *lb0 = (double *) malloc((nb[0])*sizeof(double));
    double *ub0 = (double *) malloc((nb[0])*sizeof(double));

    // pitch angle rate
    idxb0[0] = 0;
    lb0[0] = dbeta_min;
    ub0[0] = dbeta_max;

    // generator torque
    idxb0[1] = 1;
    lb0[1] = dM_gen_min;
    ub0[1] = dM_gen_max;

    // dummy state bounds
    for (int ii = 0; ii < nbx[0]; ii++)
    {
        idxb0[nbu[0]+ii] = nbu[0]+ii;
        lb0[nbu[0]+ii] = - acados_inf;
        ub0[nbu[0]+ii] =   acados_inf;
    }


    // middle stages
    int *idxb1 = (int *) malloc(nb[1]*sizeof(int));
    double *lb1 = (double *) malloc((nb[1])*sizeof(double));
    double *ub1 = (double *) malloc((nb[1])*sizeof(double));

    // pitch angle rate
    idxb1[0] = 0;
    lb1[0] = dbeta_min;
    ub1[0] = dbeta_max;

    // generator torque rate
    idxb1[1] = 1;
    lb1[1] = dM_gen_min;
    ub1[1] = dM_gen_max;

    // generator angular velocity
    idxb1[2] = 2;
    lb1[2] = OmegaR_min;
    ub1[2] = OmegaR_max;

    // pitch angle
    idxb1[3] = 8;
    lb1[3] = beta_min;
    ub1[3] = beta_max;

    // generator torque
    idxb1[4] = 9;
    lb1[4] = M_gen_min;
    ub1[4] = M_gen_max;

    // last stage
    int *idxbN = (int *) malloc(nb[NN]*sizeof(int));
    double *lbN = (double *) malloc((nb[NN])*sizeof(double));
    double *ubN = (double *) malloc((nb[NN])*sizeof(double));

    // generator angular velocity
    idxbN[0] = 0;
    lbN[0] = OmegaR_min;
    ubN[0] = OmegaR_max;

    // pitch angle
    idxbN[1] = 6;
    lbN[1] = beta_min;
    ubN[1] = beta_max;

    // generator torque
    idxbN[2] = 7;
    lbN[2] = M_gen_min;
    ubN[2] = M_gen_max;

    /* nonlinear constraints */

    // middle stages
    external_function_generic h1;
    double *lh1;
    double *uh1;
    lh1 = (double *) malloc((nh[1])*sizeof(double));
    uh1 = (double *) malloc((nh[1])*sizeof(double));

    if (nh[1] > 0)
    {
        h1.evaluate = &ext_fun_h1;

        // electric power
        lh1[0] = Pel_min;
        uh1[0] = Pel_max;
    }

    // softed
    if (ns[1] > 0)
    {
        idxs1[0] = nb[1]+ng[1];
        ls1[0] = 0.0;
        us1[0] = 0.0;
    }



    /* linear least squares */

    // output definition
    // y = {x[0], x[4]; u[0]; u[1]; u[2]};
    //   = Vx * x + Vu * u

    double *Vx = (double *) malloc((ny_*nx_)*sizeof(double));
    for (int ii = 0; ii < ny_*nx_; ii++)
        Vx[ii] = 0.0;
    Vx[0+ny_*0] = 1.0;
    Vx[1+ny_*4] = 1.0;

    double *Vu = (double *) malloc((ny_*nu_)*sizeof(double));
    for (int ii = 0; ii < ny_*nu_; ii++)
        Vu[ii] = 0.0;
    Vu[2+ny_*0] = 1.0;
    Vu[3+ny_*1] = 1.0;

    double *W = (double *) malloc((ny_*ny_)*sizeof(double));
    for (int ii = 0; ii < ny_*ny_; ii++)
        W[ii] = 0.0;
    W[0+ny_*0] = 1.5114;
    W[1+ny_*0] = -0.0649;
    W[0+ny_*1] = -0.0649;
    W[1+ny_*1] = 0.0180;
    W[2+ny_*2] = 0.01;
    W[3+ny_*3] = 0.001;

    /* slacks */

    // first stage
    double *lZ0 = (double *) malloc(ns[0]*sizeof(double));
    double *uZ0 = (double *) malloc(ns[0]*sizeof(double));
    double *lz0 = (double *) malloc(ns[0]*sizeof(double));
    double *uz0 = (double *) malloc(ns[0]*sizeof(double));

    // middle stages
    double *lZ1 = (double *) malloc(ns[1]*sizeof(double));
    double *uZ1 = (double *) malloc(ns[1]*sizeof(double));
    double *lz1 = (double *) malloc(ns[1]*sizeof(double));
    double *uz1 = (double *) malloc(ns[1]*sizeof(double));
    lZ1[0] = 1e2;
    uZ1[0] = 1e2;
    lz1[0] = 0e1;
    uz1[0] = 0e1;

    // final stage
    double *lZN = (double *) malloc(ns[NN]*sizeof(double));
    double *uZN = (double *) malloc(ns[NN]*sizeof(double));
    double *lzN = (double *) malloc(ns[NN]*sizeof(double));
    double *uzN = (double *) malloc(ns[NN]*sizeof(double));

    /************************************************
    * plan + config
    ************************************************/

    ocp_nlp_solver_plan *plan = ocp_nlp_plan_create(NN);

    plan->nlp_solver = SQP;

    for (int i = 0; i <= NN; i++)
        plan->nlp_cost[i] = LINEAR_LS;

    ocp_qp_solver_t qp_solver_type = qp_solver_en(qp_solver_str);
    plan->ocp_qp_solver_plan.qp_solver = qp_solver_type;

    sim_solver_t integrator_type = integrator_en(integrator_str);

    switch (integrator_type)
    {
        case IRK:
            for (int i = 0; i < NN; i++)
            {
                plan->nlp_dynamics[i] = CONTINUOUS_MODEL;
                plan->sim_solver_plan[i].sim_solver = IRK;
            }
            break;

        case ERK:
            for (int i = 0; i < NN; i++)
            {
                plan->nlp_dynamics[i] = CONTINUOUS_MODEL;
                plan->sim_solver_plan[i].sim_solver = ERK;
            }
            break;

        case LIFTED_IRK:
            for (int i = 0; i < NN; i++)
            {
                plan->nlp_dynamics[i] = CONTINUOUS_MODEL;
                plan->sim_solver_plan[i].sim_solver = LIFTED_IRK;
            }
            break;

        case GNSF:
            for (int i = 0; i < NN; i++)
            {
                plan->nlp_dynamics[i] = CONTINUOUS_MODEL;
                plan->sim_solver_plan[i].sim_solver = GNSF;
            }
            break;

        default:
            for (int i = 0; i < NN; i++)
            {
                plan->nlp_dynamics[i] = CONTINUOUS_MODEL;

                if (i%4 == 0)
                    plan->sim_solver_plan[i].sim_solver = IRK;
                else if (i%4 == 1)
                    plan->sim_solver_plan[i].sim_solver = ERK;
                else if (i%4 == 2)
                    plan->sim_solver_plan[i].sim_solver = LIFTED_IRK;
                else if (i%4 == 3)
                    plan->sim_solver_plan[i].sim_solver = GNSF;

            }
            break;
    }

    for (int i = 0; i <= NN; i++)
    {
        plan->nlp_constraints[i] = BGH;
    }

    ocp_nlp_solver_config *config = ocp_nlp_config_create(*plan, NN);

    /************************************************
    * ocp_nlp_dims
    ************************************************/

    ocp_nlp_dims *dims = ocp_nlp_dims_create(config);
    ocp_nlp_dims_initialize(config, nx, nu, ny, nbx, nbu, ng, nh, nq, ns, nz, dims);

    /************************************************
    * dynamics
    ************************************************/

    // explicit model
    external_function_param_casadi *expl_vde_for =
    (external_function_param_casadi *) malloc(NN*sizeof(external_function_param_casadi));
    // implicit model
    external_function_param_casadi *impl_ode_fun =
        (external_function_param_casadi *) malloc(NN*sizeof(external_function_param_casadi));
    external_function_param_casadi *impl_ode_fun_jac_x_xdot =
        (external_function_param_casadi *) malloc(NN*sizeof(external_function_param_casadi));
    external_function_param_casadi *impl_ode_jac_x_xdot_u =
        (external_function_param_casadi *) malloc(NN*sizeof(external_function_param_casadi));
    external_function_param_casadi *impl_ode_fun_jac_x_xdot_u =
        (external_function_param_casadi *) malloc(NN*sizeof(external_function_param_casadi));
    // gnsf model
    external_function_param_casadi *phi_fun =
        (external_function_param_casadi *) malloc(NN*sizeof(external_function_param_casadi));
    external_function_param_casadi *phi_fun_jac_y =
        (external_function_param_casadi *) malloc(NN*sizeof(external_function_param_casadi));
    external_function_param_casadi *phi_jac_y_uhat =
        (external_function_param_casadi *) malloc(NN*sizeof(external_function_param_casadi));
    external_function_param_casadi *f_lo_jac_x1_x1dot_u_z =
        (external_function_param_casadi *) malloc(NN*sizeof(external_function_param_casadi));

    select_dynamics_wt_casadi(NN, expl_vde_for, impl_ode_fun, impl_ode_fun_jac_x_xdot,
        impl_ode_jac_x_xdot_u, impl_ode_fun_jac_x_xdot_u, phi_fun, phi_fun_jac_y, phi_jac_y_uhat,
        f_lo_jac_x1_x1dot_u_z);

    // explicit model
    external_function_param_casadi_create_array(NN, expl_vde_for, np);
    // implicit model
    external_function_param_casadi_create_array(NN, impl_ode_fun, np);
    external_function_param_casadi_create_array(NN, impl_ode_fun_jac_x_xdot, np);
    external_function_param_casadi_create_array(NN, impl_ode_jac_x_xdot_u, np);
    external_function_param_casadi_create_array(NN, impl_ode_fun_jac_x_xdot_u, np);
    // gnsf model
    external_function_param_casadi_create_array(NN, phi_fun, np);
    external_function_param_casadi_create_array(NN, phi_fun_jac_y, np);
    external_function_param_casadi_create_array(NN, phi_jac_y_uhat, np);
    external_function_param_casadi_create_array(NN, f_lo_jac_x1_x1dot_u_z, np);

    // GNSF import matrices function
    external_function_casadi get_matrices_fun;
    get_matrices_fun.casadi_fun            = &wt_nx6p2_get_matrices_fun;
    get_matrices_fun.casadi_work           = &wt_nx6p2_get_matrices_fun_work;
    get_matrices_fun.casadi_sparsity_in    = &wt_nx6p2_get_matrices_fun_sparsity_in;
    get_matrices_fun.casadi_sparsity_out   = &wt_nx6p2_get_matrices_fun_sparsity_out;
    get_matrices_fun.casadi_n_in           = &wt_nx6p2_get_matrices_fun_n_in;
    get_matrices_fun.casadi_n_out          = &wt_nx6p2_get_matrices_fun_n_out;
    external_function_casadi_create(&get_matrices_fun);

    external_function_generic *get_model_matrices = (external_function_generic *) &get_matrices_fun;

    for (int i = 0; i < NN; i++)
    {
        if (plan->sim_solver_plan[i].sim_solver == GNSF)
        {
            /* initialize additional gnsf dimensions */
            ocp_nlp_dynamics_cont_dims *dyn_dims = (ocp_nlp_dynamics_cont_dims *) dims->dynamics[i];
            sim_gnsf_dims *gnsf_dims = (sim_gnsf_dims *) dyn_dims->sim;

            gnsf_dims->nx1 = 8;
            gnsf_dims->nz = 0;
            gnsf_dims->nz1 = 0;
            gnsf_dims->n_out = 1;
            gnsf_dims->ny = 5;
            gnsf_dims->nuhat = 0;
        }
    }


    /************************************************
    * nlp_in
    ************************************************/

    ocp_nlp_in *nlp_in = ocp_nlp_in_create(config, dims);

    // sampling times
    for (int ii = 0; ii < NN; ii++)
    {
        nlp_in->Ts[ii] = 0.2;
    }

    // output definition: y = [x; u]

    /* cost */

    // linear ls
    ocp_nlp_cost_ls_model **cost = (ocp_nlp_cost_ls_model **) nlp_in->cost;

    for (int i = 0; i <= NN; i++)
    {
        // Cyt
        blasfeo_pack_tran_dmat(ny[i], nu[i], Vu, ny_, &cost[i]->Cyt, 0, 0);
        blasfeo_pack_tran_dmat(ny[i], nx[i], Vx, ny_, &cost[i]->Cyt, nu[i], 0);

        // W
        blasfeo_pack_dmat(ny[i], ny[i], W, ny_, &cost[i]->W, 0, 0);
    }

    // slacks (middle stages)
    for (int ii = 1; ii < NN; ii++)
    {
        blasfeo_pack_dvec(ns[ii], lZ1, &cost[ii]->Z, 0);
        blasfeo_pack_dvec(ns[ii], uZ1, &cost[ii]->Z, ns[ii]);
        blasfeo_pack_dvec(ns[ii], lz1, &cost[ii]->z, 0);
        blasfeo_pack_dvec(ns[ii], uz1, &cost[ii]->z, ns[ii]);
    }


    /* dynamics */

    int set_fun_status;

    for (int i = 0; i < NN; i++)
    {
        if (plan->sim_solver_plan[i].sim_solver == ERK)
        {
            set_fun_status = nlp_set_model_in_stage(config, nlp_in, i, "expl_vde_for",
                &expl_vde_for[i]);
            REQUIRE(set_fun_status == 0);
            if (set_fun_status != 0) exit(1);
        }
        else if (plan->sim_solver_plan[i].sim_solver == IRK)
        {
            set_fun_status = nlp_set_model_in_stage(config, nlp_in, i, "impl_ode_fun",
                &impl_ode_fun[i]);
            REQUIRE(set_fun_status == 0);
            if (set_fun_status != 0) exit(1);
            set_fun_status = nlp_set_model_in_stage(config, nlp_in, i, "impl_ode_fun_jac_x_xdot",
                &impl_ode_fun_jac_x_xdot[i]);
            REQUIRE(set_fun_status == 0);
            if (set_fun_status != 0) exit(1);
            set_fun_status = nlp_set_model_in_stage(config, nlp_in, i, "impl_ode_jac_x_xdot_u",
                &impl_ode_jac_x_xdot_u[i]);
            REQUIRE(set_fun_status == 0);
            if (set_fun_status != 0) exit(1);
        }
        else if (plan->sim_solver_plan[i].sim_solver == GNSF)
        {
            set_fun_status = nlp_set_model_in_stage(config, nlp_in, i, "phi_fun", &phi_fun[i]);
            REQUIRE(set_fun_status == 0);
            if (set_fun_status != 0) exit(1);
            set_fun_status = nlp_set_model_in_stage(config, nlp_in, i, "phi_fun_jac_y",
                &phi_fun_jac_y[i]);
            REQUIRE(set_fun_status == 0);
            if (set_fun_status != 0) exit(1);
            set_fun_status = nlp_set_model_in_stage(config, nlp_in, i, "phi_jac_y_uhat",
                &phi_jac_y_uhat[i]);
            REQUIRE(set_fun_status == 0);
            if (set_fun_status != 0) exit(1);
            set_fun_status = nlp_set_model_in_stage(config, nlp_in, i, "f_lo_jac_x1_x1dot_u_z",
                &f_lo_jac_x1_x1dot_u_z[i]);
        }
        else if (plan->sim_solver_plan[i].sim_solver == LIFTED_IRK)
        {
            set_fun_status = nlp_set_model_in_stage(config, nlp_in, i, "impl_ode_fun",
                &impl_ode_fun[i]);
            REQUIRE(set_fun_status == 0);
            if (set_fun_status != 0) exit(1);
            set_fun_status = nlp_set_model_in_stage(config, nlp_in, i, "impl_ode_fun_jac_x_xdot_u",
                &impl_ode_fun_jac_x_xdot_u[i]);
            REQUIRE(set_fun_status == 0);
            if (set_fun_status != 0) exit(1);
        }
        else
        {
            printf("\nWrong sim name\n\n");
            exit(1);
        }
    }


    /* constraints */

    ocp_nlp_constraints_bgh_model **constraints =
        (ocp_nlp_constraints_bgh_model **) nlp_in->constraints;

    /* box constraints */

    // fist stage
    blasfeo_pack_dvec(nb[0], lb0, &constraints[0]->d, 0);
    blasfeo_pack_dvec(nb[0], ub0, &constraints[0]->d, nb[0]+ng[0]+nh[0]);
    for (int ii = 0; ii < nb[0]; ii++) constraints[0]->idxb[ii] = idxb0[ii];
    // middle stages
    for (int i = 1; i < NN; i++)
    {
        blasfeo_pack_dvec(nb[i], lb1, &constraints[i]->d, 0);
        blasfeo_pack_dvec(nb[i], ub1, &constraints[i]->d, nb[i]+ng[i]+nh[i]);
        for (int ii = 0; ii < nb[i]; ii++) constraints[i]->idxb[ii] = idxb1[ii];
    }
    // last stage
    blasfeo_pack_dvec(nb[NN], lbN, &constraints[NN]->d, 0);
    blasfeo_pack_dvec(nb[NN], ubN, &constraints[NN]->d, nb[NN]+ng[NN]+nh[NN]);
    for (int ii = 0; ii < nb[NN]; ii++) constraints[NN]->idxb[ii] = idxbN[ii];

    /* nonlinear constraints */

    // middle stages
    for (int i = 1; i < NN; i++)
    {
        if (nh[i] > 0)
        {
            blasfeo_pack_dvec(nh[i], lh1, &constraints[i]->d, nb[i]+ng[i]);
            blasfeo_pack_dvec(nh[i], uh1, &constraints[i]->d, 2*nb[i]+2*ng[i]+nh[i]);
            constraints[i]->h = &h1;
        }
    }

    /* soft constraints */

    // middle stages
    for (int i = 1; i < NN; i++)
    {
        if (ns[i] > 0)
        {
            blasfeo_pack_dvec(ns[i], ls1, &constraints[i]->d, 2*nb[i]+2*ng[i]+2*nh[i]);
            blasfeo_pack_dvec(ns[i], us1, &constraints[i]->d, 2*nb[i]+2*ng[i]+2*nh[i]+ns[i]);
            for (int ii = 0; ii < ns[i]; ii++) constraints[i]->idxs[ii] = idxs1[ii];
        }
    }

    /************************************************
    * sqp opts
    ************************************************/

    void *nlp_opts = ocp_nlp_opts_create(config, dims);
    ocp_nlp_sqp_opts *sqp_opts = (ocp_nlp_sqp_opts *) nlp_opts;

    for (int i = 0; i < NN; ++i)
    {
        ocp_nlp_dynamics_cont_opts *dynamics_stage_opts =
            (ocp_nlp_dynamics_cont_opts *) sqp_opts->dynamics[i];
        sim_rk_opts *sim_opts = (sim_rk_opts *) dynamics_stage_opts->sim_solver;

        if (plan->sim_solver_plan[i].sim_solver == ERK)
        {
            sim_opts->ns = 4;
            sim_opts->num_steps = 10;
        }
        else if (plan->sim_solver_plan[i].sim_solver == IRK)
        {
            sim_opts->ns = 4;
            sim_opts->num_steps = 1;
            sim_opts->jac_reuse = true;
        }
        else if (plan->sim_solver_plan[i].sim_solver == LIFTED_IRK)
        {
            sim_opts->ns = 4;
            sim_opts->num_steps = 1;
        }
        else if (plan->sim_solver_plan[i].sim_solver == GNSF)
        {
            sim_opts->ns = 4;
            sim_opts->num_steps = 1;
            sim_opts->newton_iter = 1;
            sim_opts->jac_reuse = true;
        }
    }

    sqp_opts->maxIter = MAX_SQP_ITERS;
    sqp_opts->min_res_g = 1e-6;
    sqp_opts->min_res_b = 1e-8;
    sqp_opts->min_res_d = 1e-8;
    sqp_opts->min_res_m = 1e-8;

    // partial condensing
    if (plan->ocp_qp_solver_plan.qp_solver == PARTIAL_CONDENSING_HPIPM)
    {
        ocp_nlp_sqp_opts *sqp_opts = (ocp_nlp_sqp_opts *) nlp_opts;
        ocp_qp_partial_condensing_solver_opts *pcond_solver_opts =
            (ocp_qp_partial_condensing_solver_opts *) sqp_opts->qp_solver_opts;
        pcond_solver_opts->pcond_opts->N2 = 10;
    }

    config->opts_update(config, dims, nlp_opts);

    /************************************************
    * ocp_nlp out
    ************************************************/

    ocp_nlp_out *nlp_out = ocp_nlp_out_create(config, dims);

    ocp_nlp_solver *solver = ocp_nlp_create(config, dims, nlp_opts);

    /************************************************
    * 	precomputation (after all options are set)
    ************************************************/

    for (int i = 0; i < NN; i++){
        if (plan->sim_solver_plan[i].sim_solver == GNSF)
        {
            ocp_nlp_dynamics_cont_model *dynamics =
                (ocp_nlp_dynamics_cont_model *) nlp_in->dynamics[i];
            gnsf_model* model = (gnsf_model *)dynamics->sim_model;

            // get gnsf_dims
            ocp_nlp_dynamics_cont_dims *dyn_dims = (ocp_nlp_dynamics_cont_dims *) dims->dynamics[i];
            sim_gnsf_dims *gnsf_dims = (sim_gnsf_dims *) dyn_dims->sim;

            // get sim opts
            ocp_nlp_dynamics_cont_opts *dynamics_stage_opts =
                (ocp_nlp_dynamics_cont_opts *) sqp_opts->dynamics[i];
            sim_rk_opts *sim_opts = (sim_rk_opts *) dynamics_stage_opts->sim_solver;

            // import model matrices
            sim_gnsf_import_matrices(gnsf_dims, model, get_model_matrices);

            // get sim_solver_config
            sim_solver_config *sim_sol_config =
                (sim_solver_config *) config->dynamics[i]->sim_solver;

            // get sim_solver memory
            ocp_nlp_sqp_memory *mem = (ocp_nlp_sqp_memory *) solver->mem;
            ocp_nlp_dynamics_cont_memory* dynamics_mem =
                (ocp_nlp_dynamics_cont_memory *) mem->dynamics[i];
            char *mem_ptr = (char *) dynamics_mem;

            // mem_ptr now points to the memory of the integrator;
            mem_ptr += sizeof(ocp_nlp_dynamics_cont_memory);

            // precompute
            sim_gnsf_precompute(sim_sol_config, gnsf_dims, model, sim_opts, mem_ptr, solver->work,
                nlp_in->Ts[i]);
            // NOTE; solver->work can be used, as it is for sure larger than the workspace
            //       needed to precompute, as the latter is part of the first.
        }
    }


    /************************************************
    * sqp solve
    ************************************************/

    int nmpc_problems = 40;

    int status;

    acados_timer timer;
    acados_tic(&timer);

    // warm start output initial guess of solution
    for (int i = 0; i <= NN; i++)
    {
        blasfeo_pack_dvec(2, u0_ref, nlp_out->ux+i, 0);
        blasfeo_pack_dvec(1, wind0_ref+i, nlp_out->ux+i, 2);
        blasfeo_pack_dvec(nx[i], x0_ref, nlp_out->ux+i, nu[i]);
    }

    // update x0 as box constraint
    blasfeo_pack_dvec(nx[0], x0_ref, &constraints[0]->d, nbu[0]);
    blasfeo_pack_dvec(nx[0], x0_ref, &constraints[0]->d, nb[0]+ng[0]+nh[0]+nbu[0]);

    for (int idx = 0; idx < nmpc_problems; idx++)
    {
        // update wind distrurbance as external function parameter
        for (int ii = 0; ii < NN; ii++)
        {
            if (plan->sim_solver_plan[ii].sim_solver == ERK)
            {
                expl_vde_for[ii].set_param(expl_vde_for+ii, wind0_ref+idx+ii);
            }
            else if ((plan->sim_solver_plan[ii].sim_solver == IRK) |
                     (plan->sim_solver_plan[ii].sim_solver == LIFTED_IRK))
            {
                impl_ode_fun[ii].set_param(impl_ode_fun+ii, wind0_ref+idx+ii);
                impl_ode_fun_jac_x_xdot[ii].set_param(impl_ode_fun_jac_x_xdot+ii, wind0_ref+idx+ii);
                impl_ode_jac_x_xdot_u[ii].set_param(impl_ode_jac_x_xdot_u+ii, wind0_ref+idx+ii);
                impl_ode_fun_jac_x_xdot_u[ii].set_param(impl_ode_fun_jac_x_xdot_u+ii,
                    wind0_ref+idx+ii);
            }
            else if (plan->sim_solver_plan[ii].sim_solver == GNSF)
            {
                phi_fun[ii].set_param(phi_fun+ii, wind0_ref+idx+ii);
                phi_fun_jac_y[ii].set_param(phi_fun_jac_y+ii, wind0_ref+idx+ii);
                phi_jac_y_uhat[ii].set_param(phi_jac_y_uhat+ii, wind0_ref+idx+ii);
                f_lo_jac_x1_x1dot_u_z[ii].set_param(f_lo_jac_x1_x1dot_u_z+ii, wind0_ref+idx+ii);
            }
            else
            {
                printf("\nWrong sim name\n\n");
                exit(1);
            }
        }
        // update reference
        for (int i = 0; i <= NN; i++)
        {
            BLASFEO_DVECEL(&cost[i]->y_ref, 0) = y_ref[(idx + i)*4+0];
            BLASFEO_DVECEL(&cost[i]->y_ref, 1) = y_ref[(idx + i)*4+1];
            if (i < NN)
            {
                BLASFEO_DVECEL(&cost[i]->y_ref, 2) = y_ref[(idx + i)*4+2];
                BLASFEO_DVECEL(&cost[i]->y_ref, 3) = y_ref[(idx + i)*4+3];
            }
        }

        // solve NLP
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

        // update initial condition
        // TODO(dimitris): maybe simulate system instead of passing x[1] as next state
        blasfeo_dveccp(nx_, &nlp_out->ux[1], nu_, &constraints[0]->d, nbu[0]);
        blasfeo_dveccp(nx_, &nlp_out->ux[1], nu_, &constraints[0]->d, nbu[0]+nb[0]+ng[0]);

        // print info
        printf("\nproblem #%d, status %d, iters %d\n", idx, status,
            ((ocp_nlp_sqp_memory *)solver->mem)->sqp_iter);
        printf("xsim = \n");
        blasfeo_print_tran_dvec(dims->nx[0], &nlp_out->ux[0], dims->nu[0]);
        printf("electrical power = %f\n",
            0.944*97/100*BLASFEO_DVECEL(&nlp_out->ux[0], 2)*BLASFEO_DVECEL(&nlp_out->ux[0], 7));
        printf("Max residuals = %e\n", max_res);

        REQUIRE((status == 0 || status == 1 && MAX_SQP_ITERS == 1));
        REQUIRE(max_res <= TOL);

        // shift trajectories
        blasfeo_unpack_dvec(dims->nx[NN], &nlp_out->ux[NN-1], dims->nu[NN-1], x_end);
        blasfeo_unpack_dvec(dims->nu[NN-1], &nlp_out->ux[NN-2], dims->nu[NN-2], u_end);

        shift_states(dims, nlp_out, x_end);
        shift_controls(dims, nlp_out, u_end);
    }

    double time = acados_toc(&timer);

    printf("\n\ntotal time (including printing) = %f ms (time per SQP = %f)\n\n",
        time*1e3, time*1e3/nmpc_problems);


    /************************************************
    * free memory
    ************************************************/

    // TODO(dimitris): VALGRIND!
    external_function_casadi_free(&get_matrices_fun);

    external_function_param_casadi_free(expl_vde_for);
    external_function_param_casadi_free(impl_ode_fun);
    external_function_param_casadi_free(impl_ode_fun_jac_x_xdot);
    external_function_param_casadi_free(impl_ode_jac_x_xdot_u);
    external_function_param_casadi_free(impl_ode_fun_jac_x_xdot_u);
    external_function_param_casadi_free(phi_fun);
    external_function_param_casadi_free(phi_fun_jac_y);
    external_function_param_casadi_free(phi_jac_y_uhat);
    external_function_param_casadi_free(f_lo_jac_x1_x1dot_u_z);

    free(expl_vde_for);
    free(impl_ode_fun);
    free(impl_ode_fun_jac_x_xdot);
    free(impl_ode_jac_x_xdot_u);
    free(impl_ode_fun_jac_x_xdot_u);

    free(phi_fun);
    free(phi_fun_jac_y);
    free(phi_jac_y_uhat);
    free(f_lo_jac_x1_x1dot_u_z);


    free(nlp_opts);
    free(nlp_in);
    free(nlp_out);
    free(solver);
    free(dims);
    free(config);
    free(plan);

    free(lb0);
    free(ub0);
    free(lb1);
    free(ub1);
    free(lbN);
    free(ubN);

    free(lh1);
    free(uh1);
    free(Vx);
    free(Vu);
    free(W);

    free(idxb0);
    free(idxb1);
    free(idxbN);

    free(ls0);
    free(us0);
    free(ls1);
    free(us1);
    free(lsN);
    free(usN);
    free(idxs0);
    free(idxs1);
    free(idxsN);

    free(lZ0);
    free(uZ0);
    free(lz0);
    free(uz0);
    free(lZ1);
    free(uZ1);
    free(lz1);
    free(uz1);
    free(lZN);
    free(uZN);
    free(lzN);
    free(uzN);

    free(x_end);
    free(u_end);
}



/************************************************
* TEST CASE: wind turbine
************************************************/

TEST_CASE("wind turbine nmpc", "[NLP solver]")
{
    std::vector<std::string> integrators = {"IRK", "ERK", "LIFTED_IRK", "GNSF", "MIXED"};
    std::vector<std::string> qp_solvers = { "SPARSE_HPIPM",
                                            // "SPARSE_HPMPC",
                                            // "SPARSE_QPDUNES",
                                            "DENSE_HPIPM",
                                            // "DENSE_QPOASES"
#ifdef ACADOS_WITH_OOQP
                                            // "SPARSE_OOQP",
                                            // "DENSE_OOQP"
#endif
#ifdef ACADOS_WITH_QORE
                                            // ,"DENSE_QORE"
                                            };
#else
    };
#endif

    for (std::string qp_solver_str : qp_solvers)
    {
        SECTION("QP solver: " + qp_solver_str)
        {
            for (std::string integrator_str : integrators)
            {
                SECTION("Integrator: " + integrator_str)
                {
                    setup_and_solve_nlp(integrator_str, qp_solver_str);
                }
            }
        }
    }
}
