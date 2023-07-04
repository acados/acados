/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
 */

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>

#include "test/test_utils/eigen.h"
#include "catch/include/catch.hpp"

// std
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// blasfeo
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"

// acados
#include "acados_c/external_function_interface.h"
#include "acados_c/ocp_nlp_interface.h"
#include "acados/ocp_nlp/ocp_nlp_sqp.h"

#include "acados/utils/mem.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"

// example specific
#include "examples/c/wt_model_nx6/nx6p2/wt_model.h"
#include "examples/c/wt_model_nx6/setup.c"

#define NN 40

#define MAX_SQP_ITERS 15
#define TOL 1e-6



ocp_qp_solver_t qp_solver_en(std::string const& inString)
{
    if (inString == "SPARSE_HPIPM") return PARTIAL_CONDENSING_HPIPM;
    if (inString == "DENSE_HPIPM") return FULL_CONDENSING_HPIPM;
#ifdef ACADOS_WITH_HPMPC
    if (inString == "SPARSE_HPMPC") return PARTIAL_CONDENSING_HPMPC;
#endif
#ifdef ACADOS_WITH_QPOASES
    if (inString == "DENSE_QPOASES") return FULL_CONDENSING_QPOASES;
#endif
#ifdef ACADOS_WITH_QPDUNES
    if (inString == "SPARSE_QPDUNES") return PARTIAL_CONDENSING_QPDUNES;
#endif
#ifdef ACADOS_WITH_OOQP
    if (inString == "DENSE_OOQP") return FULL_CONDENSING_OOQP;
    if (inString == "SPARSE_OOQP") return PARTIAL_CONDENSING_OOQP;
#endif
#ifdef ACADOS_WITH_OSQP
    if (inString == "SPARSE_OSQP") return PARTIAL_CONDENSING_OSQP;
#endif
#ifdef ACADOS_WITH_QORE
    if (inString == "DENSE_QORE") return FULL_CONDENSING_QORE;
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
     blasfeo_pack_dvec(dims->nx[N], x_end, 1, &out->ux[N], dims->nu[N]);
}



static void shift_controls(ocp_nlp_dims *dims, ocp_nlp_out *out, double *u_end)
{
    int N = dims->N;

    for (int i = 0; i < N-1; i++)
         blasfeo_dveccp(dims->nu[i], &out->ux[i], 0, &out->ux[i+1], 0);
     blasfeo_pack_dvec(dims->nu[N-1], u_end, 1, &out->ux[N-1], 0);
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
    // struct blasfeo_dvec *ux = in[0];
    struct blasfeo_dvec_args *x_args = (struct blasfeo_dvec_args *) in[0];
    // struct blasfeo_dvec_args *u_args = in[1];

    struct blasfeo_dvec *x = x_args->x;
    // struct blasfeo_dvec *u = u_args->x;

    int x_offset = x_args->xi;
    // int u_offset = u_args->xi;

    // h
    struct blasfeo_dvec_args *h_args = (struct blasfeo_dvec_args *) out[0];
    struct blasfeo_dvec *h = h_args->x;
    int xi = h_args->xi;
    BLASFEO_DVECEL(h, xi) = alpha * BLASFEO_DVECEL(x, x_offset) * BLASFEO_DVECEL(x, x_offset+5);

    // jac
    struct blasfeo_dmat_args *jac_args = (struct blasfeo_dmat_args *) out[1];
    struct blasfeo_dmat *jac = jac_args->A;
    int ai = jac_args->ai;
    int aj = jac_args->aj;
    blasfeo_dgese(nu+nx, nh, 0.0, jac, ai, aj);
    BLASFEO_DMATEL(jac, ai+nu+0, aj) = alpha * BLASFEO_DVECEL(x, x_offset+5);
    BLASFEO_DMATEL(jac, ai+nu+5, aj) = alpha * BLASFEO_DVECEL(x, x_offset+0);

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

    // optimization variables
    int nx[NN+1] = {}; // states
    int nu[NN+1] = {}; // inputs
    int nz[NN+1] = {}; // algebraic variables
    int ns[NN+1] = {}; // slacks
    // cost
    int ny[NN+1] = {}; // measurements
    // constraints
    int nbx[NN+1] = {}; // state bounds
    int nbu[NN+1] = {}; // input bounds
    int ng[NN+1] = {}; // general linear constraints
    int nh[NN+1] = {}; // nonlinear constraints
    int nsh[NN+1] = {}; // softed nonlinear constraints

    // TODO(dimitris): setup bounds on states and controls based on ACADO controller
    nx[0] = nx_;
    nu[0] = nu_;
    nbx[0] = nx_;
    nbu[0] = nu_;
    ng[0] = 0;
    // TODO(dimitris): add bilinear constraints later
    nh[0] = 0;
    nsh[0] = 0;
    ns[0] = nsh[0];
    ny[0] = 4;
    nz[0] = 0;

    for (int i = 1; i < NN; i++)
    {
        nx[i] = nx_;
        nu[i] = nu_;
        nbx[i] = 3;
        nbu[i] = nu_;
        ng[i] = 0;
        nh[i] = 1;
        nsh[i] = 1;
        ns[i] = nsh[i];
        ny[i] = 4;
        nz[i] = 0;
    }

    nx[NN] = nx_;
    nu[NN] = 0;
    nbx[NN] = 3;
    nbu[NN] = 0;
    ng[NN] = 0;
    nh[NN] = 0;
    nsh[NN] = 0;
    ns[NN] = nsh[NN];
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

    // middle stage
    int *idxsh1 = (int *) malloc(nsh[1]*sizeof(int));
    double *lsh1 = (double *) malloc((nsh[1])*sizeof(double));
    double *ush1 = (double *) malloc((nsh[1])*sizeof(double));


    /* box constraints */

    // acados inf
    double acados_inf = 1e8;

    // first stage

    // input bounds
    int *idxbu0 = (int *) malloc(nbu[0]*sizeof(int));
    double *lbu0 = (double *) malloc((nbu[0])*sizeof(double));
    double *ubu0 = (double *) malloc((nbu[0])*sizeof(double));

    // pitch angle rate
    idxbu0[0] = 0;
    lbu0[0] = dbeta_min;
    ubu0[0] = dbeta_max;

    // generator torque
    idxbu0[1] = 1;
    lbu0[1] = dM_gen_min;
    ubu0[1] = dM_gen_max;

    // state bounds
    int *idxbx0 = (int *) malloc(nbx[0]*sizeof(int));
    double *lbx0 = (double *) malloc((nbx[0])*sizeof(double));
    double *ubx0 = (double *) malloc((nbx[0])*sizeof(double));

    // dummy
    for (int ii=0; ii < nbx[0]; ii++)
    {
        idxbx0[ii] = ii;
        lbx0[ii] = - acados_inf;
        ubx0[ii] =   acados_inf;
    }


    // middle stages

    // input bounds
    int *idxbu1 = (int *) malloc(nbu[1]*sizeof(int));
    double *lbu1 = (double *) malloc((nbu[1])*sizeof(double));
    double *ubu1 = (double *) malloc((nbu[1])*sizeof(double));

    // pitch angle rate
    idxbu1[0] = 0;
    lbu1[0] = dbeta_min;
    ubu1[0] = dbeta_max;

    // generator torque rate
    idxbu1[1] = 1;
    lbu1[1] = dM_gen_min;
    ubu1[1] = dM_gen_max;

    // state bounds
    int *idxbx1 = (int *) malloc(nbx[1]*sizeof(int));
    double *lbx1 = (double *) malloc((nbx[1])*sizeof(double));
    double *ubx1 = (double *) malloc((nbx[1])*sizeof(double));

    // generator angular velocity
    idxbx1[0] = 0;
    lbx1[0] = OmegaR_min;
    ubx1[0] = OmegaR_max;

    // pitch angle
    idxbx1[1] = 6;
    lbx1[1] = beta_min;
    ubx1[1] = beta_max;

    // generator torque
    idxbx1[2] = 7;
    lbx1[2] = M_gen_min;
    ubx1[2] = M_gen_max;

    // last stage
    // state bounds
    int *idxbxN = (int *) malloc(nbx[NN]*sizeof(int));
    double *lbxN = (double *) malloc((nbx[NN])*sizeof(double));
    double *ubxN = (double *) malloc((nbx[NN])*sizeof(double));

    // generator angular velocity
    idxbxN[0] = 0;
    lbxN[0] = OmegaR_min;
    ubxN[0] = OmegaR_max;

    // pitch angle
    idxbxN[1] = 6;
    lbxN[1] = beta_min;
    ubxN[1] = beta_max;

    // generator torque
    idxbxN[2] = 7;
    lbxN[2] = M_gen_min;
    ubxN[2] = M_gen_max;

    // to shift
    double *specific_u = (double *) malloc(nu_*sizeof(double));
    double *specific_x = (double *) malloc(nx_*sizeof(double));

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
    if (nsh[1] > 0)
    {
        idxsh1[0] = 0;
        lsh1[0] = 0.0;
        ush1[0] = 0.0;
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

    double *VxN = (double *) malloc((ny[NN]*nx[NN])*sizeof(double));
    for (int ii=0; ii < ny[NN]*nx[NN]; ii++)
        VxN[ii] = 0.0;
    VxN[0+ny[NN]*0] = 1.0;
    VxN[1+ny[NN]*4] = 1.0;

    double *W = (double *) malloc((ny_*ny_)*sizeof(double));
    for (int ii = 0; ii < ny_*ny_; ii++)
        W[ii] = 0.0;
    W[0+ny_*0] = 1.5114;
    W[1+ny_*0] = -0.0649;
    W[0+ny_*1] = -0.0649;
    W[1+ny_*1] = 0.0180;
    W[2+ny_*2] = 0.01;
    W[3+ny_*3] = 0.001;

    double *W_N = (double *) malloc((ny[NN]*ny[NN])*sizeof(double));
    W_N[0+ny[NN]*0] = 1.5114;
    W_N[1+ny[NN]*0] = -0.0649;
    W_N[0+ny[NN]*1] = -0.0649;
    W_N[1+ny[NN]*1] = 0.0180;

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

    ocp_nlp_plan_t *plan = ocp_nlp_plan_create(NN);

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

    ocp_nlp_config *config = ocp_nlp_config_create(*plan);

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
        ocp_nlp_dims_set_cost(config, dims, i, "ny", &ny[i]);

        ocp_nlp_dims_set_constraints(config, dims, i, "nbx", &nbx[i]);
        ocp_nlp_dims_set_constraints(config, dims, i, "nbu", &nbu[i]);
        ocp_nlp_dims_set_constraints(config, dims, i, "ng", &ng[i]);
        ocp_nlp_dims_set_constraints(config, dims, i, "nh", &nh[i]);
        ocp_nlp_dims_set_constraints(config, dims, i, "nsh", &nsh[i]);

    }

    /************************************************
    * dynamics
    ************************************************/

    // explicit model
    external_function_param_casadi *expl_vde_for =
    (external_function_param_casadi *) malloc(NN*external_function_param_casadi_struct_size());
    // implicit model
    external_function_param_casadi *impl_ode_fun =
        (external_function_param_casadi *) malloc(NN*external_function_param_casadi_struct_size());
    external_function_param_casadi *impl_ode_fun_jac_x_xdot =
        (external_function_param_casadi *) malloc(NN*external_function_param_casadi_struct_size());
    external_function_param_casadi *impl_ode_jac_x_xdot_u =
        (external_function_param_casadi *) malloc(NN*external_function_param_casadi_struct_size());
    external_function_param_casadi *impl_ode_fun_jac_x_xdot_u =
        (external_function_param_casadi *) malloc(NN*external_function_param_casadi_struct_size());
    // gnsf model
    external_function_param_casadi *phi_fun =
        (external_function_param_casadi *) malloc(NN*external_function_param_casadi_struct_size());
    external_function_param_casadi *phi_fun_jac_y =
        (external_function_param_casadi *) malloc(NN*external_function_param_casadi_struct_size());
    external_function_param_casadi *phi_jac_y_uhat =
        (external_function_param_casadi *) malloc(NN*external_function_param_casadi_struct_size());
    external_function_param_casadi *f_lo_jac_x1_x1dot_u_z =
        (external_function_param_casadi *) malloc(NN*external_function_param_casadi_struct_size());

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

    /* initialize additional gnsf dimensions */
    int gnsf_nx1 = 8;
    int gnsf_nz1 = 0;
    int gnsf_nout = 1;
    int gnsf_ny = 5;
    int gnsf_nuhat = 0;

    for (int i = 0; i < NN; i++)
    {
        if (plan->sim_solver_plan[i].sim_solver == GNSF)
        {
            ocp_nlp_dims_set_dynamics(config, dims, i, "gnsf_nx1", &gnsf_nx1);
            ocp_nlp_dims_set_dynamics(config, dims, i, "gnsf_nz1", &gnsf_nz1);
            ocp_nlp_dims_set_dynamics(config, dims, i, "gnsf_nout", &gnsf_nout);
            ocp_nlp_dims_set_dynamics(config, dims, i, "gnsf_ny", &gnsf_ny);
            ocp_nlp_dims_set_dynamics(config, dims, i, "gnsf_nuhat", &gnsf_nuhat);
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
    for (int i = 0; i <= NN; i++)
    {
        // Cyt
        ocp_nlp_cost_model_set(config, dims, nlp_in, i, "Vu", Vu);
        if (i < NN)
            ocp_nlp_cost_model_set(config, dims, nlp_in, i, "Vx", Vx);
        else
            ocp_nlp_cost_model_set(config, dims, nlp_in, i, "Vx", VxN);
        // W
        ocp_nlp_cost_model_set(config, dims, nlp_in, i, "W", W);
    }
    ocp_nlp_cost_model_set(config, dims, nlp_in, NN, "W", W_N);

    // slacks (middle stages)
    for (int ii = 1; ii < NN; ii++)
    {
        ocp_nlp_cost_model_set(config, dims, nlp_in, ii, "Zl", lZ1);
        ocp_nlp_cost_model_set(config, dims, nlp_in, ii, "Zu", uZ1);
        ocp_nlp_cost_model_set(config, dims, nlp_in, ii, "zl", lz1);
        ocp_nlp_cost_model_set(config, dims, nlp_in, ii, "zu", uz1);
    }


    /* dynamics */

    int set_fun_status;

    for (int i = 0; i < NN; i++)
    {
        if (plan->sim_solver_plan[i].sim_solver == ERK)
        {
            set_fun_status = ocp_nlp_dynamics_model_set(config, dims, nlp_in, i, "expl_vde_for",
                &expl_vde_for[i]);
            REQUIRE(set_fun_status == 0);
            if (set_fun_status != 0) exit(1);
        }
        else if (plan->sim_solver_plan[i].sim_solver == IRK)
        {
            set_fun_status = ocp_nlp_dynamics_model_set(config, dims, nlp_in, i, "impl_ode_fun",
                &impl_ode_fun[i]);
            REQUIRE(set_fun_status == 0);
            if (set_fun_status != 0) exit(1);
            set_fun_status =
                ocp_nlp_dynamics_model_set(config, dims, nlp_in, i, "impl_ode_fun_jac_x_xdot",
                                           &impl_ode_fun_jac_x_xdot[i]);
            REQUIRE(set_fun_status == 0);
            if (set_fun_status != 0) exit(1);
            set_fun_status = ocp_nlp_dynamics_model_set(config, dims, nlp_in, i, "impl_ode_jac_x_xdot_u",
                &impl_ode_jac_x_xdot_u[i]);
            REQUIRE(set_fun_status == 0);
            if (set_fun_status != 0) exit(1);
        }
        else if (plan->sim_solver_plan[i].sim_solver == GNSF)
        {
            set_fun_status = ocp_nlp_dynamics_model_set(config, dims, nlp_in, i, "phi_fun", &phi_fun[i]);
            REQUIRE(set_fun_status == 0);
            if (set_fun_status != 0) exit(1);
            set_fun_status = ocp_nlp_dynamics_model_set(config, dims, nlp_in, i, "phi_fun_jac_y",
                &phi_fun_jac_y[i]);
            REQUIRE(set_fun_status == 0);
            if (set_fun_status != 0) exit(1);
            set_fun_status = ocp_nlp_dynamics_model_set(config, dims, nlp_in, i, "phi_jac_y_uhat",
                &phi_jac_y_uhat[i]);
            REQUIRE(set_fun_status == 0);
            if (set_fun_status != 0) exit(1);
            set_fun_status = ocp_nlp_dynamics_model_set(config, dims, nlp_in, i, "f_lo_jac_x1_x1dot_u_z",
                &f_lo_jac_x1_x1dot_u_z[i]);
            REQUIRE(set_fun_status == 0);
            if (set_fun_status != 0) exit(1);
            set_fun_status = ocp_nlp_dynamics_model_set(config, dims, nlp_in, i, "get_gnsf_matrices",
                                                        &get_matrices_fun);
            REQUIRE(set_fun_status == 0);
            if (set_fun_status != 0) exit(1);
        }
        else if (plan->sim_solver_plan[i].sim_solver == LIFTED_IRK)
        {
            set_fun_status = ocp_nlp_dynamics_model_set(config, dims, nlp_in, i, "impl_ode_fun",
                &impl_ode_fun[i]);
            REQUIRE(set_fun_status == 0);
            if (set_fun_status != 0) exit(1);
            set_fun_status =
                ocp_nlp_dynamics_model_set(config, dims, nlp_in, i, "impl_ode_fun_jac_x_xdot_u",
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
    ocp_nlp_constraints_model_set(config, dims, nlp_in, 0, "idxbu", idxbu0);
    ocp_nlp_constraints_model_set(config, dims, nlp_in, 0, "lbu", lbu0);
    ocp_nlp_constraints_model_set(config, dims, nlp_in, 0, "ubu", ubu0);
    ocp_nlp_constraints_model_set(config, dims, nlp_in, 0, "idxbx", idxbx0);
    ocp_nlp_constraints_model_set(config, dims, nlp_in, 0, "lbx", lbx0);
    ocp_nlp_constraints_model_set(config, dims, nlp_in, 0, "ubx", ubx0);
    // middle stages
    for (int i = 1; i < NN; i++)
    {
        ocp_nlp_constraints_model_set(config, dims, nlp_in, i, "idxbu", idxbu1);
        ocp_nlp_constraints_model_set(config, dims, nlp_in, i, "lbu", lbu1);
        ocp_nlp_constraints_model_set(config, dims, nlp_in, i, "ubu", ubu1);
        ocp_nlp_constraints_model_set(config, dims, nlp_in, i, "idxbx", idxbx1);
        ocp_nlp_constraints_model_set(config, dims, nlp_in, i, "lbx", lbx1);
        ocp_nlp_constraints_model_set(config, dims, nlp_in, i, "ubx", ubx1);
    }
    // last stage
    ocp_nlp_constraints_model_set(config, dims, nlp_in, NN, "idxbx", idxbxN);
    ocp_nlp_constraints_model_set(config, dims, nlp_in, NN, "lbx", lbxN);
    ocp_nlp_constraints_model_set(config, dims, nlp_in, NN, "ubx", ubxN);


    /* nonlinear constraints */

    // middle stages
    for (int i = 1; i < NN; i++)
    {
        if (nh[i] > 0)
        {
            ocp_nlp_constraints_model_set(config, dims, nlp_in, i, "lh", lh1);
            ocp_nlp_constraints_model_set(config, dims, nlp_in, i, "uh", uh1);
            ocp_nlp_constraints_model_set(config, dims, nlp_in, i, "nl_constr_h_fun_jac", &h1);
        }
    }

    /* soft constraints */

    // middle stages
    for (int i = 1; i < NN; i++)
    {
        if (ns[i] > 0)
        {
            ocp_nlp_constraints_model_set(config, dims, nlp_in, i, "lsh", lsh1);
            ocp_nlp_constraints_model_set(config, dims, nlp_in, i, "ush", ush1);
            ocp_nlp_constraints_model_set(config, dims, nlp_in, i, "idxsh", idxsh1);
        }
    }

    /************************************************
    * sqp opts
    ************************************************/

    void *nlp_opts = ocp_nlp_solver_opts_create(config, dims);

    // sim opts
    for (int i = 0; i < NN; ++i)
    {

        if (plan->sim_solver_plan[i].sim_solver == ERK)
        {
            int ns = 4;
            int num_steps = 10;
            ocp_nlp_solver_opts_set_at_stage(config, nlp_opts, i, "dynamics_num_steps", &num_steps);
            ocp_nlp_solver_opts_set_at_stage(config, nlp_opts, i, "dynamics_ns", &ns);
        }
        else if (plan->sim_solver_plan[i].sim_solver == IRK)
        {
            int num_steps = 1;
            int ns = 4;
            bool jac_reuse = true;

            ocp_nlp_solver_opts_set_at_stage(config, nlp_opts, i, "dynamics_num_steps", &num_steps);
            ocp_nlp_solver_opts_set_at_stage(config, nlp_opts, i, "dynamics_ns", &ns);
            ocp_nlp_solver_opts_set_at_stage(config, nlp_opts, i, "dynamics_jac_reuse", &jac_reuse);
        }
        else if (plan->sim_solver_plan[i].sim_solver == LIFTED_IRK)
        {
            int num_steps = 1;
            int ns = 4;

            ocp_nlp_solver_opts_set_at_stage(config, nlp_opts, i, "dynamics_num_steps", &num_steps);
            ocp_nlp_solver_opts_set_at_stage(config, nlp_opts, i, "dynamics_ns", &ns);
        }
        else if (plan->sim_solver_plan[i].sim_solver == GNSF)
        {
            int num_steps = 1;
            int ns = 4;
            int newton_iter = 1;
            bool jac_reuse = true;

            ocp_nlp_solver_opts_set_at_stage(config, nlp_opts, i, "dynamics_num_steps", &num_steps);
            ocp_nlp_solver_opts_set_at_stage(config, nlp_opts, i, "dynamics_ns", &ns);
            ocp_nlp_solver_opts_set_at_stage(config, nlp_opts, i, "dynamics_jac_reuse", &jac_reuse);
            ocp_nlp_solver_opts_set_at_stage(config, nlp_opts, i, "dynamics_newton_iter", &newton_iter);
        }
    }

    int max_iter = MAX_SQP_ITERS;
    double tol_stat = 1e-6;
    double tol_eq   = 1e-8;
    double tol_ineq = 1e-8;
    double tol_comp = 1e-8;

    ocp_nlp_solver_opts_set(config, nlp_opts, "max_iter", &max_iter);
    ocp_nlp_solver_opts_set(config, nlp_opts, "tol_stat", &tol_stat);
    ocp_nlp_solver_opts_set(config, nlp_opts, "tol_eq", &tol_eq);
    ocp_nlp_solver_opts_set(config, nlp_opts, "tol_ineq", &tol_ineq);
    ocp_nlp_solver_opts_set(config, nlp_opts, "tol_comp", &tol_comp);


    // partial condensing
    if (plan->ocp_qp_solver_plan.qp_solver == PARTIAL_CONDENSING_HPIPM)
    {
        int cond_N = 10;
        ocp_nlp_solver_opts_set(config, nlp_opts, "qp_cond_N", &cond_N);
    }

    config->opts_update(config, dims, nlp_opts);

    /************************************************
    * ocp_nlp out
    ************************************************/

    ocp_nlp_out *nlp_out = ocp_nlp_out_create(config, dims);

    ocp_nlp_solver *solver = ocp_nlp_solver_create(config, dims, nlp_opts);

    /************************************************
    *     precomputation (after all options are set)
    ************************************************/

    int status = ocp_nlp_precompute(solver, nlp_in, nlp_out);

    /************************************************
    * sqp solve
    ************************************************/

    int nmpc_problems = 10;

    acados_timer timer;
    acados_tic(&timer);

    // warm start output initial guess of solution
    for (int i = 0; i <= NN; i++)
    {
        blasfeo_pack_dvec(2, u0_ref, 1, nlp_out->ux+i, 0);
        blasfeo_pack_dvec(1, wind0_ref+i, 1, nlp_out->ux+i, 2);
        blasfeo_pack_dvec(nx[i], x0_ref, 1, nlp_out->ux+i, nu[i]);
    }

    // set x0 as box constraint
    ocp_nlp_constraints_model_set(config, dims, nlp_in, 0, "lbx", x0_ref);
    ocp_nlp_constraints_model_set(config, dims, nlp_in, 0, "ubx", x0_ref);

    double total_power = 0.0; // sum up to get total objective
    double electrical_power = 0.0;

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
            ocp_nlp_cost_model_set(config, dims, nlp_in, i, "yref", &y_ref[(idx + i)*4]);
        }

        // solve NLP
        status = ocp_nlp_solve(solver, nlp_in, nlp_out);

        ocp_nlp_res *residual;
        ocp_nlp_get(config, solver, "nlp_res", &residual);

        double max_res = 0.0;
        double inf_norm_res_stat = residual->inf_norm_res_stat;
        double inf_norm_res_eq = residual->inf_norm_res_eq;
        double inf_norm_res_ineq = residual->inf_norm_res_ineq;
        double inf_norm_res_comp = residual->inf_norm_res_comp;
        max_res = (inf_norm_res_stat > max_res) ? inf_norm_res_stat : max_res;
        max_res = (inf_norm_res_eq > max_res) ? inf_norm_res_eq : max_res;
        max_res = (inf_norm_res_ineq > max_res) ? inf_norm_res_ineq : max_res;
        max_res = (inf_norm_res_comp > max_res) ? inf_norm_res_comp : max_res;

        // update initial condition
        // TODO(dimitris): maybe simulate system instead of passing x[1] as next state
        ocp_nlp_out_get(config, dims, nlp_out, 1, "x", specific_x);
        ocp_nlp_constraints_model_set(config, dims, nlp_in, 0, "lbx", specific_x);
        ocp_nlp_constraints_model_set(config, dims, nlp_in, 0, "ubx", specific_x);

        int sqp_iter;
        double time_lin, time_qp_sol, time_tot;

        ocp_nlp_get(config, solver, "sqp_iter", &sqp_iter);
        ocp_nlp_get(config, solver, "time_tot", &time_tot);
        ocp_nlp_get(config, solver, "time_qp_sol", &time_qp_sol);
        ocp_nlp_get(config, solver, "time_lin", &time_lin);

        printf("\nproblem #%d, status %d, iters %d, time (total %f, lin %f, qp_sol %f) ms\n",
            idx, status, sqp_iter, time_tot*1e3, time_lin*1e3, time_qp_sol*1e3);

        // TODO(oj): rather simulate then take x from SQP solution
        printf("xsim = \n");
        ocp_nlp_out_get(config, dims, nlp_out, 0, "x", x_end);
        d_print_mat(1, nx[0], x_end, 1);

        electrical_power = 0.944*97/100* x_end[0] * x_end[5];
        total_power += electrical_power;
        printf("electrical power = %f\n", electrical_power);
        printf("Max residuals = %e\n", max_res);

        REQUIRE((status == 0 || status == 1 && MAX_SQP_ITERS == 1));
        REQUIRE(max_res <= TOL);

        // shift trajectories
        blasfeo_unpack_dvec(dims->nx[NN], &nlp_out->ux[NN-1], dims->nu[NN-1], x_end, 1);
        blasfeo_unpack_dvec(dims->nu[NN-1], &nlp_out->ux[NN-2], dims->nu[NN-2], u_end, 1);

        shift_states(dims, nlp_out, x_end);
        shift_controls(dims, nlp_out, u_end);
    }

    double time = acados_toc(&timer);

    printf("\n\ntotal time (including printing) = %f ms (time per SQP = %f)\n\n",
        time*1e3, time*1e3/nmpc_problems);
    printf("total electrical power %e\n\n", total_power);
    REQUIRE(total_power - 5.161e1 > 0); // ensure MPC has a performance close to the known optimum



    /************************************************
    * free memory
    ************************************************/

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

    ocp_nlp_solver_opts_destroy(nlp_opts);
    ocp_nlp_in_destroy(nlp_in);
    ocp_nlp_out_destroy(nlp_out);
    ocp_nlp_solver_destroy(solver);
    ocp_nlp_dims_destroy(dims);
    ocp_nlp_config_destroy(config);
    ocp_nlp_plan_destroy(plan);

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

    free(W_N);
    free(W);
    free(VxN);
    free(Vx);
    free(Vu);
    free(lh1);
    free(uh1);

    free(idxbu0);
    free(lbu0);
    free(ubu0);
    free(idxbx0);
    free(lbx0);
    free(ubx0);

    free(idxbx1);
    free(lbu1);
    free(ubu1);
    free(idxbu1);
    free(lbx1);
    free(ubx1);

    free(idxbxN);
    free(lbxN);
    free(ubxN);

    free(idxsh1);
    free(lsh1);
    free(ush1);

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
