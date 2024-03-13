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


#include <stdio.h>
#include <stdlib.h>

#include "acados/utils/print.h"
#include "acados/ocp_nlp/ocp_nlp_constraints_bgp.h"
#include "acados/ocp_nlp/ocp_nlp_cost_ls.h"
#include "acados/ocp_nlp/ocp_nlp_dynamics_common.h"
#include "acados/ocp_nlp/ocp_nlp_dynamics_cont.h"
#include "acados/ocp_nlp/ocp_nlp_sqp.h"
#include "acados/sim/sim_common.h"
#include "acados/sim/sim_irk_integrator.h"

#include "acados_c/ocp_nlp_interface.h"
#include "acados_c/external_function_interface.h"

#include "blasfeo/include/blasfeo_d_aux.h"

#include "simple_dae_model/simple_dae_model.h"
#include "simple_dae_model/simple_dae_constr.h"

#define FORMULATION 1 
// 0: without Vz*z term 
// 1: with Vz*z and without Vx*x 
// 2: same as (1) + nonlinear constraint on z: h(x,u,z(x,u)) = [2, -2] \leq [z_1, z_2] \leq [4, 2]

int main() {

	double lh[2];
	double uh[2];

	lh[0] = 2;
	lh[1] = -2;

	uh[0] = 4;
	uh[1] = 2;

    int NH, NBX;
    if (FORMULATION == 2) {
		NH  = 2;
		NBX = 0;
	} else {
		NH  = 0;
		NBX = 2;
	}

	int num_states = 2, num_controls = 2, N = 20;
	int num_alg_states = 2;
	double Tf = 1.0, R[2] = {1e-3, 1e-3}, QN[2] = {1e1, 1e1};
    double Q[2] = {1e1, 1e1};
	int idxb_0[2] = {2, 3};
	int idxb[2] = {2, 3};
	double x0[num_states];

    x0[0] =  3;  
    x0[1] =  -1.8;  

	int max_num_sqp_iterations = 100;

    int nx_ = num_states;
    int nz_ = num_alg_states;
    int nu_ = num_controls;
    int ny_;

    if (FORMULATION == 0) {
        ny_ = nu_ + nx_;
    } else {
        ny_ = nu_ + nz_;
    }

	int nx[N+1];
    int nu[N+1];
    int nbx[N+1];
    int nbu[N+1];
    int ng[N+1];
    int nh[N+1];
	int	ns[N+1];
    int nz[N+1];
    int ny[N+1];

    for(int i = 0; i < N+1; i++) {
        nx[i] = nx_;
        nu[i] = nu_;
        nbx[i] = NBX;
        nbu[i] = 0;
        ng[i] = 0;
        nh[i] = NH;
        ns[i] = 0;
        nz[i] = nz_;
        ny[i] = ny_;
    }

    nbx[0] = nx_;
    nbx[N] = 0;
    nu[N]  = 0;
    nh[N]  = 0;
    ny[N]  = nx_;
    nz[N]  = 0;

    /* linear least squares */

    // output definition
    // y  = Vx * x + Vu * u + Vz * z

    double *Vx = malloc((ny_*nx_)*sizeof(double));
    for (int ii=0; ii < ny_*nx_; ii++)
        Vx[ii] = 0.0;

    if (FORMULATION == 0) {
        Vx[0+ny_*0] = 1.0;
        Vx[1+ny_*1] = 1.0;
    } else {
        Vx[0+ny_*0] = 0.0;
        Vx[1+ny_*1] = 0.0;
    }

    double *Vu = malloc((ny_*nu_)*sizeof(double));
    for (int ii=0; ii < ny_*nu_; ii++)
        Vu[ii] = 0.0;

    Vu[2+ny_*0] = 1.0;
    Vu[3+ny_*1] = 1.0;

    double *Vz = malloc((ny_*nz_)*sizeof(double));
    for (int ii=0; ii < ny_*nz_; ii++)
        Vz[ii] = 0.0;

    if (FORMULATION == 0) {
        Vz[0+ny_*0] = 0.0;
        Vz[1+ny_*1] = 0.0;
    } else {
        Vz[0+ny_*0] = 1.0;
        Vz[1+ny_*1] = 1.0;
    }

    double *VxN = malloc((ny[N]*nx_)*sizeof(double));
    for (int ii=0; ii < ny[N]*nx_; ii++)
        VxN[ii] = 0.0;

    VxN[0+ny[N]*0] = 1.0;
    VxN[1+ny[N]*1] = 1.0;

    double *W = malloc((ny_*ny_)*sizeof(double));
    for (int ii=0; ii<ny_*ny_; ii++)
        W[ii] = 0.0;

    W[0+ny_*0] = Q[0];
    W[1+ny_*1] = Q[1];
    W[2+ny_*2] = R[0];
    W[3+ny_*3] = R[1];

    double *WN = malloc((ny[N]*ny[N])*sizeof(double));
    for (int ii=0; ii<ny[N]*ny[N]; ii++)
        WN[ii] = 0.0;

    WN[0+ny[N]*0] = QN[0];
    WN[1+ny[N]*1] = QN[1];

	// Make plan
	ocp_nlp_plan_t *plan = ocp_nlp_plan_create(N);
	plan->nlp_solver = SQP;
	plan->ocp_qp_solver_plan.qp_solver = PARTIAL_CONDENSING_HPIPM;

//	plan->ocp_qp_solver_plan.qp_solver = FULL_CONDENSING_QPOASES;

	for (int i = 0; i <= N; i++)
		plan->nlp_cost[i] = LINEAR_LS;
	for (int i = 0; i < N; i++)
	{
		plan->nlp_dynamics[i] = CONTINUOUS_MODEL;
		plan->sim_solver_plan[i].sim_solver = IRK;
	}

	for (int i = 0; i <= N; i++)
		plan->nlp_constraints[i] = BGH;

	ocp_nlp_config *config = ocp_nlp_config_create(*plan);

	ocp_nlp_dims *dims = ocp_nlp_dims_create(config);

    ocp_nlp_dims_set_opt_vars(config, dims, "nx", nx);
    ocp_nlp_dims_set_opt_vars(config, dims, "nu", nu);
    ocp_nlp_dims_set_opt_vars(config, dims, "nz", nz);
    ocp_nlp_dims_set_opt_vars(config, dims, "ns", ns);

    for (int i = 0; i <= N; i++) {
        ocp_nlp_dims_set_cost(config, dims, i, "ny", &ny[i]);
        ocp_nlp_dims_set_constraints(config, dims, i, "nbx", &nbx[i]);
        ocp_nlp_dims_set_constraints(config, dims, i, "nbu", &nbu[i]);
        ocp_nlp_dims_set_constraints(config, dims, i, "ng", &ng[i]);
        ocp_nlp_dims_set_constraints(config, dims, i, "nh", &nh[i]);
        ocp_nlp_dims_set_constraints(config, dims, i, "nz", &nz[i]);
    }

	external_function_param_casadi impl_ode_fun[N];
	external_function_param_casadi impl_ode_fun_jac_x_xdot_z[N];
	external_function_param_casadi impl_ode_jac_x_xdot_u_z[N];

	for (int ii = 0; ii < N; ++ii) {
        impl_ode_fun[ii].casadi_fun = &simple_dae_impl_ode_fun;
        impl_ode_fun[ii].casadi_work = &simple_dae_impl_ode_fun_work;
        impl_ode_fun[ii].casadi_sparsity_in = &simple_dae_impl_ode_fun_sparsity_in;
        impl_ode_fun[ii].casadi_sparsity_out = &simple_dae_impl_ode_fun_sparsity_out;
        impl_ode_fun[ii].casadi_n_in = &simple_dae_impl_ode_fun_n_in;
        impl_ode_fun[ii].casadi_n_out = &simple_dae_impl_ode_fun_n_out;

        impl_ode_fun_jac_x_xdot_z[ii].casadi_fun = &simple_dae_impl_ode_fun_jac_x_xdot_z;
        impl_ode_fun_jac_x_xdot_z[ii].casadi_work = &simple_dae_impl_ode_fun_jac_x_xdot_z_work;
        impl_ode_fun_jac_x_xdot_z[ii].casadi_sparsity_in = &simple_dae_impl_ode_fun_jac_x_xdot_z_sparsity_in;
        impl_ode_fun_jac_x_xdot_z[ii].casadi_sparsity_out = &simple_dae_impl_ode_fun_jac_x_xdot_z_sparsity_out;
        impl_ode_fun_jac_x_xdot_z[ii].casadi_n_in = &simple_dae_impl_ode_fun_jac_x_xdot_z_n_in;
        impl_ode_fun_jac_x_xdot_z[ii].casadi_n_out = &simple_dae_impl_ode_fun_jac_x_xdot_z_n_out;

        impl_ode_jac_x_xdot_u_z[ii].casadi_fun = &simple_dae_impl_ode_jac_x_xdot_u_z;
        impl_ode_jac_x_xdot_u_z[ii].casadi_work = &simple_dae_impl_ode_jac_x_xdot_u_z_work;
        impl_ode_jac_x_xdot_u_z[ii].casadi_sparsity_in = &simple_dae_impl_ode_jac_x_xdot_u_z_sparsity_in;
        impl_ode_jac_x_xdot_u_z[ii].casadi_sparsity_out = &simple_dae_impl_ode_jac_x_xdot_u_z_sparsity_out;
        impl_ode_jac_x_xdot_u_z[ii].casadi_n_in = &simple_dae_impl_ode_jac_x_xdot_u_z_n_in;
        impl_ode_jac_x_xdot_u_z[ii].casadi_n_out = &simple_dae_impl_ode_jac_x_xdot_u_z_n_out;
	}

	// impl_ode
    int	tmp_size = 0;
	for (int ii=0; ii<N; ii++)
	{
		tmp_size += external_function_param_casadi_calculate_size(impl_ode_fun+ii, 0);
	}
	void *impl_ode_casadi_mem = malloc(tmp_size);
	void *c_ptr = impl_ode_casadi_mem;
	for (int ii=0; ii<N; ii++)
	{
		external_function_param_casadi_assign(impl_ode_fun+ii, c_ptr);
		c_ptr += external_function_param_casadi_calculate_size(impl_ode_fun+ii, 0);
	}
	//
	tmp_size = 0;
	for (int ii=0; ii<N; ii++)
	{
		tmp_size += external_function_param_casadi_calculate_size(impl_ode_fun_jac_x_xdot_z+ii, 0);
	}
	void *impl_ode_fun_jac_x_xdot_z_mem = malloc(tmp_size);
	c_ptr = impl_ode_fun_jac_x_xdot_z_mem;
	for (int ii=0; ii<N; ii++)
	{
		external_function_param_casadi_assign(impl_ode_fun_jac_x_xdot_z+ii, c_ptr);
		c_ptr += external_function_param_casadi_calculate_size(impl_ode_fun_jac_x_xdot_z+ii, 0);
	}

	tmp_size = 0;
	for (int ii=0; ii<N; ii++)
	{
		tmp_size += external_function_param_casadi_calculate_size(impl_ode_jac_x_xdot_u_z+ii, 0);
	}
	void *impl_ode_jac_x_xdot_u_z_mem = malloc(tmp_size);
	c_ptr = impl_ode_jac_x_xdot_u_z_mem;
	for (int ii=0; ii<N; ii++)
	{
		external_function_param_casadi_assign(impl_ode_jac_x_xdot_u_z+ii, c_ptr);
		c_ptr += external_function_param_casadi_calculate_size(impl_ode_jac_x_xdot_u_z+ii, 0);
	}

	ocp_nlp_in *nlp_in = ocp_nlp_in_create(config, dims);

	for (int i = 0; i < N; ++i)
		nlp_in->Ts[i] = Tf/N;

	// NLP cost: linear least squares
    // C
	ocp_nlp_cost_ls_model **cost_ls = (ocp_nlp_cost_ls_model **) nlp_in->cost;

	for (int i = 0; i < N; ++i) {
        ocp_nlp_cost_model_set(config, dims, nlp_in, i, "Vx", Vx);
        ocp_nlp_cost_model_set(config, dims, nlp_in, i, "Vu", Vu);
        ocp_nlp_cost_model_set(config, dims, nlp_in, i, "Vz", Vz);
	}

    ocp_nlp_cost_model_set(config, dims, nlp_in, N, "Vx", VxN);
    // ocp_nlp_cost_model_set(config, dims, nlp_in, N, "Vz", Vz);
    
	// W
	for (int i = 0; i < N; ++i) ocp_nlp_cost_model_set(config, dims, nlp_in, i, "W", W);
    
	// WN
    ocp_nlp_cost_model_set(config, dims, nlp_in, N, "W", WN);

	// y_ref
    for (int i = 0; i <= N; ++i)
		blasfeo_dvecse(ny[i], 0.0, &cost_ls[i]->y_ref, 0);

	// NLP dynamics
	for (int i = 0; i < N; ++i) {
		ocp_nlp_dynamics_cont_model *dynamics = nlp_in->dynamics[i];
		irk_model *model = dynamics->sim_model;
		model->impl_ode_fun = (external_function_generic *) &impl_ode_fun[i];
		model->impl_ode_fun_jac_x_xdot_z = (external_function_generic *) &impl_ode_fun_jac_x_xdot_z[i];
		model->impl_ode_jac_x_xdot_u_z = (external_function_generic *) &impl_ode_jac_x_xdot_u_z[i]; 
	}

	// bounds
	ocp_nlp_constraints_bgh_model **constraints = (ocp_nlp_constraints_bgh_model **) nlp_in->constraints;

	ocp_nlp_constraints_model_set(config, dims, nlp_in, 0, "ubx", x0);
	ocp_nlp_constraints_model_set(config, dims, nlp_in, 0, "lbx", x0);
    constraints[0]->idxb = idxb_0;

    if (FORMULATION == 2) {
        external_function_param_casadi * nl_constr_h_fun_jac;
		nl_constr_h_fun_jac = (external_function_param_casadi *) malloc(sizeof(external_function_param_casadi)*N);

        for (int ii = 0; ii < N; ++ii) {
            nl_constr_h_fun_jac[ii].casadi_fun          = &simple_dae_constr_h_fun_jac_ut_xt;
            nl_constr_h_fun_jac[ii].casadi_work         = &simple_dae_constr_h_fun_jac_ut_xt_work;
            nl_constr_h_fun_jac[ii].casadi_sparsity_in  = &simple_dae_constr_h_fun_jac_ut_xt_sparsity_in;
            nl_constr_h_fun_jac[ii].casadi_sparsity_out = &simple_dae_constr_h_fun_jac_ut_xt_sparsity_out;
            nl_constr_h_fun_jac[ii].casadi_n_in         = &simple_dae_constr_h_fun_jac_ut_xt_n_in;
            nl_constr_h_fun_jac[ii].casadi_n_out        = &simple_dae_constr_h_fun_jac_ut_xt_n_out;
            external_function_param_casadi_create(&nl_constr_h_fun_jac[ii], 0);
            ocp_nlp_constraints_model_set(config, dims, nlp_in, ii, "lh", lh);
            ocp_nlp_constraints_model_set(config, dims, nlp_in, ii, "uh", uh);
			ocp_nlp_constraints_model_set(config, dims, nlp_in, ii, "nl_constr_h_fun_jac", &nl_constr_h_fun_jac[ii]);
        }
    } else { 
		for (int ii = 1; ii < N; ++ii) {
			ocp_nlp_constraints_model_set(config, dims, nlp_in, ii, "ubx", uh);
			ocp_nlp_constraints_model_set(config, dims, nlp_in, ii, "lbx", lh);
			constraints[ii]->idxb = idxb;
		}
	}


	void *nlp_opts = ocp_nlp_solver_opts_create(config, dims);
   
    bool output_z_val = true; 
    bool sens_algebraic_val = true; 
    bool reuse_val = true; 
    int num_steps_val = 5; 
    for (int i = 0; i < N; i++) ocp_nlp_solver_opts_set_at_stage(config, nlp_opts, i, "dynamics_output_z", &output_z_val);
    for (int i = 0; i < N; i++) ocp_nlp_solver_opts_set_at_stage(config, nlp_opts, i, "dynamics_sens_algebraic", &sens_algebraic_val);
    for (int i = 0; i < N; i++) ocp_nlp_solver_opts_set_at_stage(config, nlp_opts, i, "dynamics_jac_reuse", &reuse_val);
    for (int i = 0; i < N; i++) ocp_nlp_solver_opts_set_at_stage(config, nlp_opts, i, "dynamics_num_steps", &num_steps_val);

	double tol_stat = 1e-6;
	double tol_eq   = 1e-6;
	double tol_ineq = 1e-6;
	double tol_comp = 1e-6;
	ocp_nlp_solver_opts_set(config, nlp_opts, "max_iter", &max_num_sqp_iterations);
	ocp_nlp_solver_opts_set(config, nlp_opts, "tol_stat", &tol_stat);
	ocp_nlp_solver_opts_set(config, nlp_opts, "tol_eq", &tol_eq);
	ocp_nlp_solver_opts_set(config, nlp_opts, "tol_ineq", &tol_ineq);
	ocp_nlp_solver_opts_set(config, nlp_opts, "tol_comp", &tol_comp);

    int print_level = 1;
    ocp_nlp_solver_opts_set(config, nlp_opts, "print_level", &print_level);

	if(plan->ocp_qp_solver_plan.qp_solver == PARTIAL_CONDENSING_HPIPM)
	{
		int N2 = N;
		ocp_nlp_solver_opts_set(config, nlp_opts, "qp_cond_N", &N2);
	}

	ocp_nlp_out *nlp_out = ocp_nlp_out_create(config, dims);
	for (int i = 0; i <= N; ++i) {
		blasfeo_dvecse(nu[i]+nx[i], 0.0, nlp_out->ux+i, 0);
		blasfeo_dvecse(2, 0.0, nlp_out->ux+i, 0);
    }

	ocp_nlp_solver *solver = ocp_nlp_solver_create(config, dims, nlp_opts);
    ocp_nlp_precompute(solver, nlp_in, nlp_out);

	// NLP solution
    acados_timer timer;
    acados_tic(&timer);

	int solver_status = ocp_nlp_solve(solver, nlp_in, nlp_out);

    double elapsed_time = acados_toc(&timer);

	printf("\nsolution\n");
	ocp_nlp_out_print(dims, nlp_out);

    int sqp_iter;
    ocp_nlp_get(config, solver, "sqp_iter", &sqp_iter);
    printf("\n\nstatus = %i, avg time = %f ms, iters = %d\n\n", solver_status, elapsed_time, sqp_iter);

    free(Vx);
    free(Vu);
    free(Vz);
    free(VxN);
    free(W);
    free(WN);

	return solver_status;
}
