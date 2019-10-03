/*
 * Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
 * Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
 * Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
 * Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
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
#include <cstdlib>
#include <vector>

#include "acados/utils/print.h"
#include "acados/ocp_nlp/ocp_nlp_constraints_bghp.h"
#include "acados/ocp_nlp/ocp_nlp_cost_ls.h"
#include "acados/ocp_nlp/ocp_nlp_dynamics_common.h"
#include "acados/ocp_nlp/ocp_nlp_dynamics_cont.h"
#include "acados/ocp_nlp/ocp_nlp_sqp.h"
#include "acados/sim/sim_erk_integrator.h"

#include "acados_c/ocp_nlp_interface.h"

#include "blasfeo/include/blasfeo_d_aux.h"

#include "pendulum_model/pendulum_model.h"
#include "pendulum_model/constraint.h"
#include "pendulum_model/position.h"

#define PI 3.1415926535897932

int main() {

	int num_states = 4, num_controls = 1, N = 20;
	double Tf = 1.0, Q = 1e-10, R = 1e-4, QN = 1e-10;
	std::vector<int> idxb_0 {1, 2, 3, 4};
	std::vector<double> x0(num_states, 0);
	x0[2] = PI;

	double radius2 = 0.04, neg_inf = -1000000;
	int max_num_sqp_iterations = 100;

	std::vector<int> nx(N+1, num_states), nu(N+1, num_controls), nbx(N+1, 0), nbu(N+1, 0),
		nb(N+1, 0), ng(N+1, 0), nh(N+1, 0), np(N+1, 0),
		ns(N+1, 0), nz(N+1, 0), nv(N+1, num_states+num_controls), ny(N+1, num_states+num_controls);

	nbx.at(0) = num_states;
	nb.at(0) = num_states;

	nu.at(N) = 0;
	nv.at(N) = 4;
	ny.at(N) = 4;
	nh.at(N) = 1;
	np.at(N) = 2;

	// Make plan
	ocp_nlp_plan *plan = ocp_nlp_plan_create(N);
	plan->nlp_solver = SQP;
	plan->ocp_qp_solver_plan.qp_solver = PARTIAL_CONDENSING_HPIPM;
	for (int i = 0; i <= N; i++)
		plan->nlp_cost[i] = LINEAR_LS;
	for (int i = 0; i < N; i++)
	{
		plan->nlp_dynamics[i] = CONTINUOUS_MODEL;
		plan->sim_solver_plan[i].sim_solver = ERK;
	}
	for (int i = 0; i <= N; i++)
		plan->nlp_constraints[i] = BGHP;

	ocp_nlp_config *config = ocp_nlp_config_create(*plan);

	ocp_nlp_dims *dims = ocp_nlp_dims_create(config);
    ocp_nlp_dims_set_opt_vars(config, dims, "nx", nx.data());
    ocp_nlp_dims_set_opt_vars(config, dims, "nu", nu.data());
    ocp_nlp_dims_set_opt_vars(config, dims, "nz", nz.data());
    ocp_nlp_dims_set_opt_vars(config, dims, "ns", ns.data());

    for (int i = 0; i <= N; i++)
    {
        ocp_nlp_dims_set_cost(config, dims, i, "ny", ny.data()+i);

        ocp_nlp_dims_set_constraints(config, dims, i, "nbx", nbx.data()+i);
        ocp_nlp_dims_set_constraints(config, dims, i, "nbu", nbu.data()+i);
        ocp_nlp_dims_set_constraints(config, dims, i, "ng", ng.data()+i);
        ocp_nlp_dims_set_constraints(config, dims, i, "nh", nh.data()+i);
        ocp_nlp_dims_set_constraints(config, dims, i, "nh", nh.data()+i);
        ocp_nlp_dims_set_constraints(config, dims, i, "np", np.data()+i);
    }



	external_function_casadi forw_vde_casadi[N];
	for (int i = 0; i < N; ++i) {
		forw_vde_casadi[i].casadi_fun = &pendulum_ode_expl_vde_forw;
		forw_vde_casadi[i].casadi_n_in = &pendulum_ode_expl_vde_forw_n_in;
		forw_vde_casadi[i].casadi_n_out = &pendulum_ode_expl_vde_forw_n_out;
		forw_vde_casadi[i].casadi_sparsity_in = &pendulum_ode_expl_vde_forw_sparsity_in;
		forw_vde_casadi[i].casadi_sparsity_out = &pendulum_ode_expl_vde_forw_sparsity_out;
		forw_vde_casadi[i].casadi_work = &pendulum_ode_expl_vde_forw_work;
	}

	// NLP model: forward VDEs
	int function_size = 0;
	for (int i = 0; i < N; ++i)
		function_size += external_function_casadi_calculate_size(forw_vde_casadi+i);

	char *c_ptr = (char *) calloc(1, function_size);
	for (int i = 0; i < N; ++i) {
		external_function_casadi_assign(forw_vde_casadi+i, c_ptr);
		c_ptr += external_function_casadi_calculate_size(forw_vde_casadi+i);
	}

	ocp_nlp_in *nlp_in = ocp_nlp_in_create(config, dims);

	for (int i = 0; i < N; ++i)
		nlp_in->Ts[i] = Tf/N;

	// NLP cost: linear least squares
    // C
	ocp_nlp_cost_ls_model **cost_ls = (ocp_nlp_cost_ls_model **) nlp_in->cost;
	for (int i = 0; i <= N; ++i) {
		blasfeo_dgese(nv[i], ny[i], 0.0, &cost_ls[i]->Cyt, 0, 0);
        for (int j = 0; j < nu[i]; j++)
            BLASFEO_DMATEL(&cost_ls[i]->Cyt, j, nx[i]+j) = 1.0;
        for (int j = 0; j < nx[i]; j++)
            BLASFEO_DMATEL(&cost_ls[i]->Cyt, nu[i]+j, j) = 1.0;
	}

	// W
	for (int i = 0; i < N; ++i) {
		blasfeo_dgese(ny[i], ny[i], 0.0, &cost_ls[i]->W, 0, 0);
        for (int j = 0; j < nx[i]; j++)
            BLASFEO_DMATEL(&cost_ls[i]->W, j, j) = Q;
        for (int j = 0; j < nu[i]; j++)
            BLASFEO_DMATEL(&cost_ls[i]->W, nx[i]+j, nx[i]+j) = R;
	}
	// WN
	blasfeo_dgese(ny[N], ny[N], 0.0, &cost_ls[N]->W, 0, 0);
	for (int j = 0; j < nx[N]; j++)
		BLASFEO_DMATEL(&cost_ls[N]->W, j, j) = QN;

	// y_ref
    for (int i = 0; i <= N; ++i)
		blasfeo_dvecse(ny[i], 0.0, &cost_ls[i]->y_ref, 0);

	// NLP dynamics
	for (int i = 0; i < N; ++i) {
		ocp_nlp_dynamics_cont_model *dynamics = (ocp_nlp_dynamics_cont_model *) nlp_in->dynamics[i];
		erk_model *model = (erk_model *) dynamics->sim_model;
		model->expl_vde_for = (external_function_generic *) &forw_vde_casadi[i];
	}

	// convex-composite constraint
	external_function_casadi nonlinear_constraint;
	nonlinear_constraint.casadi_fun = &constraint;
	nonlinear_constraint.casadi_n_in = &constraint_n_in;
	nonlinear_constraint.casadi_n_out = &constraint_n_out;
	nonlinear_constraint.casadi_sparsity_in = &constraint_sparsity_in;
	nonlinear_constraint.casadi_sparsity_out = &constraint_sparsity_out;
	nonlinear_constraint.casadi_work = &constraint_work;

	int constraint_size = external_function_casadi_calculate_size(&nonlinear_constraint);
	void *ptr = malloc(constraint_size);
	external_function_casadi_assign(&nonlinear_constraint, ptr);

	// nonlinear part of convex-composite constraint
	external_function_casadi position_constraint;
	position_constraint.casadi_fun = &position;
	position_constraint.casadi_n_in = &position_n_in;
	position_constraint.casadi_n_out = &position_n_out;
	position_constraint.casadi_sparsity_in = &position_sparsity_in;
	position_constraint.casadi_sparsity_out = &position_sparsity_out;
	position_constraint.casadi_work = &position_work;

	constraint_size = external_function_casadi_calculate_size(&position_constraint);
	ptr = malloc(constraint_size);
	external_function_casadi_assign(&position_constraint, ptr);

	// bounds
	ocp_nlp_constraints_bghp_model **constraints = (ocp_nlp_constraints_bghp_model **) nlp_in->constraints;

    constraints[0]->idxb = idxb_0.data();
	blasfeo_pack_dvec(nb[0], x0.data(), &constraints[0]->d, 0);
	blasfeo_pack_dvec(nb[0], x0.data(), &constraints[0]->d, nb[0]+ng[0]);

	blasfeo_pack_dvec(nh[N], &neg_inf, &constraints[N]->d, nb[N]+ng[N]);
	blasfeo_pack_dvec(nh[N], &radius2, &constraints[N]->d, 2*(nb[N]+ng[N])+nh[N]);
	constraints[N]->nl_constr_h_fun_jac = (external_function_generic *) &nonlinear_constraint;
	constraints[N]->p = (external_function_generic *) &position_constraint;

	void *nlp_opts = ocp_nlp_solver_opts_create(config, dims);

    int max_iter = max_num_sqp_iterations;
    double tol_stat = 1e-9;
    double tol_eq   = 1e-9;
    double tol_ineq = 1e-9;
    double tol_comp = 1e-9;

    ocp_nlp_solver_opts_set(config, nlp_opts, "max_iter", &max_iter);
    ocp_nlp_solver_opts_set(config, nlp_opts, "tol_stat", &tol_stat);
    ocp_nlp_solver_opts_set(config, nlp_opts, "tol_eq", &tol_eq);
    ocp_nlp_solver_opts_set(config, nlp_opts, "tol_ineq", &tol_ineq);
    ocp_nlp_solver_opts_set(config, nlp_opts, "tol_comp", &tol_comp);

	int N2 = N;
    ocp_nlp_solver_opts_set(config, nlp_opts, "qp_cond_N", &N2);
    ocp_nlp_solver_opts_update(config, dims, nlp_opts);

	ocp_nlp_out *nlp_out = ocp_nlp_out_create(config, dims);
	for (int i = 0; i <= N; ++i)
		blasfeo_dvecse(nu[i]+nx[i], 0.0, nlp_out->ux+i, 0);

	// for (int i = 0; i <= N; ++i)
		// BLASFEO_DVECEL(nlp_out->ux+i, 3) = PI;

	ocp_nlp_solver *solver = ocp_nlp_solver_create(config, dims, nlp_opts);
    int solver_status = ocp_nlp_precompute(solver, nlp_in, nlp_out);

	// NLP solution
    acados_timer timer;
    acados_tic(&timer);

	solver_status = ocp_nlp_solve(solver, nlp_in, nlp_out);

    double elapsed_time = acados_toc(&timer);

	printf("\nsolution\n");
	ocp_nlp_out_print(dims, nlp_out);

    printf("\n\nstatus = %i, avg time = %f ms, iters = %d\n\n", solver_status, elapsed_time, nlp_out->sqp_iter);

	return solver_status;
}
