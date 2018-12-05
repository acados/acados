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

#include <numeric>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "acados/utils/print.h"
#include "acados/ocp_qp/ocp_qp_partial_condensing_solver.h"
#include "acados/ocp_nlp/ocp_nlp_constraints_bgh.h"
#include "acados/ocp_nlp/ocp_nlp_cost_ls.h"
#include "acados/ocp_nlp/ocp_nlp_dynamics_cont.h"
#include "acados/ocp_nlp/ocp_nlp_sqp.h"
#include "acados/sim/sim_erk_integrator.h"

#include "acados_c/ocp_nlp_interface.h"

#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"

#include "pendulum_model/pendulum_model.h"
// #include "pendulum_model/jac_constraint.h"
#include "pendulum_model/constraint.h"

std::ostream& operator<<(std::ostream& strm, std::vector<double> v) {
	for (auto e : v)
		strm << std::scientific << std::showpos << e << " ";
	strm << "\n";
	return strm;
}

int main() {

	int num_states = 4, num_controls = 1, N = 100;
	double Tf = 1.0, R = 1e-2, F_max = 80;
	std::vector<int> idxb_0 {0, 1, 2, 3, 4};
	std::vector<double> b0_l {-F_max, 0, 0, M_PI, 0}, b0_u {F_max, 0, 0, M_PI, 0};
	std::vector<double> Q {1e3, 1e-2, 1e3, 1e-2};

	int max_num_sqp_iterations = 1;

	std::vector<int> nx(N+1, num_states), nu(N+1, num_controls), nbx(N+1, 0), nbu(N+1, 0),
		nb(N+1, 0), ng(N+1, 0), nh(N+1, 0), nq(N+1, 0),
		ns(N+1, 0), nv(N+1, num_states+num_controls), ny(N+1, num_states+num_controls), nz(N+1, 0);

	nbx.at(0) = num_states;
	nbu.at(0) = num_controls;
	nb.at(0) = num_states + num_controls;

	for (int i = 1; i < N; ++i)
	{
		nbu.at(i) = num_controls;
		nb.at(i) = num_controls;
	}

	nu.at(N) = 0;
	nv.at(N) = 4;
	ny.at(N) = 4;

	// Make plan

	ocp_qp_solver_plan qp_plan = {FULL_CONDENSING_QPOASES};
	std::vector<sim_solver_plan> sim_plan(N, {ERK});
	std::vector<ocp_nlp_cost_t> cost_plan(N+1, LINEAR_LS);
	std::vector<ocp_nlp_dynamics_t> dynamics_plan(N, CONTINUOUS_MODEL);
	std::vector<ocp_nlp_constraints_t> constraints_plan(N+1, BGH);

	ocp_nlp_solver_plan plan = {
		qp_plan,
		sim_plan.data(),
		SQP,
		MIRROR,
		cost_plan.data(),
		dynamics_plan.data(),
		constraints_plan.data()
	};
	ocp_nlp_solver_config *config = ocp_nlp_config_create(plan);

	ocp_nlp_dims *dims = ocp_nlp_dims_create(config);
	ocp_nlp_dims_initialize(config, nx.data(), nu.data(), ny.data(), nbx.data(), nbu.data(), ng.data(), nh.data(), nq.data(), ns.data(), nz.data(), dims);

	external_function_casadi forw_vde_casadi[N];
	for (int i = 0; i < N; ++i) {
		forw_vde_casadi[i].casadi_fun = &pendulum_ode_expl_vde_forw;
		forw_vde_casadi[i].casadi_n_in = &pendulum_ode_expl_vde_forw_n_in;
		forw_vde_casadi[i].casadi_n_out = &pendulum_ode_expl_vde_forw_n_out;
		forw_vde_casadi[i].casadi_sparsity_in = &pendulum_ode_expl_vde_forw_sparsity_in;
		forw_vde_casadi[i].casadi_sparsity_out = &pendulum_ode_expl_vde_forw_sparsity_out;
		forw_vde_casadi[i].casadi_work = &pendulum_ode_expl_vde_forw_work;
	}

	external_function_casadi hess_vde_casadi[N];
	for (int i = 0; i < N; ++i) {
		hess_vde_casadi[i].casadi_fun = &pendulum_ode_expl_ode_hess;
		hess_vde_casadi[i].casadi_n_in = &pendulum_ode_expl_ode_hess_n_in;
		hess_vde_casadi[i].casadi_n_out = &pendulum_ode_expl_ode_hess_n_out;
		hess_vde_casadi[i].casadi_sparsity_in = &pendulum_ode_expl_ode_hess_sparsity_in;
		hess_vde_casadi[i].casadi_sparsity_out = &pendulum_ode_expl_ode_hess_sparsity_out;
		hess_vde_casadi[i].casadi_work = &pendulum_ode_expl_ode_hess_work;
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

	// NLP model: hessian VDEs
	function_size = 0;
	for (int i = 0; i < N; ++i)
		function_size += external_function_casadi_calculate_size(hess_vde_casadi+i);

	c_ptr = (char *) calloc(1, function_size);
	for (int i = 0; i < N; ++i) {
		external_function_casadi_assign(hess_vde_casadi+i, c_ptr);
		c_ptr += external_function_casadi_calculate_size(hess_vde_casadi+i);
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
            BLASFEO_DMATEL(&cost_ls[i]->W, j, j) = Q[j];
        for (int j = 0; j < nu[i]; j++)
            BLASFEO_DMATEL(&cost_ls[i]->W, nx[i]+j, nx[i]+j) = R;
	}
	// WN
	blasfeo_dgese(ny[N], ny[N], 0.0, &cost_ls[N]->W, 0, 0);
	for (int j = 0; j < nx[N]; j++)
		BLASFEO_DMATEL(&cost_ls[N]->W, j, j) = Q[j];

	// y_ref
    for (int i = 0; i <= N; ++i)
		blasfeo_dvecse(ny[i], 0.0, &cost_ls[i]->y_ref, 0);

	// NLP dynamics
	for (int i = 0; i < N; ++i) {
		ocp_nlp_dynamics_cont_model *dynamics = (ocp_nlp_dynamics_cont_model *) nlp_in->dynamics[i];
		erk_model *model = (erk_model *) dynamics->sim_model;
		model->expl_vde_for = (external_function_generic *) &forw_vde_casadi[i];
		model->expl_vde_adj = (external_function_generic *) &hess_vde_casadi[i];
		model->expl_ode_hes = (external_function_generic *) &hess_vde_casadi[i];
	}

	// NLP constraints
	ocp_nlp_constraints_bgh_model **constraints = (ocp_nlp_constraints_bgh_model **) nlp_in->constraints;

	// external_function_casadi nonlinear_constraint;
	// nonlinear_constraint.casadi_fun = &constraint;
	// nonlinear_constraint.casadi_n_in = &constraint_n_in;
	// nonlinear_constraint.casadi_n_out = &constraint_n_out;
	// nonlinear_constraint.casadi_sparsity_in = &constraint_sparsity_in;
	// nonlinear_constraint.casadi_sparsity_out = &constraint_sparsity_out;
	// nonlinear_constraint.casadi_work = &constraint_work;

	// int constraint_size = external_function_casadi_calculate_size(&nonlinear_constraint);
	// void *ptr = malloc(constraint_size);
	// external_function_casadi_assign(&nonlinear_constraint, ptr);

	// bounds
    constraints[0]->idxb = idxb_0.data();
	blasfeo_pack_dvec(nb[0], b0_l.data(), &constraints[0]->d, 0);
	blasfeo_pack_dvec(nb[0], b0_u.data(), &constraints[0]->d, nb[0]+ng[0]);

	double lbu[1] = {-F_max}, ubu[1] = {F_max};
	for (int i = 1; i < N; ++i)
	{
		constraints[i]->idxb[0] = 0;
		blasfeo_pack_dvec(nb[i], lbu, &constraints[i]->d, 0);
		blasfeo_pack_dvec(nb[i], ubu, &constraints[i]->d, nb[i]+ng[i]);
	}

    // constraints[N]->idxb = idxb_N.data();
	// blasfeo_pack_dvec(nb[N], xN.data(), &constraints[N]->d, 0);
	// blasfeo_pack_dvec(nb[N], xN.data(), &constraints[N]->d, nb[N]+ng[N]+nh[N]);

	void *nlp_opts = ocp_nlp_opts_create(config, dims);


    int maxIter = max_num_sqp_iterations;
    double min_res_g = 1e-9;
    double min_res_b = 1e-9;
    double min_res_d = 1e-9;
    double min_res_m = 1e-9;

	ocp_nlp_opts_set(config, nlp_opts, "maxIter", &maxIter);
	ocp_nlp_opts_set(config, nlp_opts, "min_res_g", &min_res_g);
	ocp_nlp_opts_set(config, nlp_opts, "min_res_b", &min_res_b);
	ocp_nlp_opts_set(config, nlp_opts, "min_res_d", &min_res_d);
	ocp_nlp_opts_set(config, nlp_opts, "min_res_m", &min_res_m);

	for (int i = 0; i < N; ++i)
	{
		sim_rk_opts *rk_opts = (sim_rk_opts *) ((ocp_nlp_dynamics_cont_opts *)sqp_opts->dynamics[i])->sim_solver;
		((ocp_nlp_dynamics_cont_opts *)sqp_opts->dynamics[i])->compute_hess = true;
		rk_opts->num_steps = 5;
		rk_opts->sens_hess = true;
		rk_opts->sens_adj = true;
	}
	// ((ocp_qp_partial_condensing_solver_opts *) sqp_opts->qp_solver_opts)->pcond_opts->N2 = 10;

	ocp_nlp_out *nlp_out = ocp_nlp_out_create(config, dims);
	for (int i = 0; i <= N; ++i)
		blasfeo_dvecse(nu[i]+nx[i], 0.0, nlp_out->ux+i, 0);

	// for (int i = 0; i <= N; ++i)
	// 	BLASFEO_DVECEL(nlp_out->ux+i, 3) = M_PI;

	ocp_nlp_solver *solver = ocp_nlp_create(config, dims, nlp_opts);
	solver_status = ocp_nlp_precompute(solver, nlp_in, nlp_out);

	std::vector<std::vector<double>> MPC_log;
	std::vector<double> timings, KKT_log;

	// NLP solution
    acados_timer timer;
	double kkt_norm_inf = __builtin_inf(), elapsed_time;
	
	int solver_status = 0, iteration_number = 0;
	// for (int j = 0; j < 100; ++j) {
		// kkt_norm_inf = __builtin_inf(), elapsed_time = 0;
		// solver_status = 0, iteration_number = 0;
		// for (int i = 0; i <= N; ++i)
			// blasfeo_dvecse(nu[i]+nx[i], 0.0, nlp_out->ux+i, 0);

		// for (int i = 0; i <= N; ++i)
			// BLASFEO_DVECEL(nlp_out->ux+i, 3) = M_PI;
		// blasfeo_pack_dvec(nb[0], x0.data(), &constraints[0]->d, 0);
		// blasfeo_pack_dvec(nb[0], x0.data(), &constraints[0]->d, nb[0]+ng[0]);

		while (kkt_norm_inf > 1e-9) {
			acados_tic(&timer);
			solver_status = ocp_nlp_solve(solver, nlp_in, nlp_out);
			elapsed_time = acados_toc(&timer);
			timings.push_back(elapsed_time);
			kkt_norm_inf = nlp_out->inf_norm_res;
			KKT_log.push_back(kkt_norm_inf);
			printf(" iteration %2d | time  %f |  KKT %e\n", iteration_number, elapsed_time, kkt_norm_inf);
			// blasfeo_print_tran_dvec(nx.at(0)+nu.at(0), &nlp_out->ux[0], 0);
			std::vector<double> ux0(nx[0]+nu[0]);
			for (int i = 0; i < N; ++i)
			{
				blasfeo_unpack_dvec(nu[i]+nx[i], &nlp_out->ux[i], 0, ux0.data());
				MPC_log.push_back(ux0);
			}
			ux0.front() = 0;
			blasfeo_unpack_dvec(nu[N]+nx[N], &nlp_out->ux[N], 0, ux0.data()+1);
			MPC_log.push_back(ux0);

			// blasfeo_dveccp(nb[0], &nlp_out->ux[1], nu.at(0), &constraints[0]->d, 0);
			// blasfeo_dveccp(nb[0], &nlp_out->ux[1], nu.at(0), &constraints[0]->d, nb[0]+ng[0]);
			iteration_number++;

			if (iteration_number >= 100)
				break;
		}
	// }

	std::ofstream ostrm("states_controls.txt");
	for (auto v : MPC_log)
		ostrm << v;
	
	std::ofstream ostrm_timings("timings.txt"), ostrm_KKT("KKT.txt");
	ostrm_timings << timings;
	ostrm_KKT << KKT_log;

	printf("\n--- solution ---\n");
	ocp_nlp_out_print(dims, nlp_out);

    printf("\n\nstatus = %i, avg time = %f ms, iters = %d\n\n", solver_status, std::accumulate(std::begin(timings), std::end(timings), 0.0)/iteration_number, nlp_out->sqp_iter);

	return solver_status;
}
