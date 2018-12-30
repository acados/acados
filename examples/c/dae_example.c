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

#include <stdio.h>

#include "acados/utils/print.h"
#include "acados/ocp_qp/ocp_qp_partial_condensing_solver.h"
#include "acados/ocp_nlp/ocp_nlp_constraints_bghp.h"
#include "acados/ocp_nlp/ocp_nlp_cost_ls.h"
#include "acados/ocp_nlp/ocp_nlp_dynamics_common.h"
#include "acados/ocp_nlp/ocp_nlp_dynamics_cont.h"
#include "acados/ocp_nlp/ocp_nlp_sqp.h"
#include "acados/sim/sim_common.h"
#include "acados/sim/sim_irk_integrator.h"

#include "acados_c/ocp_nlp_interface.h"

#include "blasfeo/include/blasfeo_d_aux.h"

#include "inverted_pendulum_model/inverted_pendulum_model_dae.h"

#define PI 3.1415926535897932

int main() {

	int num_states = 6, num_controls = 1, N = 10;
	double Tf = 1.0, R = 1e-2, QN = 1e1, F_max = 80;
    double Q[6] = {1e-3, 1e-3, 1e-3, 1e-3, 1e2, 1e-2};
	int idxb_0[6] = {1, 2, 3, 4, 5, 6};
	double x0[num_states];

    x0[0] =  1.0000;  // xpos
    x0[1] = -5.0000;  // ypos
    x0[2] =  0.1000;  // vx
    x0[3] = -0.5000;  // vy
    x0[4] =  0.1000;  // valpha
    x0[5] =  1.0000;  // alpha

	double radius2 = 0.04, neg_inf = -1000000;
	int max_num_sqp_iterations = 100;

	int nx[N+1];
    int nu[N+1];
    int nbx[N+1];
    int nbu[N+1];
	int	nb[N+1];
    int ng[N+1];
    int nh[N+1];
    int np[N+1];
	int	ns[N+1];
    int nz[N+1];
    int nv[N+1];
    int ny[N+1];

    for(int i = 0; i < N+1; i++) {
        nx[i] = num_states;
        nu[i] = num_controls;
        nbx[i] = 0;
        nbu[i] = 0;
        nb[i] = 0;
        ng[i] = 0;
        nh[i] = 0;
        np[i] = 0;
        ns[i] = 0;
        nz[i] = 5;
        nv[i] = num_states + num_controls;
        ny[i] = num_states + num_controls;
    }

    nbx[0] = num_states;
    nb[0] = num_states;
    nu[N] = 0;
    nh[N] = 0;
    np[N] = 0;
    nv[N] = num_states; 
    ny[N] = num_states;

	// Make plan
	ocp_nlp_solver_plan *plan = ocp_nlp_plan_create(N);
	plan->nlp_solver = SQP;
	// plan->ocp_qp_solver_plan.qp_solver = PARTIAL_CONDENSING_HPIPM;
	plan->ocp_qp_solver_plan.qp_solver = FULL_CONDENSING_QPOASES;
	for (int i = 0; i <= N; i++)
		plan->nlp_cost[i] = LINEAR_LS;
	for (int i = 0; i < N; i++)
	{
		plan->nlp_dynamics[i] = CONTINUOUS_MODEL;
		plan->sim_solver_plan[i].sim_solver = IRK;
	}

	for (int i = 0; i <= N; i++)
		plan->nlp_constraints[i] = BGH;

	ocp_nlp_solver_config *config = ocp_nlp_config_create(*plan);

	ocp_nlp_dims *dims = ocp_nlp_dims_create(config);
	// ocp_nlp_dims_initialize(config, nx, nu, ny, nbx, nbu, ng, nh, np, ns, nz, dims);

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
    }

	external_function_casadi impl_ode_fun[N];
	external_function_casadi impl_ode_fun_jac_x_xdot_z[N];
	// external_function_casadi impl_ode_fun_jac_x_xdot_z_u[N];
	external_function_casadi impl_ode_jac_x_xdot_z_u[N];

	for (int ii = 0; ii < N; ++ii) {
        impl_ode_fun[ii].casadi_fun = &casadi_impl_ode_fun_pendulum_dae;
        impl_ode_fun[ii].casadi_work = &casadi_impl_ode_fun_pendulum_dae_work;
        impl_ode_fun[ii].casadi_sparsity_in = &casadi_impl_ode_fun_pendulum_dae_sparsity_in;
        impl_ode_fun[ii].casadi_sparsity_out = &casadi_impl_ode_fun_pendulum_dae_sparsity_out;
        impl_ode_fun[ii].casadi_n_in = &casadi_impl_ode_fun_pendulum_dae_n_in;
        impl_ode_fun[ii].casadi_n_out = &casadi_impl_ode_fun_pendulum_dae_n_out;

        impl_ode_fun_jac_x_xdot_z[ii].casadi_fun = &casadi_impl_ode_fun_jac_x_xdot_z_pendulum_dae;
        impl_ode_fun_jac_x_xdot_z[ii].casadi_work = &casadi_impl_ode_fun_jac_x_xdot_z_pendulum_dae_work;
        impl_ode_fun_jac_x_xdot_z[ii].casadi_sparsity_in = &casadi_impl_ode_fun_jac_x_xdot_z_pendulum_dae_sparsity_in;
        impl_ode_fun_jac_x_xdot_z[ii].casadi_sparsity_out = &casadi_impl_ode_fun_jac_x_xdot_z_pendulum_dae_sparsity_out;
        impl_ode_fun_jac_x_xdot_z[ii].casadi_n_in = &casadi_impl_ode_fun_jac_x_xdot_z_pendulum_dae_n_in;
        impl_ode_fun_jac_x_xdot_z[ii].casadi_n_out = &casadi_impl_ode_fun_jac_x_xdot_z_pendulum_dae_n_out;

        // impl_ode_fun_jac_x_xdot_z_u[ii].casadi_fun = &casadi_impl_ode_fun_jac_x_xdot_z_u_pendulum_dae;
        // impl_ode_fun_jac_x_xdot_z_u[ii].casadi_work = &casadi_impl_ode_fun_jac_x_xdot_z_u_pendulum_dae_work;
        // impl_ode_fun_jac_x_xdot_z_u[ii].casadi_sparsity_in = &casadi_impl_ode_fun_jac_x_xdot_z_u_pendulum_dae_sparsity_in;
        // impl_ode_fun_jac_x_xdot_z_u[ii].casadi_sparsity_out = &casadi_impl_ode_fun_jac_x_xdot_z_u_pendulum_dae_sparsity_out;
        // impl_ode_fun_jac_x_xdot_z_u[ii].casadi_n_in = &casadi_impl_ode_fun_jac_x_xdot_z_u_pendulum_dae_n_in;
        // impl_ode_fun_jac_x_xdot_z_u[ii].casadi_n_out = &casadi_impl_ode_fun_jac_x_xdot_z_u_pendulum_dae_n_out;

        impl_ode_jac_x_xdot_z_u[ii].casadi_fun = &casadi_impl_ode_jac_x_xdot_z_u_pendulum_dae;
        impl_ode_jac_x_xdot_z_u[ii].casadi_work = &casadi_impl_ode_jac_x_xdot_z_u_pendulum_dae_work;
        impl_ode_jac_x_xdot_z_u[ii].casadi_sparsity_in = &casadi_impl_ode_jac_x_xdot_z_u_pendulum_dae_sparsity_in;
        impl_ode_jac_x_xdot_z_u[ii].casadi_sparsity_out = &casadi_impl_ode_jac_x_xdot_z_u_pendulum_dae_sparsity_out;
        impl_ode_jac_x_xdot_z_u[ii].casadi_n_in = &casadi_impl_ode_jac_x_xdot_z_u_pendulum_dae_n_in;
        impl_ode_jac_x_xdot_z_u[ii].casadi_n_out = &casadi_impl_ode_jac_x_xdot_z_u_pendulum_dae_n_out;
	}

	// impl_ode
    int	tmp_size = 0;
	for (int ii=0; ii<N; ii++)
	{
		tmp_size += external_function_casadi_calculate_size(impl_ode_fun+ii);
	}
	void *impl_ode_casadi_mem = malloc(tmp_size);
	void *c_ptr = impl_ode_casadi_mem;
	for (int ii=0; ii<N; ii++)
	{
		external_function_casadi_assign(impl_ode_fun+ii, c_ptr);
		c_ptr += external_function_casadi_calculate_size(impl_ode_fun+ii);
	}
	//
	tmp_size = 0;
	for (int ii=0; ii<N; ii++)
	{
		tmp_size += external_function_casadi_calculate_size(impl_ode_fun_jac_x_xdot_z+ii);
	}
	void *impl_ode_fun_jac_x_xdot_z_mem = malloc(tmp_size);
	c_ptr = impl_ode_fun_jac_x_xdot_z_mem;
	for (int ii=0; ii<N; ii++)
	{
		external_function_casadi_assign(impl_ode_fun_jac_x_xdot_z+ii, c_ptr);
		c_ptr += external_function_casadi_calculate_size(impl_ode_fun_jac_x_xdot_z+ii);
	}
	//
	// tmp_size = 0;
	// for (int ii=0; ii<N; ii++)
	// {
	// 	tmp_size += external_function_casadi_calculate_size(impl_ode_fun_jac_x_xdot_z_u+ii);
	// }
	// void *impl_ode_fun_jac_x_xdot_z_u_mem = malloc(tmp_size);
	// c_ptr = impl_ode_fun_jac_x_xdot_z_u_mem;
	// for (int ii=0; ii<N; ii++)
	// {
	// 	external_function_casadi_assign(impl_ode_fun_jac_x_xdot_z_u+ii, c_ptr);
	// 	c_ptr += external_function_casadi_calculate_size(impl_ode_fun_jac_x_xdot_z_u+ii);
	// }
	//
	tmp_size = 0;
	for (int ii=0; ii<N; ii++)
	{
		tmp_size += external_function_casadi_calculate_size(impl_ode_jac_x_xdot_z_u+ii);
	}
	void *impl_ode_jac_x_xdot_z_u_mem = malloc(tmp_size);
	c_ptr = impl_ode_jac_x_xdot_z_u_mem;
	for (int ii=0; ii<N; ii++)
	{
		external_function_casadi_assign(impl_ode_jac_x_xdot_z_u+ii, c_ptr);
		c_ptr += external_function_casadi_calculate_size(impl_ode_jac_x_xdot_z_u+ii);
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
		BLASFEO_DMATEL(&cost_ls[N]->W, j, j) = QN;

	// y_ref
    for (int i = 0; i <= N; ++i)
		blasfeo_dvecse(ny[i], 0.0, &cost_ls[i]->y_ref, 0);

	// NLP dynamics
	for (int i = 0; i < N; ++i) {
		ocp_nlp_dynamics_cont_model *dynamics = nlp_in->dynamics[i];
		irk_model *model = dynamics->sim_model;
		model->impl_ode_fun = (external_function_generic *) &impl_ode_fun[i];
		model->impl_ode_fun_jac_x_xdot_z = (external_function_generic *) &impl_ode_fun_jac_x_xdot_z[i];
		model->impl_ode_jac_x_xdot_z_u = (external_function_generic *) &impl_ode_jac_x_xdot_z_u[i]; // TODO(zanellia): need to swap z and u!!
	}

	// bounds
	ocp_nlp_constraints_bgh_model **constraints = (ocp_nlp_constraints_bgh_model **) nlp_in->constraints;

    constraints[0]->idxb = idxb_0;
	blasfeo_pack_dvec(nb[0], x0, &constraints[0]->d, 0);
	blasfeo_pack_dvec(nb[0], x0, &constraints[0]->d, nb[0]+ng[0]);

	void *nlp_opts = ocp_nlp_opts_create(config, dims);

	ocp_nlp_sqp_opts *sqp_opts = (ocp_nlp_sqp_opts *) nlp_opts;
    sqp_opts->maxIter = max_num_sqp_iterations;
    sqp_opts->min_res_g = 1e-9;
    sqp_opts->min_res_b = 1e-9;
    sqp_opts->min_res_d = 1e-9;
    sqp_opts->min_res_m = 1e-9;
	((ocp_qp_partial_condensing_solver_opts *) sqp_opts->qp_solver_opts)->pcond_opts->N2 = N;

	ocp_nlp_out *nlp_out = ocp_nlp_out_create(config, dims);
	for (int i = 0; i <= N; ++i) {
		blasfeo_dvecse(nu[i]+nx[i], 0.0, nlp_out->ux+i, 0);
		blasfeo_dvecse(1, 0.1, nlp_out->ux+i, 0);
    }

	ocp_nlp_solver *solver = ocp_nlp_create(config, dims, nlp_opts);

	// NLP solution
    acados_timer timer;
    acados_tic(&timer);

	int solver_status = ocp_nlp_solve(solver, nlp_in, nlp_out);

    double elapsed_time = acados_toc(&timer);

	printf("\nsolution\n");
	ocp_nlp_out_print(dims, nlp_out);

    printf("\n\nstatus = %i, avg time = %f ms, iters = %d\n\n", solver_status, elapsed_time, nlp_out->sqp_iter);

	return solver_status;
}
