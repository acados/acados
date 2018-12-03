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

#include <acados/utils/print.h>
#include <acados/ocp_qp/ocp_qp_partial_condensing_solver.h>
#include <acados/ocp_nlp/ocp_nlp_constraints_bgh.h>
#include <acados/ocp_nlp/ocp_nlp_cost_ls.h>
#include <acados/ocp_nlp/ocp_nlp_dynamics_cont.h>
#include <acados/ocp_nlp/ocp_nlp_sqp.h>
#include <acados/sim/sim_erk_integrator.h>

#include <acados_c/ocp_nlp_interface.h>

#include <blasfeo/include/blasfeo_d_aux.h>
#include <blasfeo/include/blasfeo_d_aux_ext_dep.h>

#include "{{ ra.model_name }}_model/{{ ra.model_name }}_model.h"
// #include "{{ ra.model_name }}_model/{{ ra.model_name }}_constraint.h"
// {% for item in ra.constants %}
#define {{ item.name }} {{ item.value }}
{% endfor %}
int main() {

    // solver_status = ocp_nlp_solve(solver, nlp_in, nlp_out);
    int num_states = {{ ra.dims.nx }}; 
    int num_controls = {{ ra.dims.nu }}; 
    int N = {{ ra.dims.N }};

    double Tf = 1.0, F_max = 80;
    int idxb_0[5] = {0, 1, 2, 3, 4};
    double b0_l[5] = {-F_max, 0, 0, PI, 0};
    double b0_u[5] = {F_max, 0, 0, PI, 0};
    double Q[{{ ra.dims.nx }}*{{ ra.dims.nx }}]; 
    double R[{{ ra.dims.nu }}*{{ ra.dims.nu }}]; 

    {% for i in range(ra.dims.nx): %}
        {%- for j in range(ra.dims.nx): %}
    Q[{{i}}*{{ra.dims.nx}} + {{j}}] = {{ ra.cost.Q[i,j] }}; 
        {%- endfor %}
    {%- endfor %}

    {% for i in range(ra.dims.nu): %}
        {%- for j in range(ra.dims.nu): %}
    R[{{i}}*{{ra.dims.nu}} + {{j}}] = {{ ra.cost.R[i,j] }}; 
        {%- endfor %}
    {%- endfor %}

    int max_num_sqp_iterations = 1;

    int nx[N+1];
    int nu[N+1];
    int nbx[N+1];
    int nbu[N+1];
    int nb[N+1];
    int ng[N+1];
    int nh[N+1];
    int np[N+1];
    int ns[N+1];
    int nz[N+1];
    int nv[N+1];
    int ny[N+1];

    for(int i = 0; i < N+1; i++) {
        nx[i]  = num_states;
        nu[i]  = num_controls;
        nbx[i] = 0;
        nbu[i] = num_controls;
        nb[i]  = num_controls;
        ng[i]  = 0;
        nh[i]  = 0;
        np[i]  = 0;
        ns[i]  = 0;
        nz[i]  = 0;
        nv[i]  = num_states + num_controls;
        ny[i]  = num_states + num_controls;
    }

    nbx[0] = num_states;
    nbu[0] = num_controls;
    nb[0]  = num_states + num_controls;

    nu[N]  = 0;
    nx[N]  = num_states;
    nh[N]  = 0;
    np[N]  = 0;
    nv[N]  = num_states; 
    ny[N]  = num_states;
    nbu[N] = 0;

    // Make plan
    ocp_nlp_solver_plan *plan = ocp_nlp_plan_create(N);
    plan->nlp_solver = SQP;
    plan->ocp_qp_solver_plan.qp_solver = {{ ra.solver_config.qp_solver }};
    for (int i = 0; i <= N; i++)
        plan->nlp_cost[i] = LINEAR_LS;
    for (int i = 0; i < N; i++)
    {
        plan->nlp_dynamics[i] = CONTINUOUS_MODEL;
        plan->sim_solver_plan[i].sim_solver = {{ ra.solver_config.integrator_type}};

    }

    for (int i = 0; i <= N; i++)
        plan->nlp_constraints[i] = BGH;

    plan->regularization = MIRROR;
    ocp_nlp_solver_config *config = ocp_nlp_config_create(*plan, N);

    ocp_nlp_dims *dims = ocp_nlp_dims_create(config);
    ocp_nlp_dims_initialize(config, nx, nu, ny, nbx, nbu, ng, nh, np, ns, nz, dims);

    {% if ra.solver_config.integrator_type == 'ERK': %}
        // explicit ode
        external_function_casadi forw_vde_casadi[N];
        for (int i = 0; i < N; ++i) {
            forw_vde_casadi[i].casadi_fun = &{{ ra.model_name }}_expl_vde_forw;
            forw_vde_casadi[i].casadi_n_in = &{{ ra.model_name }}_expl_vde_forw_n_in;
            forw_vde_casadi[i].casadi_n_out = &{{ ra.model_name }}_expl_vde_forw_n_out;
            forw_vde_casadi[i].casadi_sparsity_in = &{{ ra.model_name }}_expl_vde_forw_sparsity_in;
            forw_vde_casadi[i].casadi_sparsity_out = &{{ ra.model_name }}_expl_vde_forw_sparsity_out;
            forw_vde_casadi[i].casadi_work = &{{ ra.model_name }}_expl_vde_forw_work;
        }

        {% if ra.solver_config.hessian_approx == 'EXACT': %} 
            external_function_casadi hess_vde_casadi[N];
            for (int i = 0; i < N; ++i) {
                hess_vde_casadi[i].casadi_fun = &{{ ra.model_name }}_expl_ode_hess;
                hess_vde_casadi[i].casadi_n_in = &{{ ra.model_name }}_expl_ode_hess_n_in;
                hess_vde_casadi[i].casadi_n_out = &{{ ra.model_name }}_expl_ode_hess_n_out;
                hess_vde_casadi[i].casadi_sparsity_in = &{{ ra.model_name }}_expl_ode_hess_sparsity_in;
                hess_vde_casadi[i].casadi_sparsity_out = &{{ ra.model_name }}_expl_ode_hess_sparsity_out;
                hess_vde_casadi[i].casadi_work = &{{ ra.model_name }}_expl_ode_hess_work;
            }
        {% endif %}
    {% elif ra.solver_config.integrator_type == 'IRK': %}
        // implicit dae
        external_function_casadi impl_dae_fun[N];
        for (int i = 0; i < N; ++i) {
            impl_dae_fun[i].casadi_fun = &{{ ra.model_name }}_impl_dae_fun;
            impl_dae_fun[i].casadi_work = &{{ ra.model_name }}_impl_dae_fun_work;
            impl_dae_fun[i].casadi_sparsity_in = &{{ ra.model_name }}_impl_dae_fun_sparsity_in;
            impl_dae_fun[i].casadi_sparsity_out = &{{ ra.model_name }}_impl_dae_fun_sparsity_out;
            impl_dae_fun[i].casadi_n_in = &{{ ra.model_name }}_impl_dae_fun_n_in;
            impl_dae_fun[i].casadi_n_out = &{{ ra.model_name }}_impl_dae_fun_n_out;
            external_function_casadi_create(&impl_dae_fun[i]);
        }

        external_function_casadi impl_dae_fun_jac_x_xdot_z[N];
        for (int i = 0; i < N; ++i) {
            impl_dae_fun_jac_x_xdot_z[i].casadi_fun = &{{ ra.model_name }}_impl_dae_fun_jac_x_xdot_z;
            impl_dae_fun_jac_x_xdot_z[i].casadi_work = &{{ ra.model_name }}_impl_dae_fun_jac_x_xdot_z_work;
            impl_dae_fun_jac_x_xdot_z[i].casadi_sparsity_in = &{{ ra.model_name }}_impl_dae_fun_jac_x_xdot_z_sparsity_in;
            impl_dae_fun_jac_x_xdot_z[i].casadi_sparsity_out = &{{ ra.model_name }}_impl_dae_fun_jac_x_xdot_z_sparsity_out;
            impl_dae_fun_jac_x_xdot_z[i].casadi_n_in = &{{ ra.model_name }}_impl_dae_fun_jac_x_xdot_z_n_in;
            impl_dae_fun_jac_x_xdot_z[i].casadi_n_out = &{{ ra.model_name }}_impl_dae_fun_jac_x_xdot_z_n_out;
            external_function_casadi_create(&impl_dae_fun_jac_x_xdot_z[i]);
        }

        external_function_casadi impl_dae_jac_x_xdot_u_z[N];
        for (int i = 0; i < N; ++i) {
            impl_dae_jac_x_xdot_u_z[i].casadi_fun = &{{ ra.model_name }}_impl_dae_jac_x_xdot_u_z;
            impl_dae_jac_x_xdot_u_z[i].casadi_work = &{{ ra.model_name }}_impl_dae_jac_x_xdot_u_z_work;
            impl_dae_jac_x_xdot_u_z[i].casadi_sparsity_in = &{{ ra.model_name }}_impl_dae_jac_x_xdot_u_z_sparsity_in;
            impl_dae_jac_x_xdot_u_z[i].casadi_sparsity_out = &{{ ra.model_name }}_impl_dae_jac_x_xdot_u_z_sparsity_out;
            impl_dae_jac_x_xdot_u_z[i].casadi_n_in = &{{ ra.model_name }}_impl_dae_jac_x_xdot_u_z_n_in;
            impl_dae_jac_x_xdot_u_z[i].casadi_n_out = &{{ ra.model_name }}_impl_dae_jac_x_xdot_u_z_n_out;
            external_function_casadi_create(&impl_dae_jac_x_xdot_u_z[i]);
        }
    {% endif %}


    int function_size;
    char *c_ptr;
    {% if ra.solver_config.integrator_type == 'ERK': %} 
        // NLP model: forward VDEs
        function_size = 0;
        for (int i = 0; i < N; ++i)
            function_size += external_function_casadi_calculate_size(forw_vde_casadi+i);

        c_ptr = (char *) calloc(1, function_size);
        for (int i = 0; i < N; ++i) {
            external_function_casadi_assign(forw_vde_casadi+i, c_ptr);
            c_ptr += external_function_casadi_calculate_size(forw_vde_casadi+i);
        }

        {% if ra.solver_config.hessian_approx == 'EXACT': %} 
        
        // NLP model: hessian VDEs
        function_size = 0;
        for (int i = 0; i < N; ++i)
            function_size += external_function_casadi_calculate_size(hess_vde_casadi+i);

        c_ptr = (char *) calloc(1, function_size);
        for (int i = 0; i < N; ++i) {
            external_function_casadi_assign(hess_vde_casadi+i, c_ptr);
            c_ptr += external_function_casadi_calculate_size(hess_vde_casadi+i);
        }
        {% endif %}
    {% else: %}
        // NLP model: implicit dae fun
        function_size = 0;
        for (int i = 0; i < N; ++i)
            function_size += external_function_casadi_calculate_size(impl_dae_fun+i);

        c_ptr = (char *) calloc(1, function_size);
        for (int i = 0; i < N; ++i) {
            external_function_casadi_assign(impl_dae_fun+i, c_ptr);
            c_ptr += external_function_casadi_calculate_size(impl_dae_fun+i);
        }
        
        // NLP model: implicit dae jacobians
        function_size = 0;
        for (int i = 0; i < N; ++i)
            function_size += external_function_casadi_calculate_size(impl_dae_fun_jac_x_xdot_z+i);

        c_ptr = (char *) calloc(1, function_size);
        for (int i = 0; i < N; ++i) {
            external_function_casadi_assign(impl_dae_fun_jac_x_xdot_z+i, c_ptr);
            c_ptr += external_function_casadi_calculate_size(impl_dae_fun_jac_x_xdot_z+i);
        }

        function_size = 0;
        for (int i = 0; i < N; ++i)
            function_size += external_function_casadi_calculate_size(impl_dae_jac_x_xdot_u_z+i);

        c_ptr = (char *) calloc(1, function_size);
        for (int i = 0; i < N; ++i) {
            external_function_casadi_assign(impl_dae_jac_x_xdot_u_z+i, c_ptr);
            c_ptr += external_function_casadi_calculate_size(impl_dae_jac_x_xdot_u_z+i);
        }
    {% endif %}

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
            for (int k = 0; k < nx[i]; k++)
                BLASFEO_DMATEL(&cost_ls[i]->W, j, k) = Q[j*nx[i] + k];
        for (int j = 0; j < nu[i]; j++)
            for (int k = 0; k < nu[i]; k++)
                BLASFEO_DMATEL(&cost_ls[i]->W, nx[i]+j, nx[i]+k) = R[j*nu[i] + k];
    }
    // WN
    blasfeo_dgese(ny[N], ny[N], 0.0, &cost_ls[N]->W, 0, 0);
    for (int j = 0; j < nx[N]; j++)
        for (int k = 0; k < nx[N]; k++)
            BLASFEO_DMATEL(&cost_ls[N]->W, j, k) = Q[j*nx[N] + k];

    // y_ref
    for (int i = 0; i <= N; ++i)
        blasfeo_dvecse(ny[i], 0.0, &cost_ls[i]->y_ref, 0);

    // NLP dynamics
    int set_fun_status;
    for (int i = 0; i < N; ++i) {
        {% if ra.solver_config.integrator_type == 'ERK': %} 
            set_fun_status = nlp_set_model_in_stage(config, nlp_in, i, "expl_vde_for", &forw_vde_casadi[i]);
            if (set_fun_status != 0) { printf("Error while setting expl_vde_for[%i]\n", i);  exit(1); }
            {% if ra.solver_config.hessian_approx == 'EXACT': %} 
                set_fun_status = nlp_set_model_in_stage(config, nlp_in, i, "expl_ode_hes", &hess_vde_casadi[i]);
                if (set_fun_status != 0) { printf("Error while setting expl_ode_hes[%i]\n", i);  exit(1); }
            {% endif %}
        {% elif ra.solver_config.integrator_type == 'IRK': %} 
			set_fun_status = nlp_set_model_in_stage(config, nlp_in, i, "impl_ode_fun", &impl_dae_fun[i]);
			if (set_fun_status != 0) { printf("Error while setting impl_dae_fun[%i]\n", i);  exit(1); }
			set_fun_status = nlp_set_model_in_stage(config, nlp_in, i, "impl_ode_fun_jac_x_xdot", &impl_dae_fun_jac_x_xdot_z[i]);
			if (set_fun_status != 0) { printf("Error while setting impl_dae_fun_jac_x_xdot_z[%i]\n", i);  exit(1); }
			set_fun_status = nlp_set_model_in_stage(config, nlp_in, i, "impl_ode_jac_x_xdot_u", &impl_dae_jac_x_xdot_u_z[i]);
			if (set_fun_status != 0) { printf("Error while setting impl_dae_jac_x_xdot_u_z[%i]\n", i);  exit(1); }
        {% endif %}
    }

    // NLP constraints
    ocp_nlp_constraints_bgh_model **constraints = (ocp_nlp_constraints_bgh_model **) nlp_in->constraints;

    // bounds
    constraints[0]->idxb = idxb_0;
    blasfeo_pack_dvec(nb[0], b0_l, &constraints[0]->d, 0);
    blasfeo_pack_dvec(nb[0], b0_u, &constraints[0]->d, nb[0]+ng[0]);

    double lbu[1] = {-F_max}, ubu[1] = {F_max};
    for (int i = 1; i < N; ++i)
    {
        constraints[i]->idxb[0] = 0;
        blasfeo_pack_dvec(nb[i], lbu, &constraints[i]->d, 0);
        blasfeo_pack_dvec(nb[i], ubu, &constraints[i]->d, nb[i]+ng[i]);
    }

    void *nlp_opts = ocp_nlp_opts_create(config, dims);

    ocp_nlp_sqp_opts *sqp_opts = (ocp_nlp_sqp_opts *) nlp_opts;
    sqp_opts->maxIter = max_num_sqp_iterations;
    sqp_opts->min_res_g = 1e-9;
    sqp_opts->min_res_b = 1e-9;
    sqp_opts->min_res_d = 1e-9;
    sqp_opts->min_res_m = 1e-9;
    for (int i = 0; i < N; ++i)
    {
        sim_rk_opts *rk_opts = (sim_rk_opts *) ((ocp_nlp_dynamics_cont_opts *)sqp_opts->dynamics[i])->sim_solver;
        rk_opts->num_steps = 5;
        {% if ra.solver_config.hessian_approx == 'EXACT': %} 
        ((ocp_nlp_dynamics_cont_opts *)sqp_opts->dynamics[i])->compute_hess = true;
        rk_opts->sens_hess = true;
        rk_opts->sens_adj = true;
        {% endif %}
    }

    ocp_nlp_out *nlp_out = ocp_nlp_out_create(config, dims);
    for (int i = 0; i <= N; ++i)
        blasfeo_dvecse(nu[i]+nx[i], 0.0, nlp_out->ux+i, 0);

    ocp_nlp_solver *solver = ocp_nlp_create(config, dims, nlp_opts);

    // NLP solution
    acados_timer timer;
    double kkt_norm_inf = 1e12, elapsed_time;
    
    int solver_status = 0, iteration_number = 0;

    while (kkt_norm_inf > 1e-9) {
        acados_tic(&timer);
        solver_status = ocp_nlp_solve(solver, nlp_in, nlp_out);
        elapsed_time = acados_toc(&timer);
        kkt_norm_inf = nlp_out->inf_norm_res;
        printf(" iteration %2d | time  %f |  KKT %e\n", iteration_number, elapsed_time, kkt_norm_inf);
        iteration_number++;

        if (iteration_number >= 100)
            break;
    }

    printf("\n--- solution ---\n");
    ocp_nlp_out_print(dims, nlp_out);
    
    // // free memory
    // free(dims);
    // free(config);
    // free(nlp_in);
    // free(nlp_out);
    // free(nlp_opts);
    // free(solver);
    
    // // free external function 
    // external_function_casadi_free(&impl_dae_fun);
    // external_function_casadi_free(&impl_dae_fun_jac_x_xdot_z);
    // external_function_casadi_free(&impl_dae_jac_x_xdot_u_z);
    // external_function_casadi_free(&nls_cost_residual);
    // external_function_casadi_free(&nls_cost_N_residual);

    return solver_status;
}
