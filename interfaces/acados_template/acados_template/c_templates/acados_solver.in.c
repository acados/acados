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

// standard
#include <stdio.h>
#include <stdlib.h>
// acados
#include "acados/utils/print.h"
#include "acados_c/ocp_nlp_interface.h"
#include "acados_c/external_function_interface.h"

// TODO(oj): remove, when setters for Cyt,idxb available
#include "acados/ocp_nlp/ocp_nlp_constraints_bgh.h"
#include "acados/ocp_nlp/ocp_nlp_cost_ls.h"

// blasfeo
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"

// example specific
#include "{{ ra.model_name }}_model/{{ ra.model_name }}_model.h"
// #include "{{ ra.model_name }}_model/{{ ra.model_name }}_constraint.h"

#include "acados_solver_{{ra.model_name}}.h"

{% for item in ra.constants %}
#define {{ item.name }} {{ item.value }}
{% endfor %}

int acados_create() {

    // ocp_nlp_solver * nlp_solver;
    // ocp_nlp_in * nlp_in;
    // ocp_nlp_out * nlp_out;
    // void *nlp_opts;
    // ocp_nlp_config *config;
    // ocp_nlp_plan *plan;
    // ocp_nlp_dims *dims;

    int status = 0;

    int num_states = {{ ra.dims.nx }}; 
    int num_controls = {{ ra.dims.nu }}; 
    int N = {{ ra.dims.N }};

    double Tf = {{ ra.solver_config.tf }};

    int ny_ = num_controls + num_states;
    // set up bounds for stage 0 
    int idxb_0[{{ ra.dims.nbu }} + {{ ra.dims.nx }}];
    {% for i in range(ra.dims.nbu + ra.dims.nx): %}
    idxb_0[{{i}}] = {{i}};
    {%- endfor %}
    double lb0[{{ ra.dims.nbu }} + {{ ra.dims.nx }}]; 
    double ub0[{{ ra.dims.nbu }} + {{ ra.dims.nx }}];
    {% for i in range(ra.dims.nbu): %}
    lb0[{{i}}] = {{ra.constraints.lbu[i]}};
    ub0[{{i}}] = {{ra.constraints.ubu[i]}};
    {%- endfor %}

    {% for i in range(ra.dims.nbu, ra.dims.nx + ra.dims.nbu): %}
    lb0[{{i}}] = {{ra.constraints.x0[i - ra.dims.nbu]}};
    ub0[{{i}}] = {{ra.constraints.x0[i - ra.dims.nbu]}};
    {%- endfor %}

    // set up bounds for intermediate stages
    int idxb[{{ ra.dims.nbu }} + {{ ra.dims.nx }}];
    {% for i in range(ra.dims.nbu + ra.dims.nbx): %}
    idxb[{{i}}] = {{i}};
    {%- endfor %}
    double lb[{{ ra.dims.nbu }} + {{ ra.dims.nbx }}]; 
    double ub[{{ ra.dims.nbu }} + {{ ra.dims.nbx }}]; 
    {% for i in range(ra.dims.nbu): %}
    lb[{{i}}] = {{ra.constraints.lbu[i]}};
    ub[{{i}}] = {{ra.constraints.ubu[i]}};
    {%- endfor %}
    {% for i in range(ra.dims.nbu, ra.dims.nbx + ra.dims.nbu): %}
    lb[{{i}}] = {{ra.constraints.lbx[i]}};
    ub[{{i}}] = {{ra.constraints.ubx[i]}};
    {%- endfor %}

    // set up bounds for last stage
    int idxb_N[{{ ra.dims.nbx }}];
    {% for i in range(ra.dims.nbx): %}
    idxb_N[{{i}}] = {{i}};
    {%- endfor %}
    double lbN[{{ ra.dims.nbx }}]; 
    double ubN[{{ ra.dims.nbx }}]; 
    {% for i in range(ra.dims.nbx): %}
    lbN[{{i}}] = {{ra.constraints.lbx[i]}};
    ubN[{{i}}] = {{ra.constraints.ubx[i]}};
    {%- endfor %}

    double yref[{{ ra.dims.nx }}+{{ ra.dims.nu }}];
    double Q[{{ ra.dims.nx }}*{{ ra.dims.nx }}]; 
    double R[{{ ra.dims.nu }}*{{ ra.dims.nu }}];
    double W[(num_states + num_controls)*(num_states + num_controls)];

    for (int ii = 0; ii < num_controls + num_states; ii++)
        yref[ii] = 0.0;

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

    for (int ii = 0; ii < ny_ * ny_; ii++)
        W[ii] = 0.0;

    {% for j in range(ra.dims.nx): %}
        {%- for k in range(ra.dims.nx): %}
    W[{{j}}+({{ra.dims.nx}} + {{ra.dims.nu}}) * {{k}}] = {{ ra.cost.Q[j,k] }}; 
        {%- endfor %}
    {%- endfor %}

    {% for j in range(ra.dims.nx,ra.dims.nx+ra.dims.nu): %}
        {%- for k in range(ra.dims.nx,ra.dims.nx+ra.dims.nu): %}
    W[{{j}}+({{ra.dims.nx}} + {{ra.dims.nu}}) * {{k}}] = {{ ra.cost.R[j-ra.dims.nx,k-ra.dims.nx] }}; 
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
        nbx[i] = {{ra.dims.nbx}};
        nbu[i] = {{ra.dims.nbu}};
        nb[i]  = {{ra.dims.nbu}} + {{ra.dims.nbx}};
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
    nbx[N] = {{ra.dims.nbx}};
    nb[N]  = {{ra.dims.nbx}};

    // Make plan
    nlp_solver_plan = ocp_nlp_plan_create(N);
    {% if ra.solver_config.nlp_solver_type == 'SQP': %}
    nlp_solver_plan->nlp_solver = SQP;
    {% else: %}
    nlp_solver_plan->nlp_solver = SQP_RTI;
    {% endif %}
    nlp_solver_plan->ocp_qp_solver_plan.qp_solver = {{ ra.solver_config.qp_solver }};
    for (int i = 0; i <= N; i++)
        nlp_solver_plan->nlp_cost[i] = LINEAR_LS;
    for (int i = 0; i < N; i++)
    {
        nlp_solver_plan->nlp_dynamics[i] = CONTINUOUS_MODEL;
        nlp_solver_plan->sim_solver_plan[i].sim_solver = {{ ra.solver_config.integrator_type}};
    }

    for (int i = 0; i <= N; i++)
        nlp_solver_plan->nlp_constraints[i] = BGH;

    {% if ra.solver_config.hessian_approx == 'EXACT': %} 
    nlp_solver_plan->regularization = CONVEXIFICATION;
    {% endif %}
    nlp_config = ocp_nlp_config_create(*nlp_solver_plan);

    /* create and set ocp_nlp_dims */
	nlp_dims = ocp_nlp_dims_create(nlp_config);

    ocp_nlp_dims_set_opt_vars(nlp_config, nlp_dims, "nx", nx);
    ocp_nlp_dims_set_opt_vars(nlp_config, nlp_dims, "nu", nu);
    ocp_nlp_dims_set_opt_vars(nlp_config, nlp_dims, "nz", nz);
    ocp_nlp_dims_set_opt_vars(nlp_config, nlp_dims, "ns", ns);

    ocp_nlp_dims_set_cost(nlp_config, nlp_dims, "ny", ny);

    ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, "nbx", nbx);
    ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, "nbu", nbu);
    ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, "ng", ng);
    ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, "nh", nh);

    {% if ra.solver_config.integrator_type == 'ERK': %}
    // explicit ode
    external_function_casadi * forw_vde_casadi;
    forw_vde_casadi = (external_function_casadi *) malloc(sizeof(external_function_casadi)*N);
    for (int i = 0; i < N; ++i) {
        forw_vde_casadi[i].casadi_fun = &{{ ra.model_name }}_expl_vde_forw;
        forw_vde_casadi[i].casadi_n_in = &{{ ra.model_name }}_expl_vde_forw_n_in;
        forw_vde_casadi[i].casadi_n_out = &{{ ra.model_name }}_expl_vde_forw_n_out;
        forw_vde_casadi[i].casadi_sparsity_in = &{{ ra.model_name }}_expl_vde_forw_sparsity_in;
        forw_vde_casadi[i].casadi_sparsity_out = &{{ ra.model_name }}_expl_vde_forw_sparsity_out;
        forw_vde_casadi[i].casadi_work = &{{ ra.model_name }}_expl_vde_forw_work;
        external_function_casadi_create(&forw_vde_casadi[i]);
    }

    {% if ra.solver_config.hessian_approx == 'EXACT': %} 
    external_function_casadi * hess_vde_casadi;
    hess_vde_casadi = (external_function_casadi *) malloc(sizeof(external_function_casadi)*N);
    for (int i = 0; i < N; ++i) {
        hess_vde_casadi[i].casadi_fun = &{{ ra.model_name }}_expl_ode_hess;
        hess_vde_casadi[i].casadi_n_in = &{{ ra.model_name }}_expl_ode_hess_n_in;
        hess_vde_casadi[i].casadi_n_out = &{{ ra.model_name }}_expl_ode_hess_n_out;
        hess_vde_casadi[i].casadi_sparsity_in = &{{ ra.model_name }}_expl_ode_hess_sparsity_in;
        hess_vde_casadi[i].casadi_sparsity_out = &{{ ra.model_name }}_expl_ode_hess_sparsity_out;
        hess_vde_casadi[i].casadi_work = &{{ ra.model_name }}_expl_ode_hess_work;
        external_function_casadi_create(&hess_vde_casadi[i]);
    }
    {% endif %}
    {% elif ra.solver_config.integrator_type == 'IRK': %}
    // implicit dae
    external_function_casadi * impl_dae_fun;
    impl_dae_fun = (external_function_casadi *) malloc(sizeof(external_function_casadi)*N);
    for (int i = 0; i < N; ++i) {
        impl_dae_fun[i].casadi_fun = &{{ ra.model_name }}_impl_dae_fun;
        impl_dae_fun[i].casadi_work = &{{ ra.model_name }}_impl_dae_fun_work;
        impl_dae_fun[i].casadi_sparsity_in = &{{ ra.model_name }}_impl_dae_fun_sparsity_in;
        impl_dae_fun[i].casadi_sparsity_out = &{{ ra.model_name }}_impl_dae_fun_sparsity_out;
        impl_dae_fun[i].casadi_n_in = &{{ ra.model_name }}_impl_dae_fun_n_in;
        impl_dae_fun[i].casadi_n_out = &{{ ra.model_name }}_impl_dae_fun_n_out;
        external_function_casadi_create(&impl_dae_fun[i]);
    }

    external_function_casadi * impl_dae_fun_jac_x_xdot_z;
    impl_dae_fun_jac_x_xdot_z = (external_function_casadi *) malloc(sizeof(external_function_casadi)*N);
    for (int i = 0; i < N; ++i) {
        impl_dae_fun_jac_x_xdot_z[i].casadi_fun = &{{ ra.model_name }}_impl_dae_fun_jac_x_xdot_z;
        impl_dae_fun_jac_x_xdot_z[i].casadi_work = &{{ ra.model_name }}_impl_dae_fun_jac_x_xdot_z_work;
        impl_dae_fun_jac_x_xdot_z[i].casadi_sparsity_in = &{{ ra.model_name }}_impl_dae_fun_jac_x_xdot_z_sparsity_in;
        impl_dae_fun_jac_x_xdot_z[i].casadi_sparsity_out = &{{ ra.model_name }}_impl_dae_fun_jac_x_xdot_z_sparsity_out;
        impl_dae_fun_jac_x_xdot_z[i].casadi_n_in = &{{ ra.model_name }}_impl_dae_fun_jac_x_xdot_z_n_in;
        impl_dae_fun_jac_x_xdot_z[i].casadi_n_out = &{{ ra.model_name }}_impl_dae_fun_jac_x_xdot_z_n_out;
        external_function_casadi_create(&impl_dae_fun_jac_x_xdot_z[i]);
    }

    external_function_casadi * impl_dae_jac_x_xdot_u_z;
    impl_dae_jac_x_xdot_u_z = (external_function_casadi *) malloc(sizeof(external_function_casadi)*N);
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

    nlp_in = ocp_nlp_in_create(nlp_config, nlp_dims);

    for (int i = 0; i < N; ++i)
        nlp_in->Ts[i] = Tf/N;

    // NLP cost: linear least squares
    // C  // TODO(oj): this can be done using
    // // ocp_nlp_cost_set_model(nlp_config, nlp_dims, nlp_in, i, "Cyt", Cyt);
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
        ocp_nlp_cost_set_model(nlp_config, nlp_dims, nlp_in, i, "W", W);
    }
    // WN
    ocp_nlp_cost_set_model(nlp_config, nlp_dims, nlp_in, N, "W", Q);

    // y_ref
    for (int i = 0; i <= N; ++i)
        ocp_nlp_cost_set_model(nlp_config, nlp_dims, nlp_in, i, "yref", yref);

    // NLP dynamics
    int set_fun_status;
    for (int i = 0; i < N; ++i) {
    {% if ra.solver_config.integrator_type == 'ERK': %} 
        set_fun_status = ocp_nlp_dynamics_set_model(nlp_config, nlp_in, i, "expl_vde_for", &forw_vde_casadi[i]);
        if (set_fun_status != 0) { printf("Error while setting expl_vde_for[%i]\n", i);  exit(1); }
        {% if ra.solver_config.hessian_approx == 'EXACT': %} 
            set_fun_status = ocp_nlp_dynamics_set_model(nlp_config, nlp_in, i, "expl_ode_hes", &hess_vde_casadi[i]);
            if (set_fun_status != 0) { printf("Error while setting expl_ode_hes[%i]\n", i);  exit(1); }
        {% endif %}
    {% elif ra.solver_config.integrator_type == 'IRK': %} 
        set_fun_status = ocp_nlp_dynamics_set_model(nlp_config, nlp_in, i, "impl_ode_fun", &impl_dae_fun[i]);
        if (set_fun_status != 0) { printf("Error while setting impl_dae_fun[%i]\n", i);  exit(1); }
        set_fun_status = ocp_nlp_dynamics_set_model(nlp_config, nlp_in, i, "impl_ode_fun_jac_x_xdot", &impl_dae_fun_jac_x_xdot_z[i]);
        if (set_fun_status != 0) { printf("Error while setting impl_dae_fun_jac_x_xdot_z[%i]\n", i);  exit(1); }
        set_fun_status = ocp_nlp_dynamics_set_model(nlp_config, nlp_in, i, "impl_ode_jac_x_xdot_u", &impl_dae_jac_x_xdot_u_z[i]);
        if (set_fun_status != 0) { printf("Error while setting impl_dae_jac_x_xdot_u_z[%i]\n", i);  exit(1); }
    {% endif %}
    }

    // NLP constraints
    // TODO(oj): remove this when idxb setter available
    ocp_nlp_constraints_bgh_model **constraints = (ocp_nlp_constraints_bgh_model **) nlp_in->constraints;
	ocp_nlp_constraints_bgh_dims **constraints_dims = (ocp_nlp_constraints_bgh_dims **) nlp_dims->constraints;

    // bounds
    for (int i = 0; i < nb[0]; ++i)
        constraints[0]->idxb[i] = idxb_0[i];
	ocp_nlp_constraints_bounds_set(nlp_config, nlp_dims, nlp_in, 0, "lb", lb0);
	ocp_nlp_constraints_bounds_set(nlp_config, nlp_dims, nlp_in, 0, "ub", ub0);

    for (int i = 1; i < N; ++i)
    {
        for (int j = 0; j < nb[i]; ++j)
            constraints[i]->idxb[j] = idxb[j];
        ocp_nlp_constraints_bounds_set(nlp_config, nlp_dims, nlp_in, i, "lb", lb);
        ocp_nlp_constraints_bounds_set(nlp_config, nlp_dims, nlp_in, i, "ub", ub);
    }

    for (int j = 0; j < nb[N]; ++j)
        constraints[N]->idxb[j] = idxb[j];
    ocp_nlp_constraints_bounds_set(nlp_config, nlp_dims, nlp_in, N, "lb", lbN);
	ocp_nlp_constraints_bounds_set(nlp_config, nlp_dims, nlp_in, N, "ub", ubN);

    nlp_opts = ocp_nlp_opts_create(nlp_config, nlp_dims);

    {% if ra.solver_config.nlp_solver_type == 'SQP': %}

    int maxIter = max_num_sqp_iterations;
    double min_res_g = 1e-9;
    double min_res_b = 1e-9;
    double min_res_d = 1e-9;
    double min_res_m = 1e-9;

    ocp_nlp_opts_set(nlp_config, nlp_opts, "maxIter", &maxIter);
    ocp_nlp_opts_set(nlp_config, nlp_opts, "min_res_g", &min_res_g);
    ocp_nlp_opts_set(nlp_config, nlp_opts, "min_res_b", &min_res_b);
    ocp_nlp_opts_set(nlp_config, nlp_opts, "min_res_d", &min_res_d);
    ocp_nlp_opts_set(nlp_config, nlp_opts, "min_res_m", &min_res_m);


    {% else: %}
    ocp_nlp_sqp_rti_opts *sqp_opts = (ocp_nlp_sqp_rti_opts *) nlp_opts;
    {% endif %}
    for (int i = 0; i < N; ++i)
    {
        {% if ra.solver_config.hessian_approx == 'EXACT': %}
        // TODO(oj): is the following needed, and what does it do? do we
        // ((ocp_nlp_dynamics_cont_opts *) sqp_opts->dynamics[i])->compute_hess = true;
        int num_steps = 5;
        bool sens_hess = true;
        bool sens_adj = true;

        ocp_nlp_dynamics_opts_set(nlp_config, nlp_opts, i, "num_steps", &num_steps);
        // ocp_nlp_dynamics_opts_set(nlp_config, nlp_opts, i, "ns", &ns);
        ocp_nlp_dynamics_opts_set(nlp_config, nlp_opts, i, "sens_hess", &sens_hess);
        ocp_nlp_dynamics_opts_set(nlp_config, nlp_opts, i, "sens_adj", &sens_adj);
        {% endif %}
    }

    nlp_out = ocp_nlp_out_create(nlp_config, nlp_dims);
    for (int i = 0; i <= N; ++i)
        blasfeo_dvecse(nu[i]+nx[i], 0.0, nlp_out->ux+i, 0);

    nlp_solver = ocp_nlp_solver_create(nlp_config, nlp_dims, nlp_opts);

    // *_nlp_solver = nlp_solver;
    // *_nlp_in = nlp_in; 
    // *_nlp_out = nlp_out;
    // *_nlp_opts = nlp_opts;
    // *_nlp_config = nlp_config;
    // *_solver_plan = plan;
    // *_nlp_dims = nlp_dims;

    return status;
}

int acados_solve() {

    // NLP solution
    acados_timer timer;
    double kkt_norm_inf = 1e12, elapsed_time;
    
    int solver_status = 0, iteration_number = 0;

    {% if ra.solver_config.nlp_solver_type == 'SQP': %}
    while (kkt_norm_inf > 1e-9) {
        acados_tic(&timer);
        solver_status = ocp_nlp_solve(nlp_solver, nlp_in, nlp_out);
        elapsed_time = acados_toc(&timer);
        kkt_norm_inf = nlp_out->inf_norm_res;
        printf(" iteration %2d | time  %f |  KKT %e\n", iteration_number, elapsed_time, kkt_norm_inf);
        iteration_number++;

        if (iteration_number >= 100)
            break;
    }
    {% else: %}
    solver_status = ocp_nlp_solve(nlp_solver, nlp_in, nlp_out);
    {% endif %}

    printf("\n--- solution ---\n");
    ocp_nlp_out_print(nlp_solver->dims, nlp_out);

    return solver_status;
}

int acados_free() {

    // free memory
    ocp_nlp_opts_free(nlp_opts);
    ocp_nlp_in_free(nlp_in);
    ocp_nlp_out_free(nlp_out);
    ocp_nlp_free(nlp_solver);
    ocp_nlp_dims_free(nlp_dims);
    ocp_nlp_config_free(nlp_solver_plan, nlp_config);
    ocp_nlp_plan_free(nlp_solver_plan);

    // // free external function 
    // {% if ra.solver_config.integrator_type == 'IRK': %}
    // for(int i = 0; i < N; i++) {
    //     external_function_casadi_free(&impl_dae_fun[i]);
    //     external_function_casadi_free(&impl_dae_fun_jac_x_xdot_z[i]);
    //     external_function_casadi_free(&impl_dae_jac_x_xdot_u_z[i]);
    // }
    // {% else: %}
    // for(int i = 0; i < N; i++) {
    //     external_function_casadi_free(&forw_vde_casadi[i]);
    // {% if ra.solver_config.hessian_approx == 'EXACT': %}
    //     external_function_casadi_free(&hess_vde_casadi[i]);
    // {% endif %}
    // }
    // {% endif %}
    
    return 0;
}

ocp_nlp_in * acados_get_nlp_in() { return  nlp_in; }
ocp_nlp_out * acados_get_nlp_out() { return  nlp_out; }
ocp_nlp_solver * acados_get_nlp_solver() { return  nlp_solver; }
ocp_nlp_solver_config * acados_get_nlp_config() { return  nlp_config; }
void * acados_get_nlp_opts() { return  nlp_opts; }
ocp_nlp_dims * acados_get_nlp_dims() { return  nlp_dims; }
