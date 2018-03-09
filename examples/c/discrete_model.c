/*
 *    This file is part of acados.
 *
 *    acados is free software; you can redistribute it and/or
 *    modify it under the terms of thegnU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    acados is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See thegnU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of thegnU Lesser General Public
 *    License along with acados; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "acados/sim/allocate_sim.h"
#include "acados/sim/sim_casadi_wrapper.h"
#include "acados/sim/sim_discrete_model.h"
#include "acados/sim/sim_common.h"
#include "acados/ocp_nlp/allocate_ocp_nlp.h"
#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/ocp_nlp/ocp_nlp_sm_gn.h"
#include "acados/ocp_nlp/ocp_nlp_sqp.h"
#include "acados/utils/print.h"

#include "discrete_model/discrete_model.h"
#include "discrete_model/discrete_model_cost.h"
#include "discrete_model/discrete_model_costN.h"

int main() {

    // Problem-specific data
    int_t N = 1, NX = 1, NU = 1;
    real_t W[] = {1.0, 0.0, 0.0, 1.0};
    real_t ref[] = {0.0, 0.0};
    real_t x0[] = {0.8};
    int_t num_sqp_iterations = 1;

    // Problem dimensions
    int_t *nx, *nu, *nb, *ng, *ny;
    nx = calloc(N+1, sizeof(int_t));
    nu = calloc(N+1, sizeof(int_t));
    nb = calloc(N+1, sizeof(int_t));
    ng = calloc(N+1, sizeof(int_t));
    ny = calloc(N+1, sizeof(int_t));

    for (int_t k = 0; k <= N; k++) {
        nx[k] = NX;
        nu[k] = NU;
        nb[k] = 0;
        ng[k] = 0;
        ny[k] = NX+NU;
    }
    nb[0] = NX;
    nu[N] = 0;
    ny[N] = NX;

    ocp_nlp_in nlp;
    allocate_ocp_nlp_in(N, nx, nu, nb, ng, 0, &nlp);
    int_t idxb0[] = {0};
    nlp.idxb[0] = idxb0;

    ocp_nlp_ls_cost ls_cost;
    allocate_ls_cost(N, nx, nu, ny, &ls_cost);
    nlp.cost = &ls_cost;
    // Cost function
    for (int_t k = 0; k <= N; k++) {
        memcpy(ls_cost.W[k], W, (nx[k]+nu[k])*(nx[k]+nu[k])*sizeof(real_t));
        memcpy(ls_cost.y_ref[k], ref, (nx[k]+nu[k])*sizeof(real_t));
        ls_cost.fun[k]->args->fun = &discrete_model_cost;
        ls_cost.fun[k]->args->dims = &discrete_model_cost_work;
        ls_cost.fun[k]->args->sparsity = &discrete_model_cost_sparsity_out;
        casadi_wrapper_initialize(ls_cost.fun[k]->in, ls_cost.fun[k]->args, &ls_cost.fun[k]->work);
    }
    ls_cost.fun[N]->args->fun = &discrete_model_costN;
    ls_cost.fun[N]->args->dims = &discrete_model_costN_work;
    ls_cost.fun[N]->args->sparsity = &discrete_model_costN_sparsity_out;
    casadi_wrapper_initialize(ls_cost.fun[N]->in, ls_cost.fun[N]->args, &ls_cost.fun[N]->work);

    // Dynamics
    sim_solver **simulators = (sim_solver **) nlp.sim;
    for (int_t k = 0; k < N; k++) {
        sim_in *sim = simulators[k]->in;
        sim->nx = NX;
        sim->nu = NU;
        sim->discrete_model = discrete_model;
        simulators[k]->fun = sim_discrete_model;
        simulators[k]->args = NULL;
        simulators[k]->mem = NULL;
        simulators[k]->work = NULL;
    }

    ocp_nlp_out output;
    allocate_ocp_nlp_out(&nlp, &output);

    ocp_nlp_sm sensitivity_method;
    sensitivity_method.fun = &ocp_nlp_sm_gn;
    sensitivity_method.initialize = &ocp_nlp_sm_gn_initialize;
    sensitivity_method.destroy = &ocp_nlp_sm_gn_destroy;

    ocp_nlp_sqp_args *nlp_args = ocp_nlp_sqp_create_arguments();
    nlp_args->maxIter = num_sqp_iterations;
    nlp_args->sensitivity_method = &sensitivity_method;
    ocp_qp_in *qp_in = ocp_qp_in_create(N, nx, nu, nb, ng);
    nlp_args->qp_solver = create_ocp_qp_solver(qp_in, "condensing_qpoases", NULL);
    ocp_nlp_sqp_memory *nlp_mem = ocp_nlp_sqp_create_memory(&nlp, nlp_args);
    void *nlp_sqp_work = ocp_nlp_sqp_create_workspace(&nlp, nlp_args);

    ocp_nlp_sqp_initialize(&nlp, nlp_args, (void **)&nlp_mem, &nlp_sqp_work);

    // Initial guess
    for (int_t k = 0; k <= N; k++) {
        for (int_t j = 0; j < nx[k]; j++)
            nlp_mem->common->x[k][j] = 0.0;
        for (int_t j = 0; j < nu[k]; j++)
            nlp_mem->common->u[k][j] = 0.0;
    }

    nlp.lb[0] = x0;
    nlp.ub[0] = x0;
    int_t status = ocp_nlp_sqp(&nlp, &output, nlp_args, nlp_mem, nlp_sqp_work);
    printf("\n\nstatus = %i\n\n", status);

    for (int_t k = 0; k <= N; k++) {
        char states_name[MAX_STR_LEN], controls_name[MAX_STR_LEN];
        snprintf(states_name, sizeof(states_name), "x%d.txt", k);
        print_matrix_name("stdout", states_name, output.x[k], 1, nx[k]);
        snprintf(controls_name, sizeof(controls_name), "u%d.txt", k);
        print_matrix_name("stdout", controls_name, output.u[k], 1, nu[k]);
    }
}
