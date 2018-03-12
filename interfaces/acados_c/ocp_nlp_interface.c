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

#include "acados_c/ocp_nlp_interface.h"

//external
#include <stdio.h>
#include <stdlib.h>

ocp_nlp_solver_config *ocp_nlp_config_create(ocp_nlp_solver_plan *plan, int N)
{
    int bytes = ocp_nlp_solver_config_calculate_size(N);
	void *config_mem = calloc(1, bytes);
	ocp_nlp_solver_config *config = ocp_nlp_solver_config_assign(N, config_mem);

    if (plan->nlp_solver == SQP_GN)
    {
        // QP solver
        config->qp_solver = ocp_qp_config_create(plan->ocp_qp_solver_plan);
        
        // LS cost
        for (int i = 0; i <= N; ++i)
        {
            ocp_nlp_cost_ls_config_initialize_default(config->cost[i]);
        }

        // Dynamics
        for (int i = 0; i < N; ++i)
        {
		    ocp_nlp_dynamics_config_initialize_default(config->dynamics[i]);
		    config->dynamics[i]->sim_solver = sim_config_create(plan->sim_solver_plan[i]);
        }

        // Constraints
        for (int i = 0; i <= N; ++i)
		    ocp_nlp_constraints_config_initialize_default(config->constraints[i]);
    }
    else
    {
        printf("Solver not available!\n");
        exit(1);
    }
    return config;
}



ocp_nlp_dims *ocp_nlp_dims_create(int N)
{
    return NULL;
}



ocp_nlp_in *ocp_nlp_in_create(ocp_nlp_dims *dims)
{
    return NULL;
}



ocp_nlp_out *ocp_nlp_out_create(ocp_nlp_dims *dims)
{
    return NULL;
}



void *ocp_nlp_opts_create(ocp_nlp_solver_config *plan, ocp_nlp_dims *dims)
{
    return NULL;
}



ocp_nlp_solver *ocp_nlp_create(ocp_nlp_solver_config *plan, ocp_nlp_dims *dims, void *opts_)
{
    return NULL;
}



int ocp_nlp_solve(ocp_nlp_solver *solver, ocp_nlp_in *qp_in, ocp_nlp_out *qp_out)
{
    return 0;
}
