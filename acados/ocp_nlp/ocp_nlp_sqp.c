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

#include "acados/ocp_nlp/ocp_nlp_sqp.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"

// Simple fixed-step Gauss-Newton based SQP routine
int_t ocp_nlp_sqp(const ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out,
                  void *nlp_args_, void *nlp_mem_, void *nlp_work_) {
    int_t return_status = 0;

    // Initialize
    ocp_nlp_sqp_memory *sqp_mem = (ocp_nlp_sqp_memory *)nlp_mem_;

    // SQP iterations
    int_t max_sqp_iterations = ((ocp_nlp_sqp_args *)nlp_args_)->maxIter;
    for (int_t sqp_iter = 0; sqp_iter < max_sqp_iterations; sqp_iter++) {
        // Compute/update quadratic approximation
        // sqp_mem->sensitivity_method->fun(.....);

        // Prepare QP

        // Solve QP
        int_t qp_status = sqp_mem->qp_solver->fun(
            sqp_mem->qp_solver->qp_in, sqp_mem->qp_solver->qp_out,
            sqp_mem->qp_solver->args, sqp_mem->qp_solver->mem,
            sqp_mem->qp_solver->work);
        if (qp_status) {
            return_status = qp_status;
        }

        // Update optimization variables (globalization)

    }

    // Post-process solution

    return return_status;
}

void ocp_nlp_sqp_create_memory(const ocp_nlp_in *in, void *args_,
                               void *memory_) {

}

void ocp_nlp_sqp_free_memory(void *mem_) {

}
