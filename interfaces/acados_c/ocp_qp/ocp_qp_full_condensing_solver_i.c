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

#include "acados_c/ocp_qp/ocp_qp_full_condensing_solver.h"

#include "acados_c/ocp_qp/ocp_qp_full_condensing.h"



void *ocp_qp_full_condensing_solver_copy_args(ocp_qp_dims *dims, void *raw_memory, void *source_)
{
    ocp_qp_full_condensing_solver_opts *source = (ocp_qp_full_condensing_solver_opts *) source_;
    ocp_qp_full_condensing_solver_opts *dest;

    dense_qp_dims ddims;
    compute_dense_qp_dims(dims, &ddims);

    dest = ocp_qp_full_condensing_solver_assign_args(dims, source->solver, raw_memory);

    ocp_qp_full_condensing_copy_args(dims, dest->cond_opts, source->cond_opts);

    // dest->solver->copy_args(&ddims, dest->qp_solver_opts, source->qp_solver_opts);

    return (void*)dest;
}