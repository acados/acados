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

#ifndef TEST_OCP_QP_CONDENSING_TEST_HELPER_H_
#define TEST_OCP_QP_CONDENSING_TEST_HELPER_H_

#include <string>

#include "acados/ocp_qp/condensing.h"
#include "test/test_utils/eigen.h"

void calculate_num_state_bounds(const ocp_qp_in *in, condensing_workspace *work);

int_t get_num_condensed_vars(const ocp_qp_in *in);

int_t get_num_constraints(const ocp_qp_in *in, condensing_workspace *work);

void allocateUnconstrainedQPData(int_t N, int_t nx, int_t nu, ocp_qp_in * const qp);

void allocateCondensingData(const ocp_qp_in * const qp_in, condensing_in *in,
    condensing_out *out, condensing_workspace *work);

void fillWithUnconstrainedData(ocp_qp_in *qp, Eigen::VectorXd *x0, std::string scenario);

void fillWithBoundsData(ocp_qp_in *qp, int_t N, int_t nx, int_t nu, std::string scenario);

void fillWithGeneralConstraintsData(ocp_qp_in *qp, int_t N, int_t nx, int_t nu,
    std::string scenario);

#endif  // TEST_OCP_QP_CONDENSING_TEST_HELPER_H_
