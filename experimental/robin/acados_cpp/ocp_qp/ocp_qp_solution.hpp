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

#ifndef INTERFACES_ACADOS_CPP_OCP_QP_OCP_QP_SOLUTION_HPP_
#define INTERFACES_ACADOS_CPP_OCP_QP_OCP_QP_SOLUTION_HPP_

#include <memory>
#include <vector>

#include "acados/ocp_qp/ocp_qp_common.h"

using std::vector;

namespace acados
{
class ocp_qp_solution
{
 public:
    explicit ocp_qp_solution(std::unique_ptr<ocp_qp_out> solution);

    ocp_qp_solution(ocp_qp_solution&& other);

    ocp_qp_solution(const ocp_qp_solution& other);

    std::vector<std::vector<double>> states();
    std::vector<std::vector<double>> controls();
    std::vector<std::vector<double>> lag_mul_dynamics();
    std::vector<std::vector<double>> lag_mul_bounds();
    std::vector<std::vector<double>> lag_mul_constraints();
    ocp_qp_info info();

    const int N;

 private:
    std::unique_ptr<ocp_qp_out> qp_out;
};

}  // namespace acados

#endif  // INTERFACES_ACADOS_CPP_OCP_QP_OCP_QP_SOLUTION_HPP_
