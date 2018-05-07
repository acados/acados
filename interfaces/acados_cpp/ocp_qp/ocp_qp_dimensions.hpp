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

#ifndef INTERFACES_ACADOS_CPP_OCP_QP_OCP_QP_DIMENSIONS_HPP_
#define INTERFACES_ACADOS_CPP_OCP_QP_OCP_QP_DIMENSIONS_HPP_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "acados_c/ocp_qp_interface.h"

namespace acados
{
const std::vector<std::string> dimension_keys{"nx", "nu", "nbx", "nbu", "ng"};

bool valid_dimensions(std::map<std::string, std::vector<uint>> dims);

std::unique_ptr<ocp_qp_dims> make_dimensions_ptr(std::map<std::string, std::vector<uint>> dims);

}  // namespace acados

#endif  // INTERFACES_ACADOS_CPP_OCP_QP_OCP_QP_DIMENSIONS_HPP_
