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

#include "acados_cpp/ocp_qp/ocp_qp_dimensions.hpp"

#include <algorithm>

using std::map;
using std::string;
using std::vector;

namespace acados
{
bool valid_dimensions(map<string, vector<uint>> dims)
{
    auto expected_size = dims["nx"].size();
    if (expected_size <= 0) return false;

    for (auto dim : dimension_keys)
    {
        if (!dims.count(dim)) return false;
        if (dims[dim].size() != expected_size) return false;
    }
    return true;
}

std::unique_ptr<ocp_qp_dims> make_dimensions_ptr(map<string, vector<uint>> dims)
{
    if (!valid_dimensions(dims)) throw std::invalid_argument("Invalid dimensions container.");

    int N = dims["nx"].size() - 1;

    auto dim = std::unique_ptr<ocp_qp_dims>(ocp_qp_dims_create(N));

    std::copy_n(std::begin(dims["nx"]), N + 1, dim->nx);
    std::copy_n(std::begin(dims["nu"]), N + 1, dim->nu);
    std::copy_n(std::begin(dims["nbx"]), N + 1, dim->nbx);
    std::copy_n(std::begin(dims["nbu"]), N + 1, dim->nbu);
    std::copy_n(std::begin(dims["ng"]), N + 1, dim->ng);

    return dim;
}

}  // namespace acados
