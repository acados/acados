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

#ifndef INTERFACES_ACADOS_CPP_OCP_QP_UTILS_HPP_
#define INTERFACES_ACADOS_CPP_OCP_QP_UTILS_HPP_

#include <string>
#include <utility>
#include <vector>

namespace std
{
std::string to_string(std::pair<uint, uint> p)
{
    return "( " + std::to_string(p.first) + ", " + std::to_string(p.second) + " )";
}

template <typename T>
std::string to_string(std::vector<T> v)
{
    std::string result_string = " vector of length " + std::to_string(v.size()) + ": [\n ";
    for (auto it : v)
    {
        result_string += std::to_string(it) + ", ";
    }
    return result_string + "]\n";
}
}  // namespace std

#endif  // INTERFACES_ACADOS_CPP_OCP_QP_UTILS_HPP_
