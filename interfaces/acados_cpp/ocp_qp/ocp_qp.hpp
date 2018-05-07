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

#ifndef INTERFACES_ACADOS_CPP_OCP_QP_OCP_QP_HPP_
#define INTERFACES_ACADOS_CPP_OCP_QP_OCP_QP_HPP_

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "acados/ocp_qp/ocp_qp_common.h"

#include "acados_c/ocp_qp_interface.h"
#include "acados_cpp/ocp_qp/ocp_qp_solution.hpp"
#include "acados_cpp/ocp_qp/options.hpp"

namespace acados
{
class ocp_qp
{
 public:
    ocp_qp(std::vector<uint> nx, std::vector<uint> nu, std::vector<uint> nbx, std::vector<uint> nbu,
           std::vector<uint> ng);

    explicit ocp_qp(std::map<std::string, std::vector<uint>>);

    ocp_qp(uint N, uint nx, uint nu, uint nbx = 0, uint nbu = 0, uint ng = 0);

    void set(std::string field, uint stage, std::vector<double> v);
    void set(std::string field, std::vector<double> v);

    void initialize_solver(std::string solver_name, std::map<std::string, option_t *> options = {});

    ocp_qp_solution solve();

    std::vector<std::vector<double>> extract(std::string field);

    std::map<std::string, std::vector<uint>> dimensions();

    std::pair<uint, uint> shape_of(std::string field, uint stage);

    void set_bounds_indices(std::string name, uint stage, std::vector<uint> v);

    std::vector<std::vector<uint>> bounds_indices(std::string name);

    const uint N;

 private:
    vector<uint> idxb(vector<double> lower_bound, vector<double> upper_bound);

    void fill_in_bounds();

    void squeeze_dimensions();

    void expand_dimensions();

    void check_range(std::string field, uint stage);

    void check_num_elements(std::string, uint stage, uint nb_elems);

    void flatten(std::map<std::string, option_t *> &input,
                 std::map<std::string, option_t *> &output);

    std::vector<uint> nx();
    std::vector<uint> nu();
    std::vector<uint> nbx();
    std::vector<uint> nbu();
    std::vector<uint> ng();

    std::map<std::string, std::vector<std::vector<double>>> cached_bounds;

    std::unique_ptr<ocp_qp_in> qp;

    std::unique_ptr<ocp_qp_solver> solver;

    std::unique_ptr<ocp_qp_xcond_solver_config> config;

    std::unique_ptr<void, void (*)(void *)> args{nullptr, std::free};

    std::string cached_solver;

    bool needs_initializing = true;

    static std::map<std::string, std::function<void(int, ocp_qp_in *, double *)>> extract_functions;

    friend std::ostream &operator<<(std::ostream &oss, const ocp_qp &qp);
};

}  // namespace acados

#endif  // INTERFACES_ACADOS_CPP_OCP_QP_OCP_QP_HPP_
