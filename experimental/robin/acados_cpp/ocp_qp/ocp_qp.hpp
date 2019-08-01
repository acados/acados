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
#include "acados_cpp/ocp.hpp"
#include "acados_cpp/ocp_qp/ocp_qp_solution.hpp"
#include "acados_cpp/options.hpp"

namespace acados
{
class ocp_qp : private ocp
{
 public:
    ocp_qp(std::vector<int> nx, std::vector<int> nu, std::vector<int> ng, std::vector<int> ns);

    ocp_qp(int N, int nx, int nu, int ng = 0, int ns = 0);

    std::vector<std::vector<double>> get_field(std::string field);

    // Update all fields with the same values. Matrices are passed in column-major ordering.
    void set_field(std::string field, int stage, std::vector<double> v);

    // Update a single field. Matrices are passed in column-major ordering.
    void set_field(std::string field, std::vector<double> v);

    std::pair<int, int> shape_of_field(std::string field, int stage);

    void initialize_solver(std::string solver_name, std::map<std::string, option_t *> options = {});

    ocp_qp_solution solve();

    std::map<std::string, std::vector<int>> dimensions();

    const std::vector<std::string> fields{"Q",   "S",   "R",   "q",   "r", "A", "B",  "b",
                                          "lbx", "ubx", "lbu", "ubu", "C", "D", "lg", "ug"};

    int num_stages() override;

 private:
    void reset_bounds();

    bool in_range(std::string field, int stage);

    void squeeze_dimensions(std::map<std::string, std::vector<std::vector<double>>> bounds) override
    {
        ocp::squeeze_dimensions(bounds);
    }

    void fill_bounds(std::map<std::string, std::vector<std::vector<double>>> bounds) override
    {
        ocp::fill_bounds(bounds);
    }

    void change_bound_dimensions(std::vector<int> nbx, std::vector<int> nbu) override;

    void set_bound(std::string bound, int stage, std::vector<double> new_bound) override;

    std::vector<int> get_bound_indices(std::string name, int stage) override;

    void set_bound_indices(std::string name, int stage, std::vector<int> v) override;

    bool needs_initializing() override;

    void needs_initializing(bool) override;

    std::vector<int> nx();
    std::vector<int> nu();
    std::vector<int> nbx();
    std::vector<int> nbu();
    std::vector<int> ng();

    const int N;

    std::unique_ptr<ocp_qp_in> qp;

    std::unique_ptr<ocp_qp_solver> solver;

    const std::map<std::string, ocp_qp_solver_plan> available_solvers = {
        {"condensing_hpipm", {FULL_CONDENSING_HPIPM}},
        {"sparse_hpipm", {PARTIAL_CONDENSING_HPIPM}},
#ifdef ACADOS_WITH_HPMPC
        {"hpmpc", {PARTIAL_CONDENSING_HPMPC}},
#endif
#ifdef ACADOS_WITH_OOQP
        {"ooqp", {PARTIAL_CONDENSING_OOQP}},
#endif
#ifdef ACADOS_WITH_OSQP
        {"osqp", {PARTIAL_CONDENSING_OSQP}},
#endif
#ifdef ACADOS_WITH_QPDUNES
        {"qpdunes", {PARTIAL_CONDENSING_QPDUNES}},
#endif
#ifdef ACADOS_WITH_QPOASES
        {"qpoases", {FULL_CONDENSING_QPOASES}},
#endif
#ifdef ACADOS_WITH_QORE
        {"qore", {FULL_CONDENSING_QORE}},
#endif
    };

    std::unique_ptr<ocp_qp_xcond_solver_config> config;

    std::unique_ptr<void, void (*)(void *)> args{nullptr, std::free};

    std::string cached_solver;

    std::map<std::string, std::vector<std::vector<double>>> cached_bounds;

    bool needs_initializing_;

    static std::map<std::string, std::function<void(int, ocp_qp_in *, double *)>> extract_functions;
};

}  // namespace acados

#endif  // INTERFACES_ACADOS_CPP_OCP_QP_OCP_QP_HPP_
