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

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <numeric>
#include <stdexcept>

#include "acados_cpp/ocp_qp/ocp_qp.hpp"

#include "acados/utils/print.h"
#include "acados_c/ocp_qp_interface.h"
#include "acados_c/options_interface.h"

#include "acados_cpp/ocp_bounds.hpp"
#include "acados_cpp/ocp_dimensions.hpp"
#include "acados_cpp/ocp_qp/hpipm_helper.hpp"
#include "acados_cpp/utils.hpp"

using std::map;
using std::string;
using std::vector;

namespace acados
{
ocp_qp::ocp_qp(std::vector<int> nx, std::vector<int> nu, std::vector<int> ng, std::vector<int> ns)
    : N(nx.size() - 1), qp(nullptr), solver(nullptr), needs_initializing_(true)
{
    // Number of controls on last stage should be zero;
    if (!nu.empty()) nu.back() = 0;

    auto dim = create_ocp_qp_dimensions_ptr(
        {{"nx", nx}, {"nu", nu}, {"nbx", nx}, {"nbu", nu}, {"ng", ng}, {"ns", ns}});

    qp = std::unique_ptr<ocp_qp_in>(ocp_qp_in_create(NULL, dim.get()));

    for (int stage = 0; stage <= N; ++stage)
    {
        cached_bounds["lbx"].push_back(vector<double>(qp->dim->nx[stage], -INFINITY));
        cached_bounds["ubx"].push_back(vector<double>(qp->dim->nx[stage], +INFINITY));
        cached_bounds["lbu"].push_back(vector<double>(qp->dim->nu[stage], -INFINITY));
        cached_bounds["ubu"].push_back(vector<double>(qp->dim->nu[stage], +INFINITY));
    }

    reset_bounds();
}

ocp_qp::ocp_qp(int N, int nx, int nu, int ng, int ns)
    : ocp_qp(std::vector<int>(N + 1, nx), std::vector<int>(N + 1, nu), std::vector<int>(N + 1, ng),
             std::vector<int>(N + 1, ns))
{
}

void ocp_qp::set_field(string field, vector<double> v)
{
    int last_stage = N;
    if (field == "A" || field == "B" || field == "b") last_stage = N - 1;

    for (int stage = 0; stage <= last_stage; ++stage) set_field(field, stage, v);
}

/*
 * Update one field with some values. Matrices are passed in column-major ordering.
 */
void ocp_qp::set_field(string field, int stage, std::vector<double> v)
{
    if (!in_range(field, stage)) throw std::out_of_range("Stage index should be in [0, N].");
    if (!match(shape_of_field(field, stage), v.size()))
        throw std::invalid_argument("I need " + to_string(shape_of_field(field, stage)) +
                                    " elements but got " + std::to_string(v.size()) + ".");

    if (field == "lbx" || field == "ubx" || field == "lbu" || field == "ubu")
    {
        cached_bounds.at(field).at(stage) = v;

        bool needs_resetting = false;
        if (field == "lbx" || field == "ubx")
        {
            if (get_bound_indices("lbx", stage) !=
                calculate_idxb(cached_bounds["lbx"].at(stage), cached_bounds["ubx"].at(stage)))
                needs_resetting = true;
        }
        else if (field == "lbu" || field == "ubu")
            if (get_bound_indices("lbu", stage) !=
                calculate_idxb(cached_bounds["lbu"].at(stage), cached_bounds["ubu"].at(stage)))
                needs_resetting = true;
        if (needs_resetting) reset_bounds();
    }
    else if (field == "Q")
    {
        d_cvt_colmaj_to_ocp_qp_Q(stage, v.data(), qp.get());
    }
    else if (field == "S")
    {
        d_cvt_colmaj_to_ocp_qp_S(stage, v.data(), qp.get());
    }
    else if (field == "R")
    {
        d_cvt_colmaj_to_ocp_qp_R(stage, v.data(), qp.get());
    }
    else if (field == "q")
    {
        d_cvt_colmaj_to_ocp_qp_q(stage, v.data(), qp.get());
    }
    else if (field == "r")
    {
        d_cvt_colmaj_to_ocp_qp_r(stage, v.data(), qp.get());
    }
    else if (field == "A")
    {
        d_cvt_colmaj_to_ocp_qp_A(stage, v.data(), qp.get());
    }
    else if (field == "B")
    {
        d_cvt_colmaj_to_ocp_qp_B(stage, v.data(), qp.get());
    }
    else if (field == "b")
    {
        d_cvt_colmaj_to_ocp_qp_b(stage, v.data(), qp.get());
    }
    else if (field == "C")
    {
        d_cvt_colmaj_to_ocp_qp_C(stage, v.data(), qp.get());
    }
    else if (field == "D")
    {
        d_cvt_colmaj_to_ocp_qp_D(stage, v.data(), qp.get());
    }
    else if (field == "lg")
    {
        d_cvt_colmaj_to_ocp_qp_lg(stage, v.data(), qp.get());
    }
    else if (field == "ug")
    {
        d_cvt_colmaj_to_ocp_qp_ug(stage, v.data(), qp.get());
    }
    else
    {
        throw std::invalid_argument("OCP QP does not contain field " + field);
    }
}

void ocp_qp::initialize_solver(string solver_name, map<string, option_t *> options)
{
    // check if solver is available
    ocp_qp_solver_plan_t plan;
    try
    {
        plan = available_solvers.at(solver_name);
    }
    catch (std::exception e)
    {
        throw std::invalid_argument("QP solver '" + solver_name + "' is not available.");
    }

    squeeze_dimensions(cached_bounds);

    config.reset(ocp_qp_config_create(plan));
    args.reset(ocp_qp_opts_create(config.get(), qp->dim));
    process_options(solver_name, options, args.get());

    solver.reset(ocp_qp_create(config.get(), qp->dim, args.get()));

    needs_initializing(false);
    cached_solver = solver_name;
}

void ocp_qp::reset_bounds()
{
    d_change_bounds_dimensions_ocp_qp(qp->dim->nu, qp->dim->nx, qp.get());

    for (int stage = 0; stage <= N; ++stage)
    {
        std::vector<int> idx_states(nbx().at(stage));
        std::iota(std::begin(idx_states), std::end(idx_states), 0);
        set_bound_indices("x", stage, idx_states);

        std::vector<int> idx_controls(nbu().at(stage));
        std::iota(std::begin(idx_controls), std::end(idx_controls), 0);
        set_bound_indices("u", stage, idx_controls);
    }

    needs_initializing(true);
}

void ocp_qp::set_bound(string bound, int stage, vector<double> new_bound)
{
    if (bound == "lbx")
        d_cvt_colmaj_to_ocp_qp_lbx(stage, new_bound.data(), qp.get());
    else if (bound == "ubx")
        d_cvt_colmaj_to_ocp_qp_ubx(stage, new_bound.data(), qp.get());
    else if (bound == "lbu")
        d_cvt_colmaj_to_ocp_qp_lbu(stage, new_bound.data(), qp.get());
    else if (bound == "ubu")
        d_cvt_colmaj_to_ocp_qp_ubu(stage, new_bound.data(), qp.get());
    else
        throw std::invalid_argument("Expected one of {'lbx', 'ubx', 'lbu', 'ubu'} but got " +
                                    bound + ".");
}

ocp_qp_solution ocp_qp::solve()
{
    if (needs_initializing()) throw std::runtime_error("Initialize solver before calling 'solve'.");

    fill_bounds(cached_bounds);

    auto result = std::unique_ptr<ocp_qp_out>(ocp_qp_out_create(NULL, qp->dim));

    int_t return_code = ocp_qp_solve(solver.get(), qp.get(), result.get());

    if (return_code != ACADOS_SUCCESS)
    {
        if (return_code == ACADOS_MAXITER)
            throw std::runtime_error("QP solver " + cached_solver +
                                     " reached maximum number of iterations.");
        else if (return_code == ACADOS_MINSTEP)
            throw std::runtime_error("QP solver " + cached_solver + " reached minimum step size.");
        else
            throw std::runtime_error("QP solver " + cached_solver +
                                     " failed with solver-specific error code " +
                                     std::to_string(return_code));
    }
    return ocp_qp_solution(std::move(result));
}

vector<int> ocp_qp::get_bound_indices(string name, int stage)
{
    vector<int> idxb;
    if (name == "lbx" || name == "ubx")
    {
        for (int i = 0; i < qp->dim->nb[stage]; ++i)
            if (qp->idxb[stage][i] >= qp->dim->nu[stage])
                idxb.push_back(qp->idxb[stage][i] - qp->dim->nu[stage]);
    }
    else if (name == "lbu" || name == "ubu")
    {
        for (int i = 0; i < qp->dim->nb[stage]; ++i)
            if (qp->idxb[stage][i] < qp->dim->nu[stage]) idxb.push_back(qp->idxb[stage][i]);
    }
    else
        throw std::invalid_argument("Expected bounds name 'lbx', 'ubx', 'lbu' or 'ubu', but got '" +
                                    name + "' instead.");

    return idxb;
}

void ocp_qp::set_bound_indices(string name, int stage, vector<int> v)
{
    int nb_bounds;
    if (name == "x")
        nb_bounds = qp->dim->nbx[stage];
    else if (name == "u")
        nb_bounds = qp->dim->nbu[stage];
    else
        throw std::invalid_argument("Can only set bounds on x and u, you gave: '" + name + "'.");

    if ((size_t) nb_bounds != v.size())
        throw std::invalid_argument("I need " + std::to_string(nb_bounds) + " indices, you gave " +
                                    std::to_string(v.size()) + ".");
    for (int i = 0; i < nb_bounds; ++i)
        if (name == "x")
            qp->idxb[stage][qp->dim->nbu[stage] + i] = qp->dim->nu[stage] + v.at(i);
        else if (name == "u")
            qp->idxb[stage][i] = v.at(i);
}

map<string, std::function<void(int, ocp_qp_in *, double *)>> ocp_qp::extract_functions = {
    {"Q", d_cvt_ocp_qp_to_colmaj_Q},     {"S", d_cvt_ocp_qp_to_colmaj_S},
    {"R", d_cvt_ocp_qp_to_colmaj_R},     {"q", d_cvt_ocp_qp_to_colmaj_q},
    {"r", d_cvt_ocp_qp_to_colmaj_r},     {"A", d_cvt_ocp_qp_to_colmaj_A},
    {"B", d_cvt_ocp_qp_to_colmaj_B},     {"b", d_cvt_ocp_qp_to_colmaj_b},
    {"lbx", d_cvt_ocp_qp_to_colmaj_lbx}, {"ubx", d_cvt_ocp_qp_to_colmaj_ubx},
    {"lbu", d_cvt_ocp_qp_to_colmaj_lbu}, {"ubu", d_cvt_ocp_qp_to_colmaj_ubu},
    {"C", d_cvt_ocp_qp_to_colmaj_C},     {"D", d_cvt_ocp_qp_to_colmaj_D},
    {"lg", d_cvt_ocp_qp_to_colmaj_lg},   {"ug", d_cvt_ocp_qp_to_colmaj_ug}};

vector<vector<double>> ocp_qp::get_field(string field)
{
    int last_index = N;
    if (field == "A" || field == "B" || field == "b") last_index = N - 1;
    vector<vector<double>> result;
    for (int i = 0; i <= last_index; i++)
    {
        auto dims = shape_of_field(field, i);
        vector<double> v(dims.first * dims.second);
        extract_functions[field](i, qp.get(), v.data());
        result.push_back(v);
    }
    return result;
}

map<string, vector<int>> ocp_qp::dimensions()
{
    return {{"nx", nx()}, {"nu", nu()}, {"nbx", nbx()}, {"nbu", nbu()}, {"ng", ng()}};
}

std::vector<int> ocp_qp::nx()
{
    std::vector<int> tmp(N + 1);
    std::copy_n(qp->dim->nx, N + 1, tmp.begin());
    return tmp;
}

std::vector<int> ocp_qp::nu()
{
    std::vector<int> tmp(N + 1);
    std::copy_n(qp->dim->nu, N + 1, tmp.begin());
    return tmp;
}

std::vector<int> ocp_qp::nbx()
{
    std::vector<int> tmp(N + 1);
    std::copy_n(qp->dim->nbx, N + 1, tmp.begin());
    return tmp;
}

std::vector<int> ocp_qp::nbu()
{
    std::vector<int> tmp(N + 1);
    std::copy_n(qp->dim->nbu, N + 1, tmp.begin());
    return tmp;
}

std::vector<int> ocp_qp::ng()
{
    std::vector<int> tmp(N + 1);
    std::copy_n(qp->dim->ng, N + 1, tmp.begin());
    return tmp;
}

bool ocp_qp::in_range(string field, int stage)
{
    return (field == "A" || field == "B" || field == "b") ? (stage < N) : (stage <= N);
}

std::pair<int, int> ocp_qp::shape_of_field(string field, int stage)
{
    if (!in_range(field, stage)) throw std::out_of_range("Stage index should be in [0, N].");

    if (field == "Q")
        return std::make_pair(num_rows_Q(stage, qp->dim), num_cols_Q(stage, qp->dim));
    else if (field == "S")
        return std::make_pair(num_rows_S(stage, qp->dim), num_cols_S(stage, qp->dim));
    else if (field == "R")
        return std::make_pair(num_rows_R(stage, qp->dim), num_cols_R(stage, qp->dim));
    else if (field == "q")
        return std::make_pair(num_elems_q(stage, qp->dim), 1);
    else if (field == "r")
        return std::make_pair(num_elems_r(stage, qp->dim), 1);
    else if (field == "A")
        return std::make_pair(num_rows_A(stage, qp->dim), num_cols_A(stage, qp->dim));
    else if (field == "B")
        return std::make_pair(num_rows_B(stage, qp->dim), num_cols_B(stage, qp->dim));
    else if (field == "b")
        return std::make_pair(num_elems_b(stage, qp->dim), 1);
    else if (field == "lbx")
        return std::make_pair(qp->dim->nx[stage], 1);
    else if (field == "ubx")
        return std::make_pair(qp->dim->nx[stage], 1);
    else if (field == "lbu")
        return std::make_pair(qp->dim->nu[stage], 1);
    else if (field == "ubu")
        return std::make_pair(qp->dim->nu[stage], 1);
    else if (field == "C")
        return std::make_pair(num_rows_C(stage, qp->dim), num_cols_C(stage, qp->dim));
    else if (field == "D")
        return std::make_pair(num_rows_D(stage, qp->dim), num_cols_D(stage, qp->dim));
    else if (field == "lg")
        return std::make_pair(num_elems_lg(stage, qp->dim), 1);
    else if (field == "ug")
        return std::make_pair(num_elems_ug(stage, qp->dim), 1);
    else
        throw std::invalid_argument("OCP QP does not contain field " + field);
}

void ocp_qp::change_bound_dimensions(vector<int> nbx, vector<int> nbu)
{
    d_change_bounds_dimensions_ocp_qp(nbu.data(), nbx.data(), qp.get());
}

void ocp_qp::needs_initializing(bool flag) { needs_initializing_ = flag; }

bool ocp_qp::needs_initializing() { return needs_initializing_; }

int ocp_qp::num_stages() { return N; }

}  // namespace acados
