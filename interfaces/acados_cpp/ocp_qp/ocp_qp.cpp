
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <stdexcept>
#include <cstdlib>
#include <cmath>

#include "acados_cpp/ocp_qp/ocp_qp.hpp"

#include "acados/utils/print.h"
#include "acados_c/ocp_qp_interface.h"
#include "acados_c/options.h"

#include "acados_cpp/ocp_qp/hpipm_helper.hpp"
#include "acados_cpp/ocp_qp/utils.hpp"

using std::map;
using std::string;
using std::vector;

namespace acados {

ocp_qp::ocp_qp(std::vector<uint> nx, std::vector<uint> nu, std::vector<uint> nbx,
             std::vector<uint> nbu, std::vector<uint> ng) : N(nx.size()-1), qp(nullptr), solver(nullptr) {

    if (N <= 0) throw std::invalid_argument("Number of stages must be positive");

    if (nu.at(N) != 0 || nbu.at(N) != 0) {
        nu.at(N) = 0;
        nbu.at(N) = 0;
    }

    uint expected_size = nx.size();
    bool is_valid_nu = (nu.size() == expected_size || nu.size() == expected_size-1);
    bool is_valid_nbu = (nbu.size() == expected_size || nbu.size() == expected_size-1);
    if (!is_valid_nu || nbx.size() != expected_size || !is_valid_nbu || ng.size() != expected_size)
        throw std::invalid_argument("All dimensions should have length N+1");

    auto dim = std::unique_ptr<ocp_qp_dims>(ocp_qp_dims_create(N));

    // states
    std::copy_n(std::begin(nx), N+1, dim->nx);

    // controls
    std::copy_n(std::begin(nu), N+1, dim->nu);

    // bounds
    std::copy_n(std::begin(nbx), N+1, dim->nbx);
    std::copy_n(std::begin(nbu), N+1, dim->nbu);

    std::vector<uint> nb(N+1);
    for(uint i = 0; i <= N; i++)
        nb.at(i) = nbx.at(i) + nbu.at(i);
    std::copy_n(std::begin(nb), N+1, dim->nb);

    // constraints
    std::copy_n(std::begin(ng), N+1, dim->ng);

    // slacks
    std::vector<uint> ns(N+1, 0);
    std::copy_n(std::begin(ns), N+1, dim->ns);

    qp = std::unique_ptr<ocp_qp_in>(ocp_qp_in_create(NULL, dim.get()));

    for (uint i = 0; i < N; ++i) {
        std::vector<uint> idx_states(nx.at(i));
        std::iota(std::begin(idx_states), std::end(idx_states), 0);
        set_bounds_indices("x", i, idx_states);

        std::vector<uint> idx_controls(nu.at(i));
        std::iota(std::begin(idx_controls), std::end(idx_controls), 0);
        set_bounds_indices("u", i, idx_controls);
    }
    std::vector<uint> idx_states(nx.at(N));
    std::iota(std::begin(idx_states), std::end(idx_states), 0);
    set_bounds_indices("x", N, idx_states);

    for (uint stage = 0; stage <= N; ++stage) {
        cached_bounds["lbx"].push_back(vector<double>(qp->dim->nx[stage], -INFINITY));
        cached_bounds["ubx"].push_back(vector<double>(qp->dim->nx[stage], +INFINITY));
        cached_bounds["lbu"].push_back(vector<double>(qp->dim->nu[stage], -INFINITY));
        cached_bounds["ubu"].push_back(vector<double>(qp->dim->nu[stage], +INFINITY));
    }
}

ocp_qp::ocp_qp(map<string, vector<uint>> dims)
    : ocp_qp(dims["nx"], dims["nu"], dims["nbx"], dims["nbu"], dims["ng"]) {}

ocp_qp::ocp_qp(uint N, uint nx, uint nu, uint nbx, uint nbu, uint ng)
    : ocp_qp(std::vector<uint>(N+1, nx), std::vector<uint>(N+1, nu), std::vector<uint>(N+1, nbx),
      std::vector<uint>(N+1, nbu), std::vector<uint>(N+1, ng)) {}

void ocp_qp::set(std::string field, uint stage, std::vector<double> v) {

    check_num_elements(field, stage, v.size());

    if (field == "lbx" || field == "ubx" || field == "lbu" || field == "ubu") {
        cached_bounds.at(field).at(stage) = v;
        auto idxb_stage = bounds_indices(string(1, field.back())).at(stage);
        vector<uint> sq_ids;
        if (field.front() == 'l') {
            string other(field);
            other.front() = 'u';
            sq_ids = idxb(v, cached_bounds.at(other).at(stage));
        }
        else if (field.front() == 'u') {
            string other(field);
            other.front() = 'l';
            sq_ids = idxb(cached_bounds.at(other).at(stage), v);
        }

        if (idxb_stage != sq_ids) {
            expand_dimensions();
            needs_initializing = true;
        }
    } else if (field == "Q") {
        d_cvt_colmaj_to_ocp_qp_Q(stage, v.data(), qp.get());
    } else if (field == "S") {
        d_cvt_colmaj_to_ocp_qp_S(stage, v.data(), qp.get());
    } else if (field == "R") {
        d_cvt_colmaj_to_ocp_qp_R(stage, v.data(), qp.get());
    } else if (field == "q") {
        d_cvt_colmaj_to_ocp_qp_q(stage, v.data(), qp.get());
    } else if (field == "r") {
        d_cvt_colmaj_to_ocp_qp_r(stage, v.data(), qp.get());
    } else if (field == "A") {
        d_cvt_colmaj_to_ocp_qp_A(stage, v.data(), qp.get());
    } else if (field == "B") {
        d_cvt_colmaj_to_ocp_qp_B(stage, v.data(), qp.get());
    } else if (field == "b") {
        d_cvt_colmaj_to_ocp_qp_b(stage, v.data(), qp.get());
    } else if (field == "C") {
        d_cvt_colmaj_to_ocp_qp_C(stage, v.data(), qp.get());
    } else if (field == "D") {
        d_cvt_colmaj_to_ocp_qp_D(stage, v.data(), qp.get());
    } else if (field == "lg") {
        d_cvt_colmaj_to_ocp_qp_lg(stage, v.data(), qp.get());
    } else if (field == "ug") {
        d_cvt_colmaj_to_ocp_qp_ug(stage, v.data(), qp.get());
    } else {
        throw std::invalid_argument("OCP QP does not contain field " + field);
    }
}

void ocp_qp::set(string field, vector<double> v) {
    uint last_stage = N;
    if (field == "A" || field == "B" || field == "b")
        last_stage = N-1;

    for (uint stage = 0; stage <= last_stage; ++stage)
        set(field, stage, v);
}

static ocp_qp_solver_plan string_to_plan(string solver);

void ocp_qp::initialize_solver(string solver_name, map<string, option_t *> options) {

    squeeze_dimensions();
    cached_solver = solver_name;
    ocp_qp_solver_plan plan = string_to_plan(solver_name);

    config.reset(ocp_qp_config_create(plan));

    args.reset(ocp_qp_opts_create(config.get(), qp->dim));

    map<string, option_t *> solver_options;
    auto nested_options = std::make_unique<option<map<string, option_t *>>>(options);
    solver_options[solver_name] = nested_options.get();

    auto flattened_options = map<string, option_t *>();
    flatten(solver_options, flattened_options);

    for (auto opt : flattened_options) {
        string option_name = opt.first;
        option_t *opt_p = opt.second;
        bool found = set_option_int(args.get(), option_name.c_str(), std::to_int(opt_p));
        found |= set_option_double(args.get(), option_name.c_str(), std::to_double(opt_p));
        if (!found)
            throw std::invalid_argument("Option " + option_name + " not known.");
    }
    solver.reset(ocp_qp_create(config.get(), qp->dim, args.get()));
    needs_initializing = false;
}

vector<uint> ocp_qp::idxb(vector<double> lower_bound, vector<double> upper_bound) {
    vector<uint> bound_indices;
    if (lower_bound.size() != upper_bound.size())
        throw std::invalid_argument("Lower bound must have same shape as upper bound.");

    for (uint idx = 0; idx < lower_bound.size(); ++idx)
        if (lower_bound.at(idx) != -INFINITY || upper_bound.at(idx) != +INFINITY)
            bound_indices.push_back(idx); // there is a double-sided bound at this index

    return bound_indices;
}

void ocp_qp::squeeze_dimensions() {
    map<string, vector<vector<uint>>> idxb_new;
    map<string, vector<uint>> nb;
    map<string, vector<vector<double>>> lb;
    map<string, vector<vector<double>>> ub;
    for (string bound : {"x", "u"}) {
        for (uint stage = 0; stage <= N; ++stage) {
            auto idxb_stage = idxb(cached_bounds.at("lb" + bound).at(stage), cached_bounds.at("ub" + bound).at(stage));
            idxb_new[bound].push_back(idxb_stage);
            nb[bound].push_back(idxb_new.at(bound).at(stage).size());
        }
    }

    d_change_bounds_dimensions_ocp_qp(reinterpret_cast<int *>(nb.at("u").data()),
                                      reinterpret_cast<int *>(nb.at("x").data()),
                                      qp.get());

    needs_initializing = true;

    for (string bound : {"x", "u"})
        for (uint stage = 0; stage <= N; ++stage)
            set_bounds_indices(bound, stage, idxb_new.at(bound).at(stage));
}

void ocp_qp::expand_dimensions() {

    map<string, vector<vector<double>>> lb;
    map<string, vector<vector<double>>> ub;

    for (string bound : {"x", "u"}) {
        for (uint stage = 0; stage <= N; ++stage) {
            vector<double> lb_old = extract("lb" + bound).at(stage), ub_old = extract("ub" + bound).at(stage);
            vector<uint> idxb = bounds_indices(bound).at(stage);
            vector<double> lb_new, ub_new;
            for (uint idx = 0, bound_index = 0; idx < dimensions()["n" + bound].at(stage); ++idx) {
                if (bound_index < dimensions()["nb" + bound].at(stage) && idx == idxb.at(bound_index)) {
                    lb_new.push_back(lb_old.at(bound_index));
                    ub_new.push_back(ub_old.at(bound_index));
                    ++bound_index;
                } else {
                    lb_new.push_back(-INFINITY);
                    ub_new.push_back(+INFINITY);
                }
            }
            lb[bound].push_back(lb_new);
            ub[bound].push_back(ub_new);
        }
    }

    d_change_bounds_dimensions_ocp_qp(qp->dim->nu, qp->dim->nx, qp.get());

    needs_initializing = true;

    for (uint i = 0; i < N; ++i) {
        std::vector<uint> idx_states(dimensions()["nx"].at(i));
        std::iota(std::begin(idx_states), std::end(idx_states), 0);
        set_bounds_indices("x", i, idx_states);

        std::vector<uint> idx_controls(dimensions()["nu"].at(i));
        std::iota(std::begin(idx_controls), std::end(idx_controls), 0);
        set_bounds_indices("u", i, idx_controls);
    }
    std::vector<uint> idx_states(dimensions()["nx"].at(N));
    std::iota(std::begin(idx_states), std::end(idx_states), 0);
    set_bounds_indices("x", N, idx_states);
}

void ocp_qp::fill_in_bounds() {
    for (auto it : cached_bounds) {
        for (uint stage = 0; stage <= N; ++stage) {
            auto idxb_stage = bounds_indices(std::string(1, it.first.back())).at(stage);
            auto stage_bound = it.second.at(stage);
            vector<double> new_bound;
            // std::cout << "stage_bound: " << std::to_string(stage_bound);
            for (uint idx = 0; idx < idxb_stage.size(); ++idx) {
                if (it.first.front() == 'l')
                    new_bound.push_back(std::isfinite(stage_bound.at(idxb_stage.at(idx))) ? stage_bound.at(idxb_stage.at(idx)) : ACADOS_NEG_INFTY);
                else if (it.first.front() == 'u')
                    new_bound.push_back(std::isfinite(stage_bound.at(idxb_stage.at(idx))) ? stage_bound.at(idxb_stage.at(idx)) : ACADOS_POS_INFTY);
            }
            if (it.first == "lbx")
                d_cvt_colmaj_to_ocp_qp_lbx(stage, new_bound.data(), qp.get());
            else if (it.first == "ubx")
                d_cvt_colmaj_to_ocp_qp_ubx(stage, new_bound.data(), qp.get());
            else if (it.first == "lbu")
                d_cvt_colmaj_to_ocp_qp_lbu(stage, new_bound.data(), qp.get());
            else if (it.first == "ubu")
                d_cvt_colmaj_to_ocp_qp_ubu(stage, new_bound.data(), qp.get());
        }
    }
}

ocp_qp_solution ocp_qp::solve() {

    if (needs_initializing)
        throw std::runtime_error("Reinitialize solver");

    fill_in_bounds();

    auto result = std::unique_ptr<ocp_qp_out>(ocp_qp_out_create(NULL, qp->dim));

    int_t return_code = ocp_qp_solve(solver.get(), qp.get(), result.get());

    if (return_code != ACADOS_SUCCESS) {
        if (return_code == ACADOS_MAXITER)
            throw std::runtime_error("QP solver " + cached_solver + " reached maximum number of iterations.");
        else if (return_code == ACADOS_MINSTEP)
            throw std::runtime_error("QP solver " + cached_solver + " reached minimum step size.");
        else
            throw std::runtime_error("QP solver " + cached_solver + " failed with solver-specific error code " + std::to_string(return_code));
    }
    return ocp_qp_solution(std::move(result));
}



void ocp_qp::flatten(map<string, option_t *>& input, map<string, option_t *>& output) {
    for (auto opt : input) {
        if (opt.second->nested()) {
            for (auto nested_opt : *opt.second) {
                input.erase(opt.first);
                input[opt.first + "." + nested_opt.first] = nested_opt.second;
                flatten(input, output);
            }
        } else {
            output[opt.first] = opt.second;
        }
    }
}

vector<vector<uint>> ocp_qp::bounds_indices(string name) {

    vector<vector<uint>> idxb;
    if (name == "x") {
        for (uint stage = 0; stage <= N; ++stage) {
            idxb.push_back(vector<uint>());
            for(int i = 0; i < qp->dim->nb[stage]; ++i)
                if (qp->idxb[stage][i] >= qp->dim->nu[stage])
                    idxb.at(stage).push_back(qp->idxb[stage][i] - qp->dim->nu[stage]);
        }
    } else if (name == "u") {
        for (uint stage = 0; stage <= N; ++stage) {
            idxb.push_back(vector<uint>());
            for(int i = 0; i < qp->dim->nb[stage]; ++i)
                if (qp->idxb[stage][i] < qp->dim->nu[stage])
                    idxb.at(stage).push_back(qp->idxb[stage][i]);
        }
    } else throw std::invalid_argument("Can only get bounds from x and u, you gave: '" + name + "'.");
    return idxb;
}

void ocp_qp::set_bounds_indices(string name, uint stage, vector<uint> v) {
    uint nb_bounds;
    if (name == "x")
        nb_bounds = qp->dim->nbx[stage];
    else if (name == "u")
        nb_bounds = qp->dim->nbu[stage];
    else
        throw std::invalid_argument("Can only set bounds on x and u, you gave: '" + name + "'.");

    if (nb_bounds != v.size())
        throw std::invalid_argument("I need " + std::to_string(nb_bounds) + " indices, you gave " + std::to_string(v.size()) + ".");
    for (uint i = 0; i < nb_bounds; ++i)
        if (name == "x")
            qp->idxb[stage][qp->dim->nbu[stage]+i] = qp->dim->nu[stage]+v.at(i);
        else if (name == "u")
            qp->idxb[stage][i] = v.at(i);
}

map<string, std::function<void(int, ocp_qp_in *, double *)>> ocp_qp::extract_functions = {
        {"Q", d_cvt_ocp_qp_to_colmaj_Q},
        {"S", d_cvt_ocp_qp_to_colmaj_S},
        {"R", d_cvt_ocp_qp_to_colmaj_R},
        {"q", d_cvt_ocp_qp_to_colmaj_q},
        {"r", d_cvt_ocp_qp_to_colmaj_r},
        {"A", d_cvt_ocp_qp_to_colmaj_A},
        {"B", d_cvt_ocp_qp_to_colmaj_B},
        {"b", d_cvt_ocp_qp_to_colmaj_b},
        {"lbx", d_cvt_ocp_qp_to_colmaj_lbx},
        {"ubx", d_cvt_ocp_qp_to_colmaj_ubx},
        {"lbu", d_cvt_ocp_qp_to_colmaj_lbu},
        {"ubu", d_cvt_ocp_qp_to_colmaj_ubu},
        {"C", d_cvt_ocp_qp_to_colmaj_C},
        {"D", d_cvt_ocp_qp_to_colmaj_D},
        {"lg", d_cvt_ocp_qp_to_colmaj_lg},
        {"ug", d_cvt_ocp_qp_to_colmaj_ug}
    };


vector< vector<double> > ocp_qp::extract(std::string field) {
    uint last_index = N;
    if (field == "A" || field == "B" || field == "b")
        last_index = N-1;
    vector< vector<double> > result;
    for (uint i = 0; i <= last_index; i++) {
        auto dims = shape_of(field, i);
        vector<double> v(dims.first * dims.second);
        extract_functions[field](i, qp.get(), v.data());
        result.push_back(v);
    }
    return result;
}


map<string, vector<uint>> ocp_qp::dimensions() {
    return {{"nx", nx()}, {"nu", nu()}, {"nbx", nbx()}, {"nbu", nbu()}, {"ng", ng()}};
}


std::vector<uint> ocp_qp::nx() {
    std::vector<uint> tmp(N+1);
    std::copy_n(qp->dim->nx, N+1, tmp.begin());
    return tmp;
}

std::vector<uint> ocp_qp::nu() {
    std::vector<uint> tmp(N+1);
    std::copy_n(qp->dim->nu, N+1, tmp.begin());
    return tmp;
}

std::vector<uint> ocp_qp::nbx() {
    std::vector<uint> tmp(N+1);
    std::copy_n(qp->dim->nbx, N+1, tmp.begin());
    return tmp;
}

std::vector<uint> ocp_qp::nbu() {
    std::vector<uint> tmp(N+1);
    std::copy_n(qp->dim->nbu, N+1, tmp.begin());
    return tmp;
}

std::vector<uint> ocp_qp::ng() {
    std::vector<uint> tmp(N+1);
    std::copy_n(qp->dim->ng, N+1, tmp.begin());
    return tmp;
}

std::ostream& operator<<(std::ostream& oss, const ocp_qp& qp) {
    static char a[1000000];
    print_ocp_qp_in_to_string(a, qp.qp.get());
    oss << a;
    return oss;
}

void ocp_qp::check_range(std::string field, uint stage) {
    uint lower_bound = 0;
    uint upper_bound;
    if (field == "A" || field == "B" || field == "b") {
        upper_bound = N-1;
    } else {
        upper_bound = N;
    }
    if(stage < lower_bound || stage > upper_bound)
        throw std::out_of_range(std::to_string(stage) + " must be in range [" +
              std::to_string(lower_bound) + ", " + std::to_string(upper_bound) + "].");
}

std::pair<uint, uint> ocp_qp::shape_of(std::string field, uint stage) {

    check_range(field, stage);

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

static bool match(std::pair<uint, uint> dims, uint nb_elems) {
    uint nb_expected_elems = dims.first * dims.second;
    if (nb_expected_elems == 0 || nb_expected_elems == nb_elems)
        return true;
    return false;
}

void ocp_qp::check_num_elements(std::string field, uint stage, uint nb_elems) {
    if (!match(shape_of(field, stage), nb_elems))
        throw std::invalid_argument("I need " + std::to_string(shape_of(field, stage)) + " elements but got " + std::to_string(nb_elems) + ".");
}

ocp_qp_solver_plan string_to_plan(string solver) {

    ocp_qp_solver_plan plan;

    if (solver == "condensing_hpipm") {
        plan.qp_solver = FULL_CONDENSING_HPIPM;
    } else if (solver == "sparse_hpipm") {
        plan.qp_solver = PARTIAL_CONDENSING_HPIPM;
    } else if (solver == "hpmpc") {
#ifdef ACADOS_WITH_HPMPC
        plan.qp_solver = PARTIAL_CONDENSING_HPMPC;
#else
        throw std::invalid_argument("Acados compiled without solver HPMPC.");
#endif
    } else if (solver == "ooqp") {
#ifdef ACADOS_WITH_OOQP
        plan.qp_solver = PARTIAL_CONDENSING_OOQP;
#else
        throw std::invalid_argument("Acados compiled without solver OOQP.");
#endif
    } else if (solver == "qpdunes") {
#ifdef ACADOS_WITH_QPDUNES
        plan.qp_solver = PARTIAL_CONDENSING_QPDUNES;
#else 
        throw std::invalid_argument("Acados compiled without solver qpDUNES.");
#endif
    } else if (solver == "qpoases") {
#ifdef ACADOS_WITH_QPOASES
        plan.qp_solver = FULL_CONDENSING_QPOASES;
#else
        throw std::invalid_argument("Acados compiled without solver qpOASES.");
#endif
    } else if (solver == "qore") {
#ifdef ACADOS_WITH_QORE
        plan.qp_solver = FULL_CONDENSING_QORE;
#else
        throw std::invalid_argument("Acados compiled without solver QORE.");
#endif
    } else {
        throw std::invalid_argument("Solver name '" + solver + "' not known.");
    }
    return plan;
}

}  // namespace acados
