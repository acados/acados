
#include <algorithm>
#include <iostream>
#include <iterator>
#include <stdexcept>

#include "acados_cpp/ocp_qp.hpp"

#include "acados/utils/print.h"
#include "acados_c/ocp_qp.h"
#include "acados_c/options.h"

#include "acados_cpp/hpipm_helper.hpp"
#include "acados_cpp/pair.hpp"

namespace acados {

ocp_qp::ocp_qp(std::vector<uint> nx, std::vector<uint> nu, std::vector<uint> nbx, 
             std::vector<uint> nbu, std::vector<uint> ng) : N(nx.size()-1), qp(nullptr), solver(nullptr) {

    if (N <= 0) throw std::invalid_argument("Number of stages must be positive");

    uint expected_size = nx.size();
    bool is_valid_nu = (nu.size() == expected_size || nu.size() == expected_size-1);
    if (!is_valid_nu || nbx.size() != expected_size || nbu.size() != expected_size
        || ng.size() != expected_size)
            throw std::invalid_argument("Number of stages should be N");

    dim = std::unique_ptr<ocp_qp_dims>(create_ocp_qp_dims(N));

    // states
    std::copy_n(std::begin(nx), N+1, dim->nx);

    // controls
    if (nu.size() == (N+1)) nu.at(N) = 0;
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

    qp = std::unique_ptr<ocp_qp_in>(create_ocp_qp_in(dim.get()));
}

ocp_qp::ocp_qp(map<string, vector<uint>> dims)
    : ocp_qp(dims["nx"], dims["nu"], dims["nbx"], dims["nbu"], dims["ng"]) {}

ocp_qp::ocp_qp(uint N, uint nx, uint nu, uint nbx, uint nbu, uint ng)
    : ocp_qp(std::vector<uint>(N+1, nx), std::vector<uint>(N+1, nu), std::vector<uint>(N+1, nbx),
      std::vector<uint>(N+1, nbu), std::vector<uint>(N+1, ng)) {}

void ocp_qp::set(std::string field, uint stage, std::vector<double> v) {

    check_range(field, stage);
    check_nb_elements(field, stage, v.size());

    if (field == "Q")
        d_cvt_colmaj_to_ocp_qp_Q(stage, v.data(), qp.get());
    else if (field == "S")
        d_cvt_colmaj_to_ocp_qp_S(stage, v.data(), qp.get());
    else if (field == "R")
        d_cvt_colmaj_to_ocp_qp_R(stage, v.data(), qp.get());
    else if (field == "q")
        d_cvt_colmaj_to_ocp_qp_q(stage, v.data(), qp.get());
    else if (field == "r")
        d_cvt_colmaj_to_ocp_qp_r(stage, v.data(), qp.get());
    else if (field == "A")
        d_cvt_colmaj_to_ocp_qp_A(stage, v.data(), qp.get());
    else if (field == "B")
        d_cvt_colmaj_to_ocp_qp_B(stage, v.data(), qp.get());
    else if (field == "b")
        d_cvt_colmaj_to_ocp_qp_b(stage, v.data(), qp.get());
    else if (field == "lbx")
        d_cvt_colmaj_to_ocp_qp_lbx(stage, v.data(), qp.get());
    else if (field == "ubx")
        d_cvt_colmaj_to_ocp_qp_ubx(stage, v.data(), qp.get());
    else if (field == "lbu")
        d_cvt_colmaj_to_ocp_qp_lbu(stage, v.data(), qp.get());
    else if (field == "ubu")
        d_cvt_colmaj_to_ocp_qp_ubu(stage, v.data(), qp.get());
    else if (field == "C")
        d_cvt_colmaj_to_ocp_qp_C(stage, v.data(), qp.get());
    else if (field == "D")
        d_cvt_colmaj_to_ocp_qp_D(stage, v.data(), qp.get());
    else if (field == "lg")
        d_cvt_colmaj_to_ocp_qp_lg(stage, v.data(), qp.get());
    else if (field == "ug")
        d_cvt_colmaj_to_ocp_qp_ug(stage, v.data(), qp.get());
    else
        throw std::invalid_argument("OCP QP does not contain field " + field);
}

void ocp_qp::set(string field, vector<double> v) {
    uint last_stage = N;
    if (field == "A" || field == "B" || field == "b")
        last_stage = N-1;
    
    for (uint stage = 0; stage <= last_stage; ++stage)
        set(field, stage, v);
}

static ocp_qp_solver_plan string_to_plan(string solver);

ocp_qp_solution ocp_qp::solve(string solver_name, map<string, option_t *> options) {

    if (cached_solver != solver_name) {
        cached_solver = solver_name;
        ocp_qp_solver_plan plan = string_to_plan(solver_name);
        solver.reset(ocp_qp_create(&plan, dim.get(), ocp_qp_create_args(&plan, dim.get())));
    }

    solver->fcn_ptrs->initialize_default_args(solver->args);

    map<string, option_t *> solver_options;
    solver_options[solver_name] = new option<map<string, option_t *>>(options);

    auto flattened_options = map<string, option_t *>();
    flatten(solver_options, flattened_options);

    for (auto opt : flattened_options)
        update_option(opt.first, opt.second);
    delete solver_options[solver_name];
    auto result = std::unique_ptr<ocp_qp_out>(create_ocp_qp_out(dim.get()));
    int_t return_code = ocp_qp_solve(solver.get(), qp.get(), result.get());
    if (return_code != 0)
        throw std::runtime_error("qp solver failed with error code " + std::to_string(return_code));
    return ocp_qp_solution(std::move(result));
}


void ocp_qp::update_option(string option_name, option_t *opt_p) {
    bool found = set_option_int(solver->args, option_name.c_str(), std::to_int(opt_p));
    found |= set_option_double(solver->args, option_name.c_str(), std::to_double(opt_p));
    // updated |= set_option_int_array(solver.c_str(), solver->args, option.c_str(), std::to_int(opt_p));
    // updated |= set_option_double_array(solver.c_str(), solver->args, option.c_str(), std::to_int(opt_p));
    if (!found)
        throw std::invalid_argument("Option " + option_name + " not known.");
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

void ocp_qp::state_bounds_indices(uint stage, vector<uint> v) {
    uint nb_state_bounds = qp->dim->nbx[stage];
    if (nb_state_bounds != v.size())
        throw std::invalid_argument("I need " + std::to_string(nb_state_bounds) + " indices.");
    for (uint i = 0; i < nb_state_bounds; ++i)
        qp->idxb[stage][qp->dim->nbu[stage]+i] = qp->dim->nu[stage]+v.at(i);
}

void ocp_qp::control_bounds_indices(uint stage, vector<uint> v) {
    uint nb_control_bounds = qp->dim->nbu[stage];
    if (nb_control_bounds != v.size())
        throw std::invalid_argument("I need " + std::to_string(nb_control_bounds) + " indices.");
    for (uint i = 0; i < nb_control_bounds; ++i)
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
        auto dims = dimensions(field, i);
        vector<double> v(dims.first * dims.second);
        extract_functions[field](i, qp.get(), v.data());
        result.push_back(v);
    }
    return result;
}


const map<string, vector<uint>> ocp_qp::dimensions() const {
    return {{"nx", nx()}, {"nu", nu()}, {"nbx", nbx()}, {"nbu", nbu()}, {"ng", ng()}};
}


std::vector<uint> ocp_qp::nx() const {
    std::vector<uint> tmp(N+1);
    std::copy_n(qp->dim->nx, N+1, tmp.begin());
    return tmp;
}

std::vector<uint> ocp_qp::nu() const {
    std::vector<uint> tmp(N+1);
    std::copy_n(qp->dim->nu, N+1, tmp.begin());
    return tmp;
}

std::vector<uint> ocp_qp::nbx() const {
    std::vector<uint> tmp(N+1);
    std::copy_n(qp->dim->nbx, N+1, tmp.begin());
    return tmp;
}

std::vector<uint> ocp_qp::nbu() const {
    std::vector<uint> tmp(N+1);
    std::copy_n(qp->dim->nbu, N+1, tmp.begin());
    return tmp;
}

std::vector<uint> ocp_qp::ng() const {
    std::vector<uint> tmp(N+1);
    std::copy_n(qp->dim->ng, N+1, tmp.begin());
    return tmp;
}

std::ostream& operator<<(std::ostream& oss, const ocp_qp& qp) {
    static char a[10000];
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

std::pair<uint, uint> ocp_qp::dimensions(std::string field, uint stage) {

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
        return std::make_pair(num_elems_lbx(stage, qp->dim), 1);
    else if (field == "ubx")
        return std::make_pair(num_elems_ubx(stage, qp->dim), 1);
    else if (field == "lbu")
        return std::make_pair(num_elems_lbu(stage, qp->dim), 1);
    else if (field == "ubu")
        return std::make_pair(num_elems_ubu(stage, qp->dim), 1);
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

void ocp_qp::check_nb_elements(std::string field, uint stage, uint nb_elems) {
    if (!match(dimensions(field, stage), nb_elems))
        throw std::invalid_argument("Need " + std::to_string(dimensions(field, stage)) + " elements.");
}

ocp_qp_solver_plan string_to_plan(string solver) {

    ocp_qp_solver_plan plan;

    if (solver == "condensing_hpipm") {
        plan.qp_solver = FULL_CONDENSING_HPIPM;
    } else if (solver == "sparse_hpipm") {
        plan.qp_solver = PARTIAL_CONDENSING_HPIPM;
    } else if (solver == "hpmpc") {
        plan.qp_solver = PARTIAL_CONDENSING_HPMPC;
    } else if (solver == "ooqp") {
        plan.qp_solver = PARTIAL_CONDENSING_OOQP;
    } else if (solver == "qpdunes") {
        plan.qp_solver = PARTIAL_CONDENSING_QPDUNES;
    } else if (solver == "qpoases") {
        plan.qp_solver = FULL_CONDENSING_QPOASES;
    } else if (solver == "qore") {
        plan.qp_solver = FULL_CONDENSING_QORE;
    } else {
        throw std::invalid_argument("Solver not known.");
    }
    return plan;
}

}  // namespace acados
