
#include <algorithm>
#include <iostream>
#include <iterator>
#include <string>
#include <stdexcept>

#include "acados_cpp/ocp_qp.hpp"

#include "acados_c/ocp_qp.h"
#include "acados/utils/print.h"
#include "blasfeo/include/blasfeo_d_aux.h"

namespace acados {

ocp_qp::ocp_qp(int N, std::vector<int> nx, std::vector<int> nu, std::vector<int> nbx, 
             std::vector<int> nbu, std::vector<int> ng) : N(N), qp(nullptr) {

    if (N <= 0) throw std::invalid_argument("Number of stages must be positive");

    ocp_qp_dims *dim = create_ocp_qp_dims(N);
    dim->N = N;

    // states
    write_dimensions(nx, dim->nx);

    // controls
    nu.at(N) = 0;
    write_dimensions(nu, dim->nu);

    // bounds
    write_dimensions(nbx, dim->nbx);
    write_dimensions(nbu, dim->nbu);

    std::vector<int> nb(N+1);
    for(int i = 0; i <= N; i++)
        nb.at(i) = nbx.at(i) + nbu.at(i);
    write_dimensions(nb, dim->nb);

    // constraints
    write_dimensions(ng, dim->ng);

    // slacks
    std::vector<int> ns(N+1, 0);
    write_dimensions(ns, dim->ns);

    qp = create_ocp_qp_in(dim);
}

ocp_qp::ocp_qp(int N, int nx, int nu, int nbx, int nbu, int ng)
    : ocp_qp(N, std::vector<int>(N+1, nx), std::vector<int>(N+1, nu), std::vector<int>(N+1, nbx),
      std::vector<int>(N+1, nbu), std::vector<int>(N+1, ng)) {}

void ocp_qp::update(std::string field, int stage, std::vector<double> v) {

    check_range(field, stage);
    check_nb_elements(field, stage, v.size());

    if (field == "Q")
        d_update_Q(stage, v.data(), qp);
    else if (field == "S")
        d_update_S(stage, v.data(), qp);
    else if (field == "R")
        d_update_R(stage, v.data(), qp);
    else if (field == "q")
        d_update_q(stage, v.data(), qp);
    else if (field == "r")
        d_update_r(stage, v.data(), qp);
    else if (field == "A")
        d_update_A(stage, v.data(), qp);
    else if (field == "B")
        d_update_B(stage, v.data(), qp);
    else if (field == "b")
        d_update_b(stage, v.data(), qp);
    else if (field == "lbx")
        d_update_lbx(stage, v.data(), qp);
    else if (field == "ubx")
        d_update_ubx(stage, v.data(), qp);
    else if (field == "lbu")
        d_update_lbu(stage, v.data(), qp);
    else if (field == "ubu")
        d_update_ubu(stage, v.data(), qp);
    else if (field == "C")
        d_update_C(stage, v.data(), qp);
    else if (field == "D")
        d_update_D(stage, v.data(), qp);
    else if (field == "lg")
        d_update_lg(stage, v.data(), qp);
    else if (field == "ug")
        d_update_ug(stage, v.data(), qp);
    else
        throw std::invalid_argument("OCP QP does not contain field " + field);
}

std::map<string, std::function<void(int, double *, ocp_qp_in *)>> ocp_qp::extract_functions = {
        {"Q", d_copy_Q},
        {"S", d_copy_S},
        {"R", d_copy_R},
        {"q", d_copy_q},
        {"r", d_copy_r},
        {"A", d_copy_A},
        {"B", d_copy_B},
        {"b", d_copy_b},
        {"lbx", d_copy_lbx},
        {"ubx", d_copy_ubx},
        {"lbu", d_copy_lbu},
        {"ubu", d_copy_ubu},
        {"C", d_copy_C},
        {"D", d_copy_D},
        {"lg", d_copy_lg},
        {"ug", d_copy_ug}
    };


vector< vector<double> > ocp_qp::extract(std::string field) {
    int last_index = N;
    if (field == "A" || field == "B" || field == "b")
        last_index = N-1;
    vector< vector<double> > result;
    for (int i = 0; i <= last_index; i++) {
        auto dims = dimensions(field, i);
        vector<double> v(dims.first * dims.second);
        extract_functions[field](i, v.data(), qp);
        result.push_back(v);
    }
    return result;
}


std::vector<int> ocp_qp::nx() {
    std::vector<int> tmp(N+1);
    std::copy_n(qp->dim->nx, N+1, tmp.begin());
    return tmp;
}

std::vector<int> ocp_qp::nu() {
    std::vector<int> tmp(N+1);
    std::copy_n(qp->dim->nu, N+1, tmp.begin());
    return tmp;
}

std::vector<int> ocp_qp::nbx() {
    std::vector<int> tmp(N+1);
    std::copy_n(qp->dim->nbx, N+1, tmp.begin());
    return tmp;
}

std::vector<int> ocp_qp::nbu() {
    std::vector<int> tmp(N+1);
    std::copy_n(qp->dim->nbu, N+1, tmp.begin());
    return tmp;
}

std::vector<int> ocp_qp::ng() {
    std::vector<int> tmp(N+1);
    std::copy_n(qp->dim->ng, N+1, tmp.begin());
    return tmp;
}

std::ostream& operator<<(std::ostream& oss, const ocp_qp& qp) {
    print_ocp_qp_in(qp.qp);
    return oss;
}

// std::vector<double> ocp_qp::getA(int i) {
//     int num_rows = dimensions->nx[i+1], num_cols = dimensions->nx[i];
//     double A_col_major[num_rows * num_cols];
//     blasfeo_unpack_tran_dmat(num_rows, num_cols, &(qp->BAbt[i]), dimensions->nu[i], 0,
//                              A_col_major, num_rows);
//     std::vector<double> A_raveled(num_rows * num_cols);
//     std::copy_n(A_col_major, num_rows * num_cols, A_raveled.begin());
//     return A_raveled;
// }

void ocp_qp::write_dimensions(std::vector<int> dims, int *ptr) {
    if (dims.size() != N+1)
        throw std::invalid_argument("Dimensions must be defined for all N+1 nodes.");
    if (std::any_of(dims.begin(), dims.end(), [](int a) {return a < 0;}))
        throw std::invalid_argument("Dimension must be non-negative for all N+1 nodes.");
    std::copy(dims.begin(), dims.end(), ptr);
}

void ocp_qp::check_range(std::string field, int stage) {
    int lower_bound = 0;
    int upper_bound;
    if (field == "A" || field == "B" || field == "b") {
        upper_bound = N-1;
    } else {
        upper_bound = N;
    }
    if(stage < lower_bound || stage > upper_bound)
        throw std::out_of_range(std::to_string(stage) + " must be in range [" +
              std::to_string(lower_bound) + ", " + std::to_string(upper_bound) + "].");
}

std::pair<int, int> ocp_qp::dimensions(std::string field, int stage) {

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

void ocp_qp::check_nb_elements(std::string field, int stage, int nb_elems) {
    auto dims = dimensions(field, stage);
    auto nb_expected_elems = dims.first * dims.second;
    if (nb_elems != nb_expected_elems)
        throw std::invalid_argument("Vector needs " + std::to_string(nb_expected_elems)
                                    + " elements.");
}

}  // namespace acados
