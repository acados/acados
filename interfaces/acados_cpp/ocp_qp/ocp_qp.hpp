
#ifndef ACADOS_INTERFACES_ACADOS_CPP_OCP_QP_HPP_
#define ACADOS_INTERFACES_ACADOS_CPP_OCP_QP_HPP_

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <utility>

#include "acados/ocp_qp/ocp_qp_common.h"

#include "acados_c/ocp_qp_interface.h"
#include "acados_cpp/ocp_qp/ocp_qp_solution.hpp"
#include "acados_cpp/options.hpp"

namespace acados {

class ocp_qp {

public:

    ocp_qp(std::vector<int> nx, std::vector<int> nu, std::vector<int> ng, std::vector<int> ns);

    ocp_qp(int N, int nx, int nu, int ng = 0, int ns = 0);

    std::vector< std::vector<double> > get_field(std::string field);
    
    void set_field(std::string field, int stage, std::vector<double> v);
    
    void set_field(std::string field, std::vector<double> v);
    
    std::pair<int, int> shape_of_field(std::string field, int stage);

    void initialize_solver(std::string solver_name, std::map<std::string, option_t *> options = {});

    ocp_qp_solution solve();

    std::map<std::string, std::vector<int>> dimensions();

    const std::vector<std::string> fields {
        "Q", "S", "R", "q", "r", "A", "B", "b", "lbx", "ubx", "lbu", "ubu", "C", "D", "lg", "ug"
    };

    const int N;

private:

    void reset_bounds();

    void fill_in_bounds();

    void squeeze_dimensions();

    void expand_dimensions();

    std::vector<std::vector<int>> get_bounds_indices(std::string name);

    void set_bounds_indices(std::string name, int stage, std::vector<int> v);

    bool in_range(std::string field, int stage);

    std::vector<int> nx();
    std::vector<int> nu();
    std::vector<int> nbx();
    std::vector<int> nbu();
    std::vector<int> ng();

    std::map<std::string, std::vector<std::vector<double>>> cached_bounds;

    std::unique_ptr<ocp_qp_in> qp;

    std::unique_ptr<ocp_qp_solver> solver;

    const std::map<std::string, ocp_qp_solver_plan> available_solvers = {
        {"condensing_hpipm", {FULL_CONDENSING_HPIPM}},
        {"sparse_hpipm", {PARTIAL_CONDENSING_HPIPM}},
#ifdef ACADOS_WITH_HPMPC
        {"hpmpc", {PARTIAL_CONDENSING_HPMPC}},
#endif
#ifdef ACADOS_WITH_OOQP
        {"ooqp", {PARTIAL_CONDENSING_OOQP}}
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

    std::unique_ptr<void, void (*)(void *)> args {nullptr, std::free};

    std::string cached_solver;

    bool needs_initializing;

    static std::map<std::string, std::function<void(int, ocp_qp_in *, double *)>> extract_functions;
};


}  // namespace acados

#endif  // ACADOS_INTERFACES_ACADOS_CPP_OCP_QP_SOLUTION_HPP_
