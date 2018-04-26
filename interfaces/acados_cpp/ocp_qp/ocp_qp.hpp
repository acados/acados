
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

    ocp_qp(std::vector<uint> nx, std::vector<uint> nu, std::vector<uint> nbx, std::vector<uint> nbu,
           std::vector<uint> ng, std::vector<uint> ns);

    ocp_qp(uint N, uint nx, uint nu, uint nbx = 0, uint nbu = 0, uint ng = 0, uint ns = 0);

    std::vector< std::vector<double> > get_field(std::string field);
    void set_field(std::string field, uint stage, std::vector<double> v);
    void set_field(std::string field, std::vector<double> v);
    std::pair<uint, uint> shape_of_field(std::string field, uint stage);

    void initialize_solver(std::string solver_name, std::map<std::string, option_t *> options = {});

    ocp_qp_solution solve();

    std::map<std::string, std::vector<uint>> dimensions();

    const std::vector<std::string> fields {
        "Q", "S", "R", "q", "r", "A", "B", "b", "lbx", "ubx", "lbu", "ubu", "C", "D", "lg", "ug"
    };

    const uint N;

private:

    vector<uint> idxb(vector<double> lower_bound, vector<double> upper_bound);

    void fill_in_bounds();

    void squeeze_dimensions();

    void expand_dimensions();

    std::vector<std::vector<uint>> bounds_indices(std::string name);

    void set_bounds_indices(std::string name, uint stage, std::vector<uint> v);

    bool in_range(std::string field, uint stage);

    std::vector<uint> nx();
    std::vector<uint> nu();
    std::vector<uint> nbx();
    std::vector<uint> nbu();
    std::vector<uint> ng();

    std::map<std::string, std::vector<std::vector<double>>> cached_bounds;

    std::unique_ptr<ocp_qp_in> qp;

    std::unique_ptr<ocp_qp_solver> solver;

    std::map<std::string, ocp_qp_solver_plan> available_solvers;

    std::unique_ptr<ocp_qp_xcond_solver_config> config;

    std::unique_ptr<void, void (*)(void *)> args {nullptr, std::free};

    std::string cached_solver;

    bool needs_initializing;

    static std::map<std::string, std::function<void(int, ocp_qp_in *, double *)>> extract_functions;
};


}  // namespace acados

#endif  // ACADOS_INTERFACES_ACADOS_CPP_OCP_QP_SOLUTION_HPP_
