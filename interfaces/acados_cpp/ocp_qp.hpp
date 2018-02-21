
#ifndef ACADOS_INTERFACES_ACADOS_CPP_OCP_QP_HPP_
#define ACADOS_INTERFACES_ACADOS_CPP_OCP_QP_HPP_

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <utility>

#include "acados_c/ocp_qp.h"
#include "acados_cpp/ocp_qp_solution.hpp"
#include "acados_cpp/options.hpp"

namespace acados {

class ocp_qp {

public:

    ocp_qp(std::vector<uint> nx, std::vector<uint> nu, std::vector<uint> nbx, std::vector<uint> nbu, std::vector<uint> ng);

    ocp_qp(std::map<std::string, std::vector<uint>>);

    ocp_qp(uint N, uint nx, uint nu, uint nbx = 0, uint nbu = 0, uint ng = 0);

    void set(std::string field, uint stage, std::vector<double> v);
    void set(std::string field, std::vector<double> v);

    void initialize_solver(std::string solver_name, std::map<std::string, option_t *> options = {});

    ocp_qp_solution solve();

    std::vector< std::vector<double> > extract(std::string field);

    std::map<std::string, std::vector<uint>> dimensions();

    std::pair<uint, uint> dimensions(std::string field, uint stage);

    void set_bounds_indices(std::string name, uint stage, std::vector<uint> v);

    std::vector<std::vector<uint>> bounds_indices(std::string name);

    const uint N;

private:

    std::vector<uint> squeezed_bounds_indices(const std::string bound, uint stage);

    void squeeze_dimensions();

    void expand_dimensions();

    void check_range(std::string field, uint stage);
    
    void check_nb_elements(std::string, uint stage, uint nb_elems);

    void flatten(std::map<std::string, option_t *>& input, std::map<std::string, option_t *>& output);

    std::vector<uint> nx();
    std::vector<uint> nu();
    std::vector<uint> nbx();
    std::vector<uint> nbu();
    std::vector<uint> ng();

    std::unique_ptr<ocp_qp_in> qp;

    std::unique_ptr<ocp_qp_solver> solver;

    std::string cached_solver;

    bool needs_initializing = true;

    static std::map<std::string, std::function<void(int, ocp_qp_in *, double *)>> extract_functions;

    friend std::ostream& operator<<(std::ostream& oss, const ocp_qp& qp);

};


}  // namespace acados

#endif  // ACADOS_INTERFACES_ACADOS_CPP_OCP_QP_SOLUTION_HPP_
