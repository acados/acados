
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

using std::vector;
using std::string;
using std::map;

namespace acados {

class ocp_qp {

public:

    ocp_qp(vector<uint> nx, vector<uint> nu, vector<uint> nbx, vector<uint> nbu, vector<uint> ng);

    ocp_qp(map<string, vector<uint>>);

    ocp_qp(uint N, uint nx, uint nu, uint nbx = 0, uint nbu = 0, uint ng = 0);

    void set(string field, uint stage, vector<double> v);
    void set(string field, vector<double> v);

    void initialize_solver(string solver_name, map<string, option_t *> options = {});

    ocp_qp_solution solve();

    vector< vector<double> > extract(string field);

    map<string, vector<uint>> dimensions();

    std::pair<uint, uint> dimensions(string field, uint stage);

    void bounds_indices(std::string name, uint stage, vector<uint> v);

    std::vector<std::vector<uint>> bounds_indices(string name);

    const uint N;

private:
    
    void squeeze_dimensions();

    void expand_dimensions();

    void check_range(string field, uint stage);
    
    void check_nb_elements(string, uint stage, uint nb_elems);

    void flatten(map<string, option_t *>& input, map<string, option_t *>& output);

    vector<uint> nx();
    vector<uint> nu();
    vector<uint> nbx();
    vector<uint> nbu();
    vector<uint> ng();

    std::unique_ptr<ocp_qp_in> qp;

    std::unique_ptr<ocp_qp_solver> solver;

    std::string cached_solver;

    static std::map<string, std::function<void(int, ocp_qp_in *, double *)>> extract_functions;

    friend std::ostream& operator<<(std::ostream& oss, const ocp_qp& qp);

};


}  // namespace acados

#endif  // ACADOS_INTERFACES_ACADOS_CPP_OCP_QP_SOLUTION_HPP_
