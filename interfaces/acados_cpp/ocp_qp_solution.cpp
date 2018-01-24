
#include "interfaces/acados_cpp/ocp_qp_solution.hpp"

#include "hpipm/include/hpipm_d_ocp_qp_sol.h"

namespace acados {

ocp_qp_solution::ocp_qp_solution(ocp_qp_out *solution) : N(solution->dim->N), qp_out(solution) {}

vector<vector<double>> ocp_qp_solution::states() {
    vector<vector<double>> result;
    for (int stage = 0; stage <= N; ++stage) {
        vector<double> tmp(qp_out->dim->nx[stage]);
        d_cvt_ocp_qp_sol_to_colmaj_x(qp_out, tmp.data(), stage);
        result.push_back(tmp);
    }
    return result;
}
    // vector<vector<double>> controls();
    // vector<vector<double>> eq_multipliers();
    // vector<vector<double>> ineq_multipliers();
ocp_qp_info ocp_qp_solution::info() {
    return *((ocp_qp_info *) qp_out->misc);
}

}  // namespace acados
