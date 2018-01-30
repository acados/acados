
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

vector<vector<double>> ocp_qp_solution::controls() {
    vector<vector<double>> result;
    for (int stage = 0; stage <= N; ++stage) {
        vector<double> tmp(qp_out->dim->nu[stage]);
        d_cvt_ocp_qp_sol_to_colmaj_u(qp_out, tmp.data(), stage);
        result.push_back(tmp);
    }
    return result;    
}

vector<vector<double>> ocp_qp_solution::lag_mul_dynamics() {
    vector<vector<double>> result;
    for (int stage = 0; stage <= N-1; ++stage) {
        vector<double> tmp(qp_out->dim->nx[stage+1]);
        d_cvt_ocp_qp_sol_to_colmaj_pi(qp_out, tmp.data(), stage);
        result.push_back(tmp);
    }
    return result;
}

vector<vector<double>> ocp_qp_solution::lag_mul_bounds() {
    vector<vector<double>> result;
    for (int stage = 0; stage <= N; ++stage) {
        int nb = qp_out->dim->nb[stage];
        vector<double> tmp(2 * nb);
        d_cvt_ocp_qp_sol_to_colmaj_lam_lb(qp_out, tmp.data(), stage);
        d_cvt_ocp_qp_sol_to_colmaj_lam_ub(qp_out, tmp.data() + nb, stage);
        result.push_back(tmp);
    }
    return result;    
}

vector<vector<double>> ocp_qp_solution::lag_mul_constraints() {
    vector<vector<double>> result;
    for (int stage = 0; stage <= N-1; ++stage) {
        int ng = qp_out->dim->ng[stage];
        vector<double> tmp(2 * ng);
        d_cvt_ocp_qp_sol_to_colmaj_lam_lg(qp_out, tmp.data(), stage);
        d_cvt_ocp_qp_sol_to_colmaj_lam_ug(qp_out, tmp.data() + ng, stage);
        result.push_back(tmp);
    }
    return result;    
}

ocp_qp_info ocp_qp_solution::info() {
    return *((ocp_qp_info *) qp_out->misc);
}

}  // namespace acados
