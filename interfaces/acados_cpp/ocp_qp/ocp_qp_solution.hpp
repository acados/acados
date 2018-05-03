
#ifndef ACADOS_INTERFACES_ACADOS_CPP_OCP_QP_SOLUTION_HPP_
#define ACADOS_INTERFACES_ACADOS_CPP_OCP_QP_SOLUTION_HPP_

#include <memory>
#include <vector>

#include "acados/ocp_qp/ocp_qp_common.h"

using std::vector;

namespace acados {

class ocp_qp_solution {

public:

    ocp_qp_solution(std::unique_ptr<ocp_qp_out> solution);

    ocp_qp_solution(ocp_qp_solution&& other);

    ocp_qp_solution(const ocp_qp_solution& other);

    vector<vector<double>> states();
    vector<vector<double>> controls();
    vector<vector<double>> lag_mul_dynamics();
    vector<vector<double>> lag_mul_bounds();
    vector<vector<double>> lag_mul_constraints();
    ocp_qp_info info();

    const int N;

private:

    std::unique_ptr<ocp_qp_out> qp_out;

};

}  // namespace acados

#endif  // ACADOS_INTERFACES_ACADOS_CPP_OCP_QP_SOLUTION_HPP_
