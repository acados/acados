
#include <vector>

#include "acados/ocp_qp/ocp_qp_common.h"

using std::vector;

namespace acados {

class ocp_qp_solution {

public:

    ocp_qp_solution(ocp_qp_out *solution);

    vector<vector<double>> states();
    vector<vector<double>> controls();
    vector<vector<double>> lag_mul_dynamics();
    vector<vector<double>> lag_mul_bounds();
    vector<vector<double>> lag_mul_constraints();
    ocp_qp_info info();

    const int N;

    ocp_qp_out *qp_out;

};

}  // namespace acados
