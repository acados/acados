
#ifndef INTERFACES_ACADOS_CPP_OCP_NLP_OCP_NLP_SOLUTION_HPP_
#define INTERFACES_ACADOS_CPP_OCP_NLP_OCP_NLP_SOLUTION_HPP_

#include <memory>

#include "acados/ocp_nlp/ocp_nlp_common.h"

namespace acados {

struct ocp_nlp_info {
    int num_iter;
    double inf_norm_residual;
};

class ocp_nlp_solution {

public:

    ocp_nlp_solution(std::unique_ptr<ocp_nlp_out> solution);

    ocp_nlp_solution(ocp_nlp_solution&& other);

    ocp_nlp_solution(const ocp_nlp_solution& other);

    std::vector<std::vector<double>> states();
    std::vector<std::vector<double>> controls();
    std::vector<std::vector<double>> lag_mul_dynamics();
    std::vector<std::vector<double>> lag_mul_bounds();
    std::vector<std::vector<double>> lag_mul_constraints();
    ocp_nlp_info info();

    const int N;

private:

    std::unique_ptr<ocp_nlp_out> nlp_out_;

};
}  // namespace acados

#endif  // INTERFACES_ACADOS_CPP_OCP_NLP_OCP_NLP_SOLUTION_HPP_
