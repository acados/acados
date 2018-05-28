
#ifndef INTERFACES_ACADOS_CPP_OCP_NLP_OCP_NLP_SOLUTION_HPP_
#define INTERFACES_ACADOS_CPP_OCP_NLP_OCP_NLP_SOLUTION_HPP_

#include <memory>
#include <vector>

#include "acados/ocp_nlp/ocp_nlp_common.h"

namespace acados
{
struct ocp_nlp_info
{
    int num_iter;
    int status;
    double inf_norm_residual;
    double total_time;
    int num_qp_iter;
};

class ocp_nlp_solution
{
 public:
    ocp_nlp_solution(std::shared_ptr<ocp_nlp_out> solution, std::shared_ptr<ocp_nlp_dims> dims,
                     int status);

    ocp_nlp_solution(const ocp_nlp_solution& other);

    std::vector<std::vector<double>> states();
    std::vector<std::vector<double>> controls();
    std::vector<std::vector<double>> lag_mul_dynamics();
    std::vector<std::vector<double>> lag_mul_bounds();
    std::vector<std::vector<double>> lag_mul_constraints();

    ocp_nlp_info info();

    const int N;

 private:
    std::shared_ptr<ocp_nlp_out> nlp_out_;

    std::shared_ptr<ocp_nlp_dims> dims_;

    int status_;
};
}  // namespace acados

#endif  // INTERFACES_ACADOS_CPP_OCP_NLP_OCP_NLP_SOLUTION_HPP_
