
#include "acados_cpp/ocp_nlp/ocp_nlp_solution.hpp"

namespace acados {

ocp_nlp_solution::ocp_nlp_solution(std::unique_ptr<ocp_nlp_out> solution)
    : N(0), nlp_out_(nullptr) {

    if (solution == nullptr)
        throw std::invalid_argument("Null pointer passed to constructor");

    nlp_out_.reset(solution);
}

}  // namespace acados
