
#ifndef INTERFACES_ACADOS_CPP_OCP_NLP_FUNCTION_GENERATION_HPP_
#define INTERFACES_ACADOS_CPP_OCP_NLP_FUNCTION_GENERATION_HPP_

#include <string>

#include "acados_cpp/ocp_nlp/casadi_module.hpp"

#include "casadi/casadi.hpp"

namespace acados
{

casadi_module generate_forward_vde(const casadi::Function& model, std::string output_dir =
                                   "_casadi_generated");

casadi_module generate_ode_jacobian(const casadi::Function& model, std::string output_dir =
                                    "_casadi_generated");

casadi_module generate_nls_residual(const casadi::Function& residual, std::string output_dir =
                                "_casadi_generated");

}  // namespace acados

#endif  // INTERFACES_ACADOS_CPP_OCP_NLP_FUNCTION_GENERATION_HPP_
