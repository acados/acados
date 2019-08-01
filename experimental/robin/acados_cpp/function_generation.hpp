
#ifndef INTERFACES_ACADOS_CPP_FUNCTION_GENERATION_HPP_
#define INTERFACES_ACADOS_CPP_FUNCTION_GENERATION_HPP_

#include <string>

#include "acados_cpp/ocp_nlp/casadi_module.hpp"

#include "casadi/casadi.hpp"

namespace acados
{
/* IMPLICIT MODEL */
casadi_module generate_impl_ode_fun_jac_x_xdot_z(const casadi::Function& model,
                                std::string output_dir = "_autogen", const bool use_MX = false);

casadi_module generate_impl_ode_fun(const casadi::Function& model,
                                std::string output_dir = "_autogen", const bool use_MX = false);

casadi_module generate_impl_ode_fun_jac_x_xdot_u(const casadi::Function& model,
                                std::string output_dir = "_autogen", const bool use_MX = false);

casadi_module generate_impl_ode_jac_x_xdot_u_z(const casadi::Function& model,
                                std::string output_dir = "_autogen", const bool use_MX = false);

casadi_module generate_impl_ode_hess(const casadi::Function& model,
                                std::string output_dir = "_autogen", const bool use_MX = false);


/* EXPLICIT MODEL */
casadi_module generate_forward_vde(const casadi::Function& model,
                                std::string output_dir = "_autogen", const bool use_MX = false);

casadi_module generate_expl_ode_fun(const casadi::Function& model,
                                std::string output_dir = "_autogen", const bool use_MX = false);

casadi_module generate_expl_vde_adj(const casadi::Function& model,
                                std::string output_dir = "_autogen", const bool use_MX = false);

casadi_module generate_expl_ode_hess(const casadi::Function& model,
                                std::string output_dir = "_autogen", const bool use_MX = false);

/* NLP */
casadi_module generate_nls_residual(const casadi::Function& residual,
                                std::string output_dir = "_autogen", const bool use_MX = false);

}  // namespace acados

#endif  // INTERFACES_ACADOS_CPP_FUNCTION_GENERATION_HPP_
