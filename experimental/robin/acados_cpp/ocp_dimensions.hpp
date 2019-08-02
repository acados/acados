
#ifndef INTERFACES_ACADOS_CPP_OCP_DIMENSIONS_HPP_
#define INTERFACES_ACADOS_CPP_OCP_DIMENSIONS_HPP_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "acados_c/ocp_nlp_interface.h"
#include "acados_c/ocp_qp_interface.h"

namespace acados
{
/*
 * OCP dimensions should have fields
 * nx  - array with the number of states, per stage
 * nu  - array ""       ""        controls
 *
 * and can have following optional fields
 * nbx - array ""       ""        state bounds
 * nbu - array ""       ""        control bounds
 * ng  - array ""       ""        polytopic constraints
 * nh  - array ""       ""        nonlinear constraints
 * ns  - array ""       ""        slack variables
 */
bool are_valid_ocp_dimensions(const std::map<std::string, std::vector<int>>& dims,
                              std::vector<std::string> valid_names);

std::unique_ptr<ocp_qp_dims> create_ocp_qp_dimensions_ptr(
    const std::map<std::string, std::vector<int>>&);

// std::unique_ptr<ocp_nlp_dims>
// create_ocp_nlp_dimensions_ptr(const std::map<std::string, std::vector<int>>&);

}  // namespace acados

#endif  // INTERFACES_ACADOS_CPP_OCP_DIMENSIONS_HPP_
