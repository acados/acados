
#ifndef ACADOS_INTERFACES_ACADOS_CPP_OCP_QP_DIMENSIONS_HPP_
#define ACADOS_INTERFACES_ACADOS_CPP_OCP_QP_DIMENSIONS_HPP_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "acados_c/ocp_qp_interface.h"

namespace acados {

bool
are_valid_ocp_qp_dimensions(std::map<std::string, std::vector<uint>> dims);

std::unique_ptr<ocp_qp_dims>
create_ocp_qp_dimensions_ptr(const std::map<std::string, std::vector<uint>>&);

}  // namespace acados

#endif  // ACADOS_INTERFACES_ACADOS_CPP_OCP_QP_DIMENSIONS_HPP_
