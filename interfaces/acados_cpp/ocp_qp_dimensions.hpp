
#ifndef ACADOS_INTERFACES_ACADOS_CPP_OCP_QP_DIMENSIONS_HPP_
#define ACADOS_INTERFACES_ACADOS_CPP_OCP_QP_DIMENSIONS_HPP_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "acados_c/ocp_qp_interface.h"

namespace acados {

const std::vector<std::string> dimension_keys {"nx", "nu", "nbx", "nbu", "ng"};

bool valid_dimensions(std::map<std::string, std::vector<uint>> dims);

std::unique_ptr<ocp_qp_dims> make_dimensions_ptr(std::map<std::string, std::vector<uint>> dims);

}  // namespace acados

#endif  // ACADOS_INTERFACES_ACADOS_CPP_OCP_QP_DIMENSIONS_HPP_
