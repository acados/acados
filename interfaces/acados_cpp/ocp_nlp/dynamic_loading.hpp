
#ifndef INTERFACES_ACADOS_CPP_OCP_NLP_DYNAMIC_LOADING_HPP_
#define INTERFACES_ACADOS_CPP_OCP_NLP_DYNAMIC_LOADING_HPP_

#include <string>

namespace acados {

#define TMP_DIR "_casadi_gen"

void *compile_and_load_library(std::string source_name);

void *load_function(std::string function_name, void *handle);

void free_handle(void *handle);

std::string load_error_message();

void create_directory(std::string dir_name);

}  // namespace acados

#endif  // INTERFACES_ACADOS_CPP_OCP_NLP_DYNAMIC_LOADING_HPP_
