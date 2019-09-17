
#ifndef INTERFACES_ACADOS_CPP_OCP_NLP_CASADI_MODULE_HPP_
#define INTERFACES_ACADOS_CPP_OCP_NLP_CASADI_MODULE_HPP_

#include <memory>
#include <string>

#include "acados_c/external_function_interface.h"

#include "casadi/casadi.hpp"

namespace acados
{

class casadi_module {

 public:

    casadi_module();

    casadi_module(const casadi::Function& function, std::string output_folder);

    casadi_module(const casadi_module& other) = delete;

    casadi_module& operator=(const casadi_module& other) = delete;

    casadi_module(casadi_module&& other);

    casadi_module& operator=(casadi_module&& other);

    ~casadi_module();

    const external_function_casadi *as_external_function() const;

    std::string path_to_header();

    std::string name();

 private:

    void load_functions(std::string output_folder);

    void generate(std::string output_folder);

    casadi::Function function_;

    external_function_casadi external_function_;

    std::unique_ptr<void, void (*)(void *)> handle_;

    std::string generated_header_;

};

}  // namespace acados

#endif  // INTERFACES_ACADOS_CPP_OCP_NLP_CASADI_MODULE_HPP_
