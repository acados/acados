
#ifndef INTERFACES_ACADOS_CPP_CODE_GENERATOR_HPP_
#define INTERFACES_ACADOS_CPP_CODE_GENERATOR_HPP_

#include <fstream>
#include <string>

#include "acados_cpp/ocp_nlp/ocp_nlp.hpp"

namespace acados
{

class code_generator
{

 public:

    explicit code_generator(ocp_nlp *nlp);

    void generate_s_function(std::string name);

 private:

    void generate_s_function_header(std::ostream& out, std::string s_function_name);

    void generate_mdl_initialize_sizes(std::ostream& out);

    void generate_mdl_initialize_sample_times(std::ostream& out);

    void generate_mdl_start(std::ostream& out);

    void generate_mdl_outputs(std::ostream& out);

    void generate_mdl_terminate(std::ostream& out);

    void generate_s_function_footer(std::ostream& out);

    void generate_s_function_makefile(std::ostream& out, std::string name);

    void generate_dspace_makefile(std::ostream& out, std::string name);

    ocp_nlp *nlp_;
};

}  // namespace acados

#endif  // INTERFACES_ACADOS_CPP_CODE_GENERATOR_HPP_
