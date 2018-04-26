
#include "acados_cpp/ocp_nlp/ocp_nlp.hpp"

#include "casadi/casadi.hpp"

namespace acados {

void ocp_nlp::set_model(const casadi::Function& f, std::map<std::string, option_t *> options) {

    casadi::Dict opts {{"with_header", true}, {"with_export", false}};
    f.generate(f.name() + ".c", opts);

};

}  // namespace acados
