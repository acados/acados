
#include "acados_cpp/ocp_nlp/ocp_nlp.hpp"

#include "casadi/casadi.hpp"

namespace acados {

void set_model(const casadi::Function& f, std::map<std::string, option_t *> options = {}) {

    casadi::Dict opts;
    opts["with_header"] = true;
    opts["with_export"] = false;

    f.generate(f.name() + ".c", opts);

};

}  // namespace acados
