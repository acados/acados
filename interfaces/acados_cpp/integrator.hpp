
#ifndef INTERFACES_ACADOS_CPP_INTEGRATOR_HPP_
#define INTERFACES_ACADOS_CPP_INTEGRATOR_HPP_

#include <map>
#include <string>
#include <vector>

#include "acados/utils/types.h"
#include "acados_c/ocp_nlp_interface.h"

#include "acados_cpp/ocp.hpp"
#include "acados_cpp/ocp_nlp/ocp_nlp_solution.hpp"
#include "acados_cpp/ocp_nlp/casadi_module.hpp"
#include "acados_cpp/options.hpp"

// TODO(oj): try to only use C interface, i.e. remove line above and use void pointers
// @tobi: or what do u think?
namespace acados
{
class integrator
{
 public:
    integrator(const casadi::Function& model, std::map<std::string, option_t*> options = {});

    std::map<std::string, option_t*> integrate(std::vector<double> x, std::vector<double> u = {});

    virtual int num_stages() = 0;

    ~integrator();

 private:
    sim_solver_config* config_;
    sim_rk_opts* opts_;
    sim_in* in_;
    sim_out* out_;
    sim_solver* solver_;
    void* dims_;
    size_t nx_;
    size_t nu_;
    size_t nz_;
};

}  // namespace acados

#endif  // INTERFACES_ACADOS_CPP_INTEGRATOR_HPP_
