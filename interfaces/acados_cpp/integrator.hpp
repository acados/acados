
#ifndef INTERFACES_ACADOS_CPP_INTEGRATOR_HPP_
#define INTERFACES_ACADOS_CPP_INTEGRATOR_HPP_

#include <map>
#include <string>
#include <vector>

#include "acados/sim/sim_common.h"

namespace acados
{
class integrator
{
 public:
    integrator(const casadi::Function &model_fun, std::map<std::string, option_t *> options = {});

    std::map<std::string, option_t *> integrate( std::map<std::string, option_t *> in );

    virtual int num_stages() = 0;

    virtual ~integrator() = 0;

private:
    sim_solver_config* _config;
    sim_rk_opts* _opts;
    sim_in* _in;
    sim_out* _out;
    sim_solver* _solver
};

}  // namespace acados

#endif  // INTERFACES_ACADOS_CPP_INTEGRATOR_HPP_
