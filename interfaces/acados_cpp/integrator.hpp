#ifndef INTERFACES_ACADOS_CPP_INTEGRATOR_HPP_
#define INTERFACES_ACADOS_CPP_INTEGRATOR_HPP_

#include <map>
#include <string>
#include <vector>


#include "acados/utils/types.h"
#include "acados_c/sim_interface.h"


#include "acados_cpp/ocp_nlp/casadi_module.hpp"
#include "acados_cpp/options.hpp"

// TODO(oj): try to only use C interface, i.e. remove line above and use void pointers
// @tobi: or what do u think?
namespace acados
{
typedef enum
{
    // TODO do we really need this (see definition of sim_solver_t)
    EXPLICIT = 0,
    IMPLICIT,
    GNSF
} model_t;


class integrator
{
 public:
    integrator(const casadi::Function& model, std::map<std::string, option_t*> options = {});

    std::vector<double> integrate(std::vector<double> x, std::vector<double> u = {});


    void set_step(const double step);

    int num_stages() { return opts_->ns; }
    double step() { return in_->T; }

    ~integrator();

 private:
    void set_model(const casadi::Function& model, std::map<std::string, option_t*> options = {});

    sim_solver_config* config_;
    sim_rk_opts* opts_;
    sim_in* in_;
    sim_out* out_;
    sim_solver* solver_;
    void* dims_;
    size_t nx_;
    size_t nu_;
    size_t nz_;

    model_t model_type_;
    bool use_MX_;

    std::map<std::string, casadi_module> module_;
};

}  // namespace acados

#endif  // INTERFACES_ACADOS_CPP_INTEGRATOR_HPP_
