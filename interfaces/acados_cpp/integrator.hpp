#ifndef INTERFACES_ACADOS_CPP_INTEGRATOR_HPP_
#define INTERFACES_ACADOS_CPP_INTEGRATOR_HPP_

#include <map>
#include <string>
#include <vector>

#include "acados_c/sim_interface.h"

#include "acados_cpp/ocp_nlp/casadi_module.hpp"
#include "acados_cpp/options.hpp"

namespace acados
{
typedef enum
{
    EXPLICIT = 0,
    IMPLICIT,
    GNSF_MODEL
} model_t;


class integrator
{
 public:
    explicit integrator
        (const casadi::Function& model, std::map<std::string, option_t*> options = {});

    ~integrator();

    std::vector<double> integrate(std::vector<double> x, std::vector<double> u = {}) const;

    void print_settings() const;

    void set_step_size(const double step_size);

    int num_stages() const { return opts_->ns; }
    double step_size() const { return in_->T; }


 private:
    void model_set(casadi::Function& model, std::map<std::string, option_t*> options = {});

    casadi::Function explicit2implicit(const casadi::Function &model);

    sim_config* config_;
    sim_opts* opts_;
    sim_in* in_;
    sim_out* out_;
    sim_solver* solver_;
    sim_solver_plan sim_plan_;

    void* dims_;

    size_t nx_;
    size_t nu_;
    size_t nz_;

    model_t model_type_;
    bool use_MX_;

    // std::map<std::string, casadi_module> module_;
};

}  // namespace acados

#endif  // INTERFACES_ACADOS_CPP_INTEGRATOR_HPP_
