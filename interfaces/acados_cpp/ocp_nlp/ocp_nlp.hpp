
#ifndef INTERFACES_ACADOS_CPP_OCP_NLP_OCP_NLP_H_
#define INTERFACES_ACADOS_CPP_OCP_NLP_OCP_NLP_H_

#include <map>
#include <string>
#include <vector>

#include "acados_c/ocp_nlp_interface.h"
#include "acados/utils/types.h"

#include "acados_cpp/ocp.hpp"
#include "acados_cpp/options.hpp"
#include "acados_cpp/ocp_nlp/ocp_nlp_solution.hpp"

#include "casadi/casadi.hpp"

namespace acados {

class ocp_nlp : private ocp {
 
 public:

    ocp_nlp(std::vector<int> nx, std::vector<int> nu, std::vector<int> ng, std::vector<int> nh,
            std::vector<int> ns);

    ocp_nlp(int N, int nx, int nu, int ng = 0, int nh = 0, int ns = 0);

    void initialize_solver(std::string solver_name, std::map<std::string, option_t *> options = {});

    ocp_nlp_solution solve();

    std::vector<std::string> fields = {"lbx", "ubx", "lbu", "ubu", "C", "D", "lg", "ug", "lh", "uh"};

    void set_field(std::string field, std::vector<double> v);

    void set_field(std::string field, int stage, std::vector<double> v);

    void set_dynamics(const casadi::Function& f, std::map<std::string, option_t *> options = {});

    // void set_stage_constraint(const casadi::Function& h);
    // void set_terminal_constraint(const casadi::Function& h_N);

    // void set_stage_cost(const casadi::Function& l);
    // void set_stage_cost(std::vector<double> C, std::vector<double> W);
    // void set_stage_cost(const casadi::Function& r, std::vector<double> W);

    // void set_terminal_cost(const casadi::Function& l_N);
    // void set_terminal_cost(std::vector<double> C_N, std::vector<double> W_N);
    // void set_terminal_cost(const casadi::Function& r_N, std::vector<double> W_N);

    const int N;

 private:

    void squeeze_dimensions(std::map<std::string, std::vector<std::vector<double>>> bounds) override {
        ocp::squeeze_dimensions(bounds);
    }

    void fill_bounds(std::map<std::string, std::vector<std::vector<double>>> bounds) override {
        ocp::fill_bounds(bounds);
    }

    void change_bound_dimensions(std::vector<int> nbx, std::vector<int> nbu) override;

    bool needs_initializing() override;

    void needs_initializing(bool) override;

    void set_bound(std::string bound, int stage, std::vector<double> new_bound) override;

    std::vector<int> get_bound_indices(std::string, int) override;

    void set_bound_indices(std::string, int, std::vector<int>) override;

    int num_stages() override;

    std::map<std::string, casadi::Function> create_explicit_ode_functions(const casadi::Function& model);

    std::unique_ptr<ocp_nlp_in> nlp_;

    std::shared_ptr<ocp_nlp_dims> dims_;

    std::unique_ptr<ocp_nlp_solver_config> config_;

    std::unique_ptr<ocp_nlp_solver> solver_;

    std::unique_ptr<void, void (*)(void *)> solver_options_ {nullptr, &std::free};

    std::map<std::string, void *> dynamics_handle;

    external_function_casadi forw_vde_;

    std::map<std::string, std::vector<int>> dimensions_;

    std::map<std::string, std::vector<std::vector<double>>> cached_bounds;

    bool needs_initializing_;

};

void *load_function(std::string function_name, void *handle);

void *compile_and_load(std::string name);

}  // namespace acados

#endif  // INTERFACES_ACADOS_CPP_OCP_NLP_OCP_NLP_H_
