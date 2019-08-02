
#ifndef INTERFACES_ACADOS_CPP_OCP_NLP_OCP_NLP_HPP_
#define INTERFACES_ACADOS_CPP_OCP_NLP_OCP_NLP_HPP_

#include <map>
#include <string>
#include <vector>

#include "acados/utils/types.h"
#include "acados_c/ocp_nlp_interface.h"

#include "acados_cpp/ocp.hpp"
#include "acados_cpp/ocp_nlp/ocp_nlp_solution.hpp"
#include "acados_cpp/ocp_nlp/casadi_module.hpp"
#include "acados_cpp/options.hpp"

namespace acados
{

class code_generator;

class ocp_nlp : private ocp
{
 public:
    ocp_nlp(std::vector<int> nx, std::vector<int> nu, std::vector<int> ng, std::vector<int> nh,
            std::vector<int> ns);

    ocp_nlp(int N, int nx, int nu, int ng = 0, int nh = 0, int ns = 0);

    ~ocp_nlp();

    void initialize_solver(std::string solver_name, std::map<std::string, option_t *> options = {});

    ocp_nlp_solution solve(std::vector<double> x_guess = {}, std::vector<double> u_guess = {});

    std::vector<std::string> fields = {"lbx", "ubx", "lbu", "ubu", "C",
                                       "D",   "lg",  "ug",  "lh",  "uh"};

    void set_field(std::string field, std::vector<double> v);

    void set_field(std::string field, int stage, std::vector<double> v);

    void set_dynamics(const casadi::Function &f, std::map<std::string, option_t *> options = {});

    // Affine path constraints
    // void set_stage_constraint(std::vector<double> C, std::vector<double> D);
    // void set_stage_constraint(int stage, std::vector<double> C, std::vector<double> D);
    // void set_terminal_constraint(std::vector<double> C, std::vector<double> D);

    // void set_stage_constraint(int stage, const casadi::Function& h);
    // void set_terminal_constraint(const casadi::Function& h);

    // Linear least squares cost
    void set_stage_cost(std::vector<double> C, std::vector<double> y_ref, std::vector<double> W);
    void set_stage_cost(int stage, std::vector<double> C, std::vector<double> y_ref,
                        std::vector<double> W);
    void set_terminal_cost(std::vector<double> C, std::vector<double> y_ref, std::vector<double> W);

    // Nonlinear least squares cost
    void set_stage_cost(const casadi::Function& r, std::vector<double> y_ref,
                        std::vector<double> W);
    void set_stage_cost(int stage, const casadi::Function& r, std::vector<double> y_ref,
                        std::vector<double> W);
    void set_terminal_cost(const casadi::Function& r, std::vector<double> y_ref,
                           std::vector<double> W);

    void generate_s_function(std::string filename);

    const int N;

 private:
    void squeeze_dimensions(std::map<std::string, std::vector<std::vector<double>>> bounds) override
    {
        ocp::squeeze_dimensions(bounds);
    }

    void fill_bounds(std::map<std::string, std::vector<std::vector<double>>> bounds) override
    {
        ocp::fill_bounds(bounds);
    }

    void change_bound_dimensions(std::vector<int> nbx, std::vector<int> nbu) override;

    bool needs_initializing() override;

    void needs_initializing(bool) override;

    void set_bound(std::string bound, int stage, std::vector<double> new_bound) override;

    std::vector<int> get_bound_indices(std::string, int) override;

    void set_bound_indices(std::string, int, std::vector<int>) override;

    int num_stages() override;

    std::unique_ptr<ocp_nlp_in> nlp_;

    std::shared_ptr<ocp_nlp_dims> dims_;

    std::unique_ptr<ocp_nlp_config> config_;

    std::unique_ptr<ocp_nlp_solver> solver_;

    std::unique_ptr<void, void (*)(void *)> solver_options_{nullptr, &std::free};

    std::map<std::string, casadi_module> module_;

    std::map<std::string, std::vector<int>> d_;

    std::map<std::string, std::vector<std::vector<double>>> cached_bounds;

    bool needs_initializing_;

    std::shared_ptr<ocp_nlp_out> result_;

    ocp_nlp_plan *plan_;

    std::string cached_model_;

    friend code_generator;
};

}  // namespace acados

#endif  // INTERFACES_ACADOS_CPP_OCP_NLP_OCP_NLP_HPP_
