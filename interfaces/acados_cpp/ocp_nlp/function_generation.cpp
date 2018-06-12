
#include "acados_cpp/ocp_nlp/function_generation.hpp"

#include <vector>

namespace acados
{

using std::string;
using std::vector;

casadi_module generate_nls_residual(const casadi::Function& residual, string output_folder)
{
    casadi::SX x = residual.sx_in(0);
    casadi::SX u = residual.sx_in(1);

    vector<casadi::SX> xu {x, u};
    vector<casadi::SX> ux {u, x};

    casadi::SX r_new = casadi::SX::vertcat(residual(xu));
    casadi::SX r_jacT = casadi::SX::jacobian(r_new, casadi::SX::vertcat(ux)).T();

    casadi::Function r_fun(residual.name() + "_nls_res", {casadi::SX::vertcat(ux)},
                                                         {r_new, r_jacT});

    return casadi_module(r_fun, output_folder);
}

casadi_module generate_forward_vde(const casadi::Function& model, string output_folder)
{
    casadi::SX x = model.sx_in(0);
    casadi::SX u = model.sx_in(1);

    int_t nx = x.size1();
    int_t nu = u.size1();

    casadi::SX rhs = casadi::SX::vertcat(model(vector<casadi::SX>({x, u})));

    casadi::SX Sx = casadi::SX::sym("Sx", nx, nx);
    casadi::SX Su = casadi::SX::sym("Su", nx, nu);

    casadi::SX vde_x = casadi::SX::jtimes(rhs, x, Sx);
    casadi::SX vde_u = casadi::SX::jacobian(rhs, u) + casadi::SX::jtimes(rhs, x, Su);

    casadi::Function vde_fun(model.name() + "_expl_vde_for", {x, Sx, Su, u}, {rhs, vde_x, vde_u});

    return casadi_module(vde_fun, output_folder);
}

casadi_module generate_ode_jacobian(const casadi::Function& model, string output_folder)
{
    casadi::SX x = model.sx_in(0);
    casadi::SX u = model.sx_in(1);

    int_t nx = x.size1();
    int_t nu = u.size1();

    casadi::SX rhs = casadi::SX::vertcat(model(vector<casadi::SX>({x, u})));

    casadi::SX Sx = casadi::SX::sym("Sx", nx, nx);
    casadi::SX Su = casadi::SX::sym("Su", nx, nu);

    casadi::SX jac_x = casadi::SX::jacobian(rhs, x);
    casadi::Function jac_fun(model.name() + "_expl_ode_jac", {x, Sx, Su, u}, {rhs, jac_x});

    return casadi_module(jac_fun, output_folder);
}

}  // namespace acados
