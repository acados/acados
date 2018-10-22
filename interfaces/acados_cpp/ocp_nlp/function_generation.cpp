
#include "acados_cpp/ocp_nlp/function_generation.hpp"

#include <vector>

namespace acados
{
using std::string;
using std::vector;

casadi_module generate_nls_residual(const casadi::Function& residual, string output_folder,
                                    const bool use_MX)
{
    casadi::Function r_fun;
    if (use_MX == false)
    {
        casadi::SX x = residual.sx_in(0);
        casadi::SX u = residual.sx_in(1);

        vector<casadi::SX> xu{x, u};
        vector<casadi::SX> ux{u, x};

        casadi::SX r_new = casadi::SX::vertcat(residual(xu));
        casadi::SX r_jacT = casadi::SX::jacobian(r_new, casadi::SX::vertcat(ux)).T();

        r_fun = casadi::Function(residual.name() + "_nls_res", {casadi::SX::vertcat(ux)},
                                 {r_new, r_jacT});
    }
    else  // MX
    {
        casadi::MX x = residual.mx_in(0);
        casadi::MX u = residual.mx_in(1);

        vector<casadi::MX> xu{x, u};
        vector<casadi::MX> ux{u, x};

        casadi::MX r_new = casadi::MX::vertcat(residual(xu));
        casadi::MX r_jacT = casadi::MX::jacobian(r_new, casadi::MX::vertcat(ux)).T();

        r_fun = casadi::Function(residual.name() + "_nls_res", {casadi::MX::vertcat(ux)},
                                 {r_new, r_jacT});
    }
    return casadi_module(r_fun, output_folder);
}

casadi_module generate_forward_vde(const casadi::Function& model, string output_folder,
                                   const bool use_MX)
{
    casadi::Function vde_fun;
    if (use_MX == false)
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

        vde_fun =
            casadi::Function(model.name() + "_expl_vde_for", {x, Sx, Su, u}, {rhs, vde_x, vde_u});
    }
    else  // MX
    {
        casadi::MX x = model.mx_in(0);
        casadi::MX u = model.mx_in(1);

        int_t nx = x.size1();
        int_t nu = u.size1();

        casadi::MX rhs = casadi::MX::vertcat(model(vector<casadi::MX>({x, u})));

        casadi::MX Sx = casadi::MX::sym("Sx", nx, nx);
        casadi::MX Su = casadi::MX::sym("Su", nx, nu);

        casadi::MX vde_x = casadi::MX::jtimes(rhs, x, Sx);
        casadi::MX vde_u = casadi::MX::jacobian(rhs, u) + casadi::MX::jtimes(rhs, x, Su);

        vde_fun =
            casadi::Function(model.name() + "_expl_vde_for", {x, Sx, Su, u}, {rhs, vde_x, vde_u});
    }

    return casadi_module(vde_fun, output_folder);
}


casadi_module generate_expl_ode_fun(const casadi::Function& model, string output_folder,
                                    const bool use_MX)
{
    casadi::Function fun;
    if (use_MX == false)
    {
        casadi::SX x = model.sx_in(0);
        casadi::SX u = model.sx_in(1);

        // int_t nx = x.size1();
        // int_t nu = u.size1();

        casadi::SX rhs = casadi::SX::vertcat(model(vector<casadi::SX>({x, u})));

        // casadi::SX Sx = casadi::SX::sym("Sx", nx, nx);
        // casadi::SX Su = casadi::SX::sym("Su", nx, nu);

        // casadi::SX vde_x = casadi::SX::jtimes(rhs, x, Sx);
        // casadi::SX vde_u = casadi::SX::jacobian(rhs, u) + casadi::SX::jtimes(rhs, x, Su);

        fun = casadi::Function(model.name() + "_expl_ode_fun", {x, u}, {rhs});
    }
    else  // MX
    {
        casadi::MX x = model.mx_in(0);
        casadi::MX u = model.mx_in(1);

        casadi::MX rhs = casadi::MX::vertcat(model(vector<casadi::MX>({x, u})));

        fun = casadi::Function(model.name() + "_expl_ode_fun", {x, u}, {rhs});
    }

    return casadi_module(fun, output_folder);
}


// TODO(oj): remove this, it is not used in erk
casadi_module generate_ode_jacobian(const casadi::Function& model, string output_folder,
                                    const bool use_MX)
{
    casadi::Function jac_fun;
    if (use_MX == false)
    {
        casadi::SX x = model.sx_in(0);
        casadi::SX u = model.sx_in(1);

        int_t nx = x.size1();
        int_t nu = u.size1();

        casadi::SX rhs = casadi::SX::vertcat(model(vector<casadi::SX>({x, u})));

        casadi::SX Sx = casadi::SX::sym("Sx", nx, nx);
        casadi::SX Su = casadi::SX::sym("Su", nx, nu);

        casadi::SX jac_x = casadi::SX::jacobian(rhs, x);
        jac_fun = casadi::Function(model.name() + "_expl_ode_jac", {x, Sx, Su, u}, {rhs, jac_x});
    }
    else  // MX
    {
        casadi::MX x = model.mx_in(0);
        casadi::MX u = model.mx_in(1);

        int_t nx = x.size1();
        int_t nu = u.size1();

        casadi::MX rhs = casadi::MX::vertcat(model(vector<casadi::MX>({x, u})));

        casadi::MX Sx = casadi::MX::sym("Sx", nx, nx);
        casadi::MX Su = casadi::MX::sym("Su", nx, nu);

        casadi::MX jac_x = casadi::MX::jacobian(rhs, x);
        jac_fun = casadi::Function(model.name() + "_expl_ode_jac", {x, Sx, Su, u}, {rhs, jac_x});
    }

    return casadi_module(jac_fun, output_folder);
}

}  // namespace acados
