
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

// TODO(oj): move the following functions one directory up
/* IMPLICIT MODEL */
casadi_module generate_impl_ode_fun(const casadi::Function& model, string output_folder,
                                    const bool use_MX)
{
    casadi::Function fun;
    if (use_MX == false)
    {
        casadi::SX x = model.sx_in("x");
        casadi::SX xdot = model.sx_in("xdot");
        casadi::SX u = model.sx_in("u");
        casadi::SX z = model.sx_in("z");

        casadi::SXDict arg_in =  {{"x", x}, {"xdot", xdot}, {"u", u}, {"z", z}};

        casadi::SXDict rhs_dict = (model(arg_in));
        casadi::SX rhs = rhs_dict.begin()->second;

        fun = casadi::Function(model.name() + "_impl_ode_fun", {x, xdot, u, z}, {rhs});
    }
    else    // MX
    {
        casadi::MX x = model.mx_in("x");
        casadi::MX xdot = model.mx_in("xdot");
        casadi::MX u = model.mx_in("u");
        casadi::MX z = model.mx_in("z");

        casadi::MXDict arg_in =  {{"x", x}, {"xdot", xdot}, {"u", u}, {"z", z}};

        casadi::MXDict rhs_dict = (model(arg_in));
        casadi::MX rhs = rhs_dict.begin()->second;

        fun = casadi::Function(model.name() + "_impl_ode_fun", {x, xdot, u, z}, {rhs});
    }
    return casadi_module(fun, output_folder);
}

casadi_module generate_impl_ode_fun_jac_x_xdot_z(const casadi::Function& model,
                        string output_folder, const bool use_MX)
{
    casadi::Function fun;
    if (use_MX == false)
    {
        casadi::SX x = model.sx_in("x");
        casadi::SX xdot = model.sx_in("xdot");
        casadi::SX u = model.sx_in("u");
        casadi::SX z = model.sx_in("z");
        // casadi::SX w = casadi::SX::vertcat(vector<casadi::SX>({x, u}));
        std::cout << "GENERATE IMPL MODEL- Got xuz" << std::endl;
        casadi::SXDict arg_in =  {{"x", x}, {"xdot", xdot}, {"u", u}, {"z", z}};

        casadi::SXDict rhs_dict = (model(arg_in));
        casadi::SX rhs = rhs_dict.begin()->second;

        casadi::SX jac_x = casadi::SX::jacobian(rhs, x);
        // casadi::SX jac_u = casadi::SX::jacobian(rhs, u);
        casadi::SX jac_z = casadi::SX::jacobian(rhs, z);
        casadi::SX jac_xdot = casadi::SX::jacobian(rhs, xdot);

        fun = casadi::Function(model.name() + "_impl_ode_fun_jac_x_xdot_z",
                            {x, xdot, u, z}, {rhs, jac_x, jac_xdot, jac_z});
    }
    else  // MX
    {
        casadi::MX x = model.mx_in("x");
        casadi::MX xdot = model.mx_in("xdot");
        casadi::MX u = model.mx_in("u");
        casadi::MX z = model.mx_in("z");
        // casadi::MX w = casadi::MX::vertcat(vector<casadi::MX>({x, u}));
        std::cout << "GENERATE IMPL MODEL- Got xuz" << std::endl;
        casadi::MXDict arg_in =  {{"x", x}, {"xdot", xdot}, {"u", u}, {"z", z}};

        casadi::MXDict rhs_dict = (model(arg_in));
        casadi::MX rhs = rhs_dict.begin()->second;

        casadi::MX jac_x = casadi::MX::jacobian(rhs, x);
        // casadi::MX jac_u = casadi::MX::jacobian(rhs, u);
        casadi::MX jac_z = casadi::MX::jacobian(rhs, z);
        casadi::MX jac_xdot = casadi::MX::jacobian(rhs, xdot);

        fun = casadi::Function(model.name() + "_impl_ode_fun_jac_x_xdot_z",
                            {x, xdot, u, z}, {rhs, jac_x, jac_xdot, jac_z});
    }
    return casadi_module(fun, output_folder);
}

/* EXPLICIT MODEL */
casadi_module generate_forward_vde(const casadi::Function& model, string output_folder,
                                   const bool use_MX)
{
    casadi::Function vde_fun;
    if (use_MX == false)
    {
        casadi::SX x = model.sx_in("x");
        casadi::SX u = model.sx_in("u");

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
        casadi::MX x = model.mx_in("x");
        casadi::MX u = model.mx_in("u");

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
        casadi::SX x = model.sx_in("x");
        casadi::SX u = model.sx_in("u");

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
        casadi::MX x = model.mx_in("x");
        casadi::MX u = model.mx_in("u");

        casadi::MX rhs = casadi::MX::vertcat(model(vector<casadi::MX>({x, u})));

        fun = casadi::Function(model.name() + "_expl_ode_fun", {x, u}, {rhs});
    }

    return casadi_module(fun, output_folder);
}


casadi_module generate_expl_vde_adj(const casadi::Function& model, string output_folder,
                                    const bool use_MX)
{
    casadi::Function fun;
    if (use_MX == false)  // SX
    {
        casadi::SX x = model.sx_in("x");
        casadi::SX u = model.sx_in("u");

        casadi::SX rhs = casadi::SX::vertcat(model(vector<casadi::SX>({x, u})));

        int_t nx = x.size1();

        // MATLAB Code:
        // adj = jtimes(f_expl,[x;u],lambdaX,true);
        // expl_vde_adj = Function([model_name,'_expl_vde_adj'],{x,lambdaX,u},{adj});
        casadi::SX lambdaX = casadi::SX::sym("lambdaX", nx, 1);
        casadi::SX w = casadi::SX::vertcat(vector<casadi::SX>({x, u}));
        casadi::SX adj = casadi::SX::jtimes(rhs, w, lambdaX, true);

        fun = casadi::Function(model.name() + "_expl_vde_adj", {x, lambdaX, u}, {adj});
    }
    else  // MX
    {
        casadi::MX x = model.mx_in("x");
        casadi::MX u = model.mx_in("u");

        casadi::MX rhs = casadi::MX::vertcat(model(vector<casadi::MX>({x, u})));

        int_t nx = x.size1();

        casadi::MX lambdaX = casadi::MX::sym("lambdaX", nx, 1);
        casadi::MX w = casadi::MX::vertcat(vector<casadi::MX>({x, u}));
        casadi::MX adj = casadi::MX::jtimes(rhs, w, lambdaX, true);

        fun = casadi::Function(model.name() + "_expl_vde_adj", {x, lambdaX, u}, {adj});
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
        casadi::SX x = model.sx_in("x");
        casadi::SX u = model.sx_in("u");

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
        casadi::MX x = model.mx_in("x");
        casadi::MX u = model.mx_in("u");

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
