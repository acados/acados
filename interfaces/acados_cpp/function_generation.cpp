
#include "acados_cpp/function_generation.hpp"

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

/************************************************
* IMPLICIT MODEL
************************************************/
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

        casadi::SXDict arg_in =  {{"x", x}, {"xdot", xdot}, {"u", u}, {"z", z}};

        casadi::SXDict rhs_dict = (model(arg_in));
        casadi::SX rhs = rhs_dict.begin()->second;

        casadi::SX jac_x = casadi::SX::jacobian(rhs, x);
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

        casadi::MXDict arg_in =  {{"x", x}, {"xdot", xdot}, {"u", u}, {"z", z}};

        casadi::MXDict rhs_dict = (model(arg_in));
        casadi::MX rhs = rhs_dict.begin()->second;

        casadi::MX jac_x = casadi::MX::jacobian(rhs, x);
        casadi::MX jac_z = casadi::MX::jacobian(rhs, z);
        casadi::MX jac_xdot = casadi::MX::jacobian(rhs, xdot);

        fun = casadi::Function(model.name() + "_impl_ode_fun_jac_x_xdot_z",
                            {x, xdot, u, z}, {rhs, jac_x, jac_xdot, jac_z});
    }
    return casadi_module(fun, output_folder);
}


casadi_module generate_impl_ode_jac_x_xdot_u_z(const casadi::Function& model,
                        string output_folder, const bool use_MX)
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

        casadi::SX jac_x = casadi::SX::jacobian(rhs, x);
        casadi::SX jac_u = casadi::SX::jacobian(rhs, u);
        casadi::SX jac_z = casadi::SX::jacobian(rhs, z);
        casadi::SX jac_xdot = casadi::SX::jacobian(rhs, xdot);

        fun = casadi::Function(model.name() + "_impl_ode_jac_x_xdot_u_z",
                            {x, xdot, u, z}, {jac_x, jac_xdot, jac_u, jac_z});
    }
    else  // MX
    {
        casadi::MX x = model.mx_in("x");
        casadi::MX xdot = model.mx_in("xdot");
        casadi::MX u = model.mx_in("u");
        casadi::MX z = model.mx_in("z");

        casadi::MXDict arg_in =  {{"x", x}, {"xdot", xdot}, {"u", u}, {"z", z}};

        casadi::MXDict rhs_dict = (model(arg_in));
        casadi::MX rhs = rhs_dict.begin()->second;

        casadi::MX jac_x = casadi::MX::jacobian(rhs, x);
        casadi::MX jac_u = casadi::MX::jacobian(rhs, u);
        casadi::MX jac_z = casadi::MX::jacobian(rhs, z);
        casadi::MX jac_xdot = casadi::MX::jacobian(rhs, xdot);

        fun = casadi::Function(model.name() + "_impl_ode_jac_x_xdot_u_z",
                            {x, xdot, u, z}, {jac_x, jac_xdot, jac_u, jac_z});
    }
    return casadi_module(fun, output_folder);
}


casadi_module generate_impl_ode_hess(const casadi::Function& model,
                        string output_folder, const bool use_MX)
{
    casadi::Function fun;
    if (use_MX == false)
    {
        casadi::SX xdot = model.sx_in("xdot");
        casadi::SX x = model.sx_in("x");
        casadi::SX u = model.sx_in("u");
        casadi::SX z = model.sx_in("z");

        int_t nx = x.size1();
        int_t nu = u.size1();
        int_t nz = z.size1();

        casadi::SX x_xdot_z_u = casadi::SX::vertcat(vector<casadi::SX>({x, xdot, z, u}));

        casadi::SX multiplier  = casadi::SX::sym("multiplier", nx + nz, 1);
        casadi::SX mult_mat  = casadi::SX::sym("mult_mat", 2*nx+nz+nu, nx + nu);
        casadi::SX HESS = casadi::SX::zeros(2*nx+nz+nu, 2*nx+nz+nu);

        casadi::SX jac_x_xdot_z;
        casadi::SX hess_x_xdot_z;

        casadi::SXDict arg_in =  {{"x", x}, {"xdot", xdot}, {"u", u}, {"z", z}};

        casadi::SXDict model_dict = (model(arg_in));
        casadi::SX model_expr = model_dict.begin()->second;

        for (int ii = 0; ii < nx+nz; ii++)
        {
            jac_x_xdot_z = jacobian(model_expr(ii), x_xdot_z_u);
            hess_x_xdot_z = jacobian(jac_x_xdot_z, x_xdot_z_u);
            HESS = HESS + multiplier(ii) * hess_x_xdot_z;
        }

        casadi::SX HESS_multiplied = mtimes(mult_mat.T() , mtimes(HESS, mult_mat));

        fun = casadi::Function(model.name() + "_impl_ode_hess",
            {x, xdot, u, z, multiplier, mult_mat}, {HESS_multiplied});
    }
    else
    {
        casadi::MX xdot = model.mx_in("xdot");
        casadi::MX x = model.mx_in("x");
        casadi::MX u = model.mx_in("u");
        casadi::MX z = model.mx_in("z");

        int_t nx = x.size1();
        int_t nu = u.size1();
        int_t nz = z.size1();

        casadi::MX x_xdot_z_u = casadi::MX::vertcat(vector<casadi::MX>({x, xdot, z, u}));

        casadi::MX multiplier  = casadi::MX::sym("multiplier", nx + nz, 1);
        casadi::MX mult_mat  = casadi::MX::sym("mult_mat", 2*nx+nz+nu, nx + nu);
        casadi::MX HESS = casadi::MX::zeros(2*nx+nz+nu, 2*nx+nz+nu);

        casadi::MX jac_x_xdot_z;
        casadi::MX hess_x_xdot_z;

        casadi::MXDict arg_in =  {{"x", x}, {"xdot", xdot}, {"u", u}, {"z", z}};

        casadi::MXDict model_dict = (model(arg_in));
        casadi::MX model_expr = model_dict.begin()->second;

        for (int ii = 0; ii < nx+nz; ii++)
        {
            jac_x_xdot_z = jacobian(model_expr(ii), x_xdot_z_u);
            hess_x_xdot_z = jacobian(jac_x_xdot_z, x_xdot_z_u);
            HESS = HESS + multiplier(ii) * hess_x_xdot_z;
        }

        // HESS = HESS.simplify();
        casadi::MX HESS_multiplied = mtimes(mult_mat.T() , mtimes(HESS, mult_mat));

        fun = casadi::Function(model.name() + "_impl_ode_hess",
            {x, xdot, u, z, multiplier, mult_mat}, {HESS_multiplied});
    }
    return casadi_module(fun, output_folder);
}


/* for LIFTED_IRK */
casadi_module generate_impl_ode_fun_jac_x_xdot_u(const casadi::Function& model,
                        string output_folder, const bool use_MX)
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

        casadi::SX jac_x = casadi::SX::jacobian(rhs, x);
        casadi::SX jac_u = casadi::SX::jacobian(rhs, u);
        casadi::SX jac_xdot = casadi::SX::jacobian(rhs, xdot);

        fun = casadi::Function(model.name() + "_impl_ode_fun_jac_x_xdot_u",
                            {x, xdot, u, z}, {rhs, jac_x, jac_xdot, jac_u});
    }
    else  // MX
    {
        casadi::MX x = model.mx_in("x");
        casadi::MX xdot = model.mx_in("xdot");
        casadi::MX u = model.mx_in("u");
        casadi::MX z = model.mx_in("z");
        casadi::MXDict arg_in =  {{"x", x}, {"xdot", xdot}, {"u", u}, {"z", z}};

        casadi::MXDict rhs_dict = (model(arg_in));
        casadi::MX rhs = rhs_dict.begin()->second;

        casadi::MX jac_x = casadi::MX::jacobian(rhs, x);
        casadi::MX jac_u = casadi::MX::jacobian(rhs, u);
        casadi::MX jac_xdot = casadi::MX::jacobian(rhs, xdot);

        fun = casadi::Function(model.name() + "_impl_ode_fun_jac_x_xdot_u",
                            {x, xdot, u, z}, {rhs, jac_x, jac_xdot, jac_u});
    }
    return casadi_module(fun, output_folder);
}



/************************************************
* EXPLICIT MODEL
************************************************/
casadi_module generate_forward_vde(const casadi::Function& model,
                             string output_folder, const bool use_MX)
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

        casadi::SX rhs = casadi::SX::vertcat(model(vector<casadi::SX>({x, u})));

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

casadi_module generate_expl_ode_hess(const casadi::Function& model, string output_folder,
                                    const bool use_MX)
{
    casadi::Function fun;
    if (use_MX == false)  // SX
    {
        casadi::SX x = model.sx_in("x");
        casadi::SX u = model.sx_in("u");

        casadi::SX rhs = casadi::SX::vertcat(model(vector<casadi::SX>({x, u})));

        int_t nx = x.size1();
        int_t nu = u.size1();

        casadi::SX Sx = casadi::SX::sym("Sx", nx, nx);
        casadi::SX Sp = casadi::SX::sym("Sp", nx, nu);
        casadi::SX lambdaX = casadi::SX::sym("lambdaX", nx, 1);


        casadi::SX w = casadi::SX::vertcat(vector<casadi::SX>({x, u}));
        casadi::SX adj = casadi::SX::jtimes(rhs, w, lambdaX, true);

        std::vector<casadi::SX> SxSp = {Sx, Sp};
        std::vector<casadi::SX> aux = {casadi::SX::zeros(nu, nx), casadi::SX::eye(nu)};

        std::vector<casadi::SX> S_forw_vec {casadi::SX::horzcat(SxSp), casadi::SX::horzcat(aux)};

        casadi::SX S_forw = casadi::SX::vertcat(S_forw_vec);
        casadi::SX hess = S_forw.T() * casadi::SX::jtimes(adj, w, S_forw);

        casadi::SX hess2 = casadi::SX::sym("hess2", 0, 0);

        for (int j = 0; j < nx+nu; j++)
        {
            for (int i = j; i < nx+nu; i++)
            {
                std::vector<casadi::SX> to_concat{hess2, hess(i, j)};
                hess2 = casadi::SX::vertcat(to_concat);
            }
        }
        fun = casadi::Function(model.name() + "_expl_ode_hess",
                                {x, Sx, Sp, lambdaX, u}, {adj, hess2});
    }
    else  // MX
    {
        casadi::MX x = model.mx_in("x");
        casadi::MX u = model.mx_in("u");

        casadi::MX rhs = casadi::MX::vertcat(model(vector<casadi::MX>({x, u})));

        int_t nx = x.size1();
        int_t nu = u.size1();

        casadi::MX Sx = casadi::MX::sym("Sx", nx, nx);
        casadi::MX Sp = casadi::MX::sym("Sp", nx, nu);
        casadi::MX lambdaX = casadi::MX::sym("lambdaX", nx, 1);


        casadi::MX w = casadi::MX::vertcat(vector<casadi::MX>({x, u}));
        casadi::MX adj = casadi::MX::jtimes(rhs, w, lambdaX, true);

        std::vector<casadi::MX> SxSp = {Sx, Sp};
        std::vector<casadi::MX> aux = {casadi::MX::zeros(nu, nx), casadi::MX::eye(nu)};

        std::vector<casadi::MX> S_forw_vec {casadi::MX::horzcat(SxSp), casadi::MX::horzcat(aux)};

        casadi::MX S_forw = casadi::MX::vertcat(S_forw_vec);
        casadi::MX hess = S_forw.T() * casadi::MX::jtimes(adj, w, S_forw);

        casadi::MX hess2 = casadi::MX::sym("hess2", 0, 0);

        for (int j = 0; j < nx+nu; j++)
        {
            for (int i = j; i < nx+nu; i++)
            {
                std::vector<casadi::MX> to_concat{hess2, hess(i, j)};
                hess2 = casadi::MX::vertcat(to_concat);
            }
        }
        fun = casadi::Function(model.name() + "_expl_ode_hess",
                             {x, Sx, Sp, lambdaX, u}, {adj, hess2});
    }

    return casadi_module(fun, output_folder);
}

}  // namespace acados
