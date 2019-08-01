
#include "acados_cpp/integrator.hpp"

#include <algorithm>
#include <exception>

#include "acados_cpp/function_generation.hpp"


namespace acados
{
using std::map;
using std::string;
using std::vector;

static bool check_model(const casadi::Function &model, model_t model_type, const bool use_MX,
                        size_t &nx, size_t &nu, size_t &nz)
{
    /* CHECK inputs */
    int ix = 0;
    int iu = 0;
    int ixdot = 0;
    int iz = 0;
    int model_n_in = model.n_in();
    if (model_n_in < 1 )
        throw std::invalid_argument(
            "An ODE model should have at least 1 input: differential states x");
    // x
    try {
        ix = model.index_in("x");
    } catch (std::exception& e) {
        throw std::invalid_argument(
            "The dynamic system model should have the differential states x as an input"
             "\n please check that this input is named 'x'");
    }

    // u
    try {
        iu = model.index_in("u");
    } catch (std::exception& e) {
        throw std::invalid_argument("The dynamic system model should have the control vector u"
                    "as an input \n please check that this input is named 'u'");
    }

    bool isDAE = false;
    if (model_type == IMPLICIT)
    {
        // xdot
        try {
            ixdot = model.index_in("xdot");
        } catch (std::exception& e) {
            throw std::invalid_argument("An IMPLICIT dynamic system model should have the"
                        "time derivatives of the differential states XDOT as an input"
                        "\n please check that this input is named 'xdot'");
        }

        // z
        try {
            iz = model.index_in("z");
            isDAE = true;
        } catch (std::exception& e) {
            isDAE = false;
        }
    }

    // check on ix, iu, iz, ixdot
    if (ix > 4 || iu > 4 || iz > 4 || ixdot > 4){
        throw std::invalid_argument("The dynamic model seems to have too many inputs");
    }

    /* GET DIMENSIONS */
    if (use_MX == false)
    {
        casadi::SX x = model.sx_in("x");
        casadi::SX u = model.sx_in("u");

        nx = x.size1();
        nu = u.size1();

        if (model_type == IMPLICIT){
            casadi::SX xdot = model.sx_in("xdot");
        }
        if (isDAE){
            casadi::SX z = model.sx_in("z");
            nz = z.size1();
        }
        else nz = 0;
    }
    else  // MX
    {
        casadi::MX x = model.mx_in("x");
        casadi::MX u = model.mx_in("u");

        nx = x.size1();
        nu = u.size1();

        if (model_type == IMPLICIT){
            casadi::MX xdot = model.mx_in("xdot");
        }
        if (isDAE){
            casadi::SX z = model.sx_in("z");
            nz = z.size1();
        }
        else nz = 0;
    }


    /* CHECK outputs */
    if (model.n_out() != 1)
        throw std::invalid_argument("An ODE model should have 1 output: the right hand side\n"
                "-- an explicit model has: xdot\n -- an implicit model has: rhs \n");

    /* CHECK input/output consistency */
    unsigned int size_model_output = model.numel_out();

    if (size_model_output != nx+nz)
        throw std::invalid_argument("The ODE model should have one output,"
               " a column vector of size nx + nz\n");

    return true;
}



/* CONSTRUCTOR */
integrator::integrator(const casadi::Function &model, std::map<std::string, option_t *> options)
{
    if (options.count("model_type"))
        model_type_ = (model_t) to_int(options.at("model_type"));  // This is not a nice interface
    else
        model_type_ = EXPLICIT;

    // check variable type to use
    use_MX_ = false;
    if (options.count("use_MX")) use_MX_ = (to_int(options.at("use_MX")) > 0);

    // check integrator type - default is ERK
    if (options.count("integrator"))
    {
        if (to_string(options.at("integrator")) == "ERK")
            sim_plan_.sim_solver = ERK;
        else if (to_string(options.at("integrator")) == "IRK")
        {
            sim_plan_.sim_solver = IRK;
        }
        else if (to_string(options.at("integrator")) == "LIFTED_IRK")
        {
            sim_plan_.sim_solver = LIFTED_IRK;
        }
        else
            throw std::invalid_argument("Invalid integrator.");
    }
    else  // default integrator
        sim_plan_.sim_solver = ERK;


    if (!check_model(model, model_type_, use_MX_, nx_, nu_, nz_))
        throw std::invalid_argument("Model is invalid.");

    if (!options.count("step_size"))
        throw std::invalid_argument("Expected 'step_size' as an option.");



    config_ = sim_config_create(sim_plan_);

    dims_ = sim_dims_create(config_);

    int n_dim = (int) nx_;
    sim_dims_set(config_, dims_, "nx", &n_dim);
    n_dim = (int) nu_;
    sim_dims_set(config_, dims_, "nu", &n_dim);
    n_dim = (int) nz_;
    sim_dims_set(config_, dims_, "nz", &n_dim);


    // sim opts
    opts_ = static_cast<sim_opts *>(sim_opts_create(config_, dims_));

    if (options.count("sens_forw")) opts_->sens_forw = (to_int(options.at("sens_forw")) >= 1);

    if (options.count("sens_adj")) opts_->sens_adj = (to_int(options.at("sens_adj")) >= 1);

    if (options.count("sens_hess")) opts_->sens_hess = (to_int(options.at("sens_hess")) >= 1);

    if (options.count("jac_reuse")) opts_->jac_reuse = (to_int(options.at("jac_reuse")) >= 1);

    if (options.count("sens_algebraic"))
        opts_->sens_algebraic = (to_int(options.at("sens_algebraic")) >= 1);

    if (options.count("output_z")) opts_->output_z = (to_int(options.at("output_z")) >= 1);

    if (options.count("newton_iter")) opts_->newton_iter = to_int(options.at("newton_iter"));

    if (options.count("num_steps")) opts_->num_steps = to_int(options.at("num_steps"));

    if (options.count("stages")) opts_->ns = to_int(options.at("stages"));

    // sim opts using a void pointer..
    // opts_ = sim_opts_create(config_, dims_);
    // if (options.count("sens_forw"))
    // {
    //     bool sens_forw = (to_int(options.at("sens_forw")) >= 1);
    //     sim_opts_set(config_, opts_, "sens_forw", &sens_forw);
    // }

    // if (options.count("sens_adj"))
    // {
    //     bool sens_adj = (to_int(options.at("sens_adj")) >= 1);
    //     sim_opts_set(config_, opts_, "sens_adj", &sens_adj);
    // }

    // if (options.count("sens_hess"))
    // {
    //     bool sens_hess = (to_int(options.at("sens_hess")) >= 1);
    //     sim_opts_set(config_, opts_, "sens_hess", &sens_hess);
    // }

    // if (options.count("jac_reuse"))
    // {
    //     bool jac_reuse = (to_int(options.at("jac_reuse")) >= 1);
    //     sim_opts_set(config_, opts_, "jac_reuse", &jac_reuse);
    // }

    // if (options.count("sens_algebraic"))
    // {
    //     bool sens_algebraic = (to_int(options.at("sens_algebraic")) >= 1);
    //     sim_opts_set(config_, opts_, "sens_algebraic", &sens_algebraic);
    // }

    // if (options.count("output_z"))
    // {
    //     bool output_z = (to_int(options.at("output_z")) >= 1);
    //     sim_opts_set(config_, opts_, "output_z", &output_z);
    // }


    // if (options.count("newton_iter"))
    // {
    //     int newton_iter = to_int(options.at("newton_iter"));
    //     sim_opts_set(config_, opts_, "newton_iter", &newton_iter);
    // }

    // if (options.count("num_steps"))
    // {
    //     int num_steps = to_int(options.at("num_steps"));
    //     sim_opts_set(config_, opts_, "num_steps", &num_steps);
    // }

    // if (options.count("stages"))
    // {
    //     int stages = to_int(options.at("stages"));
    //     sim_opts_set(config_, opts_, "num_stages", &stages);
    // }

    if (options.count("model_type"))
        model_type_ = (model_t) to_int(options.at("model_type"));
    else
        model_type_ = EXPLICIT;     // default


    if (options.count("use_MX")) {
        use_MX_ = (to_int(options.at("use_MX")) > 0);
    }

    // casadi::Function new_model = model;

    casadi::Function integrator_model = model;
    // reformulate as implicit model, if IMPLICIT integrator is used
    if (model_type_ == EXPLICIT &&
             (sim_plan_.sim_solver == IRK || sim_plan_.sim_solver == LIFTED_IRK))
    {
        integrator_model = explicit2implicit(model);
    }

    // sim_in & sim_out structs in C
    in_ = sim_in_create(config_, dims_);
    out_ = sim_out_create(config_, dims_);


    // set step width
    set_step_size(to_double(options.at("step_size")));

    // generate and set model;
    model_set(integrator_model, options);

    // create the integrator
    solver_ = sim_solver_create(config_, dims_, opts_);
}


casadi::Function integrator::explicit2implicit(const casadi::Function &model)
{
    casadi::Function integrator_model;

    if (use_MX_ == false)  // SX
    {
        casadi::SX x = model.sx_in("x");
        casadi::SX u = model.sx_in("u");
        int_t nx = x.size1();

        casadi::SX rhs = casadi::SX::vertcat(model(vector<casadi::SX>({x, u})));

        casadi::SX xdot = casadi::SX::sym("xdot", nx, 1);
        casadi::SX z = casadi::SX::sym("z", 0, 1);
        rhs = rhs - xdot;

        integrator_model = casadi::Function(model.name(),
                        {x, xdot, u, z},
                        {rhs}, {"x", "xdot", "u", "z"}, {"rhs"});
    }
    else
    {
        casadi::MX x = model.mx_in("x");
        casadi::MX u = model.mx_in("u");
        int_t nx = x.size1();

        casadi::MX rhs = casadi::MX::vertcat(model(vector<casadi::MX>({x, u})));

        casadi::MX xdot = casadi::MX::sym("xdot", nx, 1);
        casadi::MX z = casadi::MX::sym("z", 0, 1);

        rhs = rhs - xdot;

        integrator_model = casadi::Function(model.name(),
                        {x, xdot, u, z},
                        {rhs}, {"x", "xdot", "u", "z"}, {"rhs"});
    }

    model_type_ = IMPLICIT;

    return integrator_model;
}


void integrator::model_set(casadi::Function &model, std::map<std::string, option_t *> options)
{

    string autogen_dir = "_autogen";

    /* generate model functions depending on integrator type and options */
    int model_set_status;
    std::map<std::string, casadi_module> module_;

    if (sim_plan_.sim_solver == IRK)
    {
        if (model_type_ == IMPLICIT)
        {
            module_["impl_ode_fun_jac_x_xdot_z"] =
                    generate_impl_ode_fun_jac_x_xdot_z(model, autogen_dir, use_MX_);
            model_set_status = sim_in_set(config_, dims_, in_,
                "impl_ode_fun_jac_x_xdot_z",
                (void *) module_["impl_ode_fun_jac_x_xdot_z"].as_external_function());

            if (model_set_status == ACADOS_FAILURE)
                    throw std::runtime_error("couldnt set integrator function"
                     " impl_ode_fun_jac_x_xdot_z correctly");

            if (opts_->jac_reuse)
            {
                module_["impl_ode_fun"] = generate_impl_ode_fun(model, autogen_dir, use_MX_);
                model_set_status = sim_in_set(config_, dims_, in_,
                    "impl_ode_fun", (void *) module_["impl_ode_fun"].as_external_function());

                if (model_set_status == ACADOS_FAILURE)
                    throw std::runtime_error("couldnt set integrator function"
                            " impl_ode_fun correctly");
            }
            if ( opts_->sens_forw || opts_->sens_hess || opts_->sens_algebraic || opts_->sens_adj )
            {
                module_["impl_ode_jac_x_xdot_u_z"] =
                            generate_impl_ode_jac_x_xdot_u_z(model, autogen_dir, use_MX_);
                model_set_status = sim_in_set(config_, dims_, in_,
                    "impl_ode_jac_x_xdot_u_z",
                    (void *) module_["impl_ode_jac_x_xdot_u_z"].as_external_function());

                if (model_set_status == ACADOS_FAILURE)
                    throw std::runtime_error("couldnt set integrator function"
                        " impl_ode_jac_x_xdot_u_z correctly");
            }
            if ( opts_->sens_hess )
            {
                module_["impl_ode_hess"] = generate_impl_ode_hess(model, autogen_dir, use_MX_);
                model_set_status = sim_in_set(config_, dims_, in_,
                    "impl_ode_hess", (void *) module_["impl_ode_hess"].as_external_function());

                if (model_set_status == ACADOS_FAILURE)
                    throw std::runtime_error("couldnt set integrator function "
                        "impl_ode_hess correctly");
            }
        }
        else
        {
            throw std::invalid_argument("IRK only supported with implicit model");
        }
    }
    else if (sim_plan_.sim_solver == LIFTED_IRK)
    {
        if (model_type_ == IMPLICIT)
        {
            module_["impl_ode_fun"] = generate_impl_ode_fun(model, autogen_dir, use_MX_);
            model_set_status = sim_in_set(config_, dims_, in_,
                "impl_ode_fun", (void *) module_["impl_ode_fun"].as_external_function());

            if (model_set_status == ACADOS_FAILURE)
                throw std::runtime_error("couldnt set integrator function impl_ode_fun correctly");

            module_["impl_ode_fun_jac_x_xdot_u"] =
                    generate_impl_ode_fun_jac_x_xdot_u(model, autogen_dir, use_MX_);
            model_set_status = sim_in_set(config_, dims_, in_,
                "impl_ode_fun_jac_x_xdot_u",
                (void *) module_["impl_ode_fun_jac_x_xdot_u"].as_external_function());

            if (model_set_status == ACADOS_FAILURE)
                throw std::runtime_error("couldnt set integrator function "
                    "impl_ode_fun_jac_x_xdot_u correctly");
        }
        else
        {
            throw std::invalid_argument("LIFTED_IRK only supported with implicit model");
        }
    }
    else if (sim_plan_.sim_solver == GNSF)
    {
        // TODO(oj): implement
        throw std::invalid_argument("GNSF is not supported integrator type.");
    }
    else  // ERK
    {
        if (model_type_ == EXPLICIT)
        {
            if (opts_->sens_forw)
            {
                module_["expl_vde_for"] = generate_forward_vde(model, autogen_dir, use_MX_);
                model_set_status = sim_in_set(config_, dims_, in_, "expl_vde_for",
                                    (void *) module_["expl_vde_for"].as_external_function());

                if (model_set_status == ACADOS_FAILURE)
                    throw std::runtime_error("couldnt set integrator function "
                        "expl_vde_for correctly");
            }
            else
            {
                module_["expl_ode_fun"] = generate_expl_ode_fun(model, autogen_dir, use_MX_);
                model_set_status = sim_in_set(config_, dims_, in_, "expl_ode_fun",
                                    (void *) module_["expl_ode_fun"].as_external_function());

                if (model_set_status == ACADOS_FAILURE)
                    throw std::runtime_error("couldnt set integrator function "
                        "expl_ode_fun correctly");
            }

            if (opts_->sens_adj && !opts_->sens_hess)
            {
                module_["expl_vde_adj"] = generate_expl_vde_adj(model);
                model_set_status = sim_in_set(config_, dims_, in_, "expl_vde_adj",
                                    (void *) module_["expl_vde_adj"].as_external_function());

                if (model_set_status == ACADOS_FAILURE)
                    throw std::runtime_error("couldnt set integrator function "
                        "expl_vde_adj correctly");
            }
            else if (opts_->sens_hess)
            {
                // throw std::invalid_argument("ERK can only be used without hessians");
                module_["expl_ode_hess"] = generate_expl_ode_hess(model);
                model_set_status = sim_in_set(config_, dims_, in_, "expl_ode_hess",
                                        (void *) module_["expl_ode_hess"].as_external_function());

                if (model_set_status == ACADOS_FAILURE)
                    throw std::runtime_error("couldnt set integrator function expl_ode_hess"
                        " correctly");
            }
        }
        else
        {
            throw std::invalid_argument("ERK can only be used with explicit model");
        }
    }
}


void integrator::print_settings() const
{
    std::cout << "\nstep_size \t: " << in_->T
              << "\nmodel_type \t: " << model_type_
              << "\nintegrator \t: " << sim_plan_.sim_solver
              << "\nuse_MX \t\t: " << use_MX_
              << "\nns \t\t: " << opts_->ns
              << "\nnum_steps \t: " << opts_->num_steps
              << "\nnewton_iter \t: " << opts_->newton_iter
              << "\nsens_forw \t: " << opts_->sens_forw
              << "\nsens_adj \t: " << opts_->sens_adj
              << "\nsens_hess \t: " << opts_->sens_hess
              << "\nsens_algebraic \t: " << opts_->sens_algebraic
              << "\njac_reuse \t: " << opts_->jac_reuse
              << std::endl;
}



void integrator::set_step_size(const double step_size) { in_->T = step_size; }


std::vector<double> integrator::integrate(std::vector<double> x, std::vector<double> u) const
{
    // TODO(mutard): use in/output dict
    /*
      if ( in.count("x") )
          double *x = to_double(in.count("x"));
      else
          throw std::invalid_argument("Missing input x.");

      if ( in.count("u") )
          double *u = to_double(in.count("u"));
      else
          throw std::invalid_argument("Missing input u.");
      */
    if (x.size() != nx_) throw std::invalid_argument("Input x has wrong size.");
    if (u.size() != nu_) throw std::invalid_argument("Input u has wrong size.");

    // mandatory parameters
    sim_in_set(config_, dims_, in_, "x", x.data());
    sim_in_set(config_, dims_, in_, "u", u.data());

    /*
    if ( options.count("sens_forw") )
        _opts->sens_forw = (to_int(options.at("sens_forw")) >= 1);

        // optional parameters
    sim_in_set(config_, dims_, in_, "xdot",  _in);
    sim_in_set(config_, dims_, in_, "Sx"  ,  _in);
    sim_in_set(config_, dims_, in_, "Su"  ,  _in);
    */

    // cast in/out?!
    // int acados_return = sim_solve(solver_, in_, out_);
    sim_solve(solver_, in_, out_);

    std::vector<double> xn(out_->xn, out_->xn + nx_);

    /*
    // new getters
    sim_out_get_xn( _config, _dims, _out, double *xn)
    sim_out_get_Sxn(_config, _dims, _out, double *Sxn)
    sim_out_get_Sun(_config, _dims, _out, double *Sun)
    */

    return xn;
}


/* DESTRUCTOR */
integrator::~integrator()
{
    sim_config_destroy(config_);
    sim_dims_destroy(dims_);
    sim_opts_destroy(opts_);
    sim_in_destroy(in_);
    sim_out_destroy(out_);
    sim_solver_destroy(solver_);
}



}  // namespace acados
