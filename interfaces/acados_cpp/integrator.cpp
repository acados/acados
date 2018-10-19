
#include "acados_cpp/integrator.hpp"

#include <algorithm>

#include "acados_cpp/ocp_nlp/function_generation.hpp"


namespace acados
{
using std::map;
using std::string;
using std::vector;

// Todo: modify this, such that only expression is taken.
static bool is_valid_model(const casadi::Function &model)
{
    if (model.n_in() != 2)
        throw std::invalid_argument("An ODE model should have 2 inputs: states and controls.");
    if (model.n_out() != 1)
        throw std::runtime_error("An ODE model should have 1 output: the right hand side");

    // casadi::SX x = model.sx_in(0);
    casadi::GenericMatrix x = model.sx_in(0);
    casadi::SX u = model.sx_in(1);
    int_t nx = x.size1();
    std::vector<casadi::SX> input{x, u};

    casadi::SX rhs = casadi::SX::vertcat(model(input));  // rhs = f_expl_expr
    if (rhs.size1() != nx)
        throw std::runtime_error("Length of right hand size should equal number of states");
    return true;
}


// TODO: generalize for MX!
/* CONSTRUCTOR */
integrator::integrator(const casadi::Function &model, std::map<std::string, option_t *> options)
{
    if (!is_valid_model(model)) throw std::invalid_argument("Model is invalid.");

    sim_solver_plan sim_plan;

    if (options.count("integrator"))
    {
        if (to_string(options.at("integrator")) == "ERK")
            sim_plan.sim_solver = ERK;
        else
            throw std::invalid_argument("Invalid integrator.");
    }
    else
        sim_plan.sim_solver = ERK;

    if (!options.count("step")) throw std::invalid_argument("Expected 'step' as an option.");

    // TODO add check: modeltype and integrator type consistent

    config_ = sim_config_create(sim_plan);

    dims_ = sim_dims_create(config_);

    // get dimensions from model
    casadi::SX x = model.sx_in(0);
    casadi::SX u = model.sx_in(1);
    nx_ = x.size1();
    nu_ = u.size1();

    // set dimensions
    config_->set_nx(dims_, nx_);
    config_->set_nu(dims_, nu_);

    // sim opts
    opts_ = static_cast<sim_rk_opts *>(sim_opts_create(config_, dims_));

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

    // sim_in & sim_out in C
    in_ = sim_in_create(config_, dims_);
    out_ = sim_out_create(config_, dims_);

    // set step width
    set_step(to_double(options.at("step")));

    // TODO: generate and set model;
    // use external_function_generic stuff
    // std::map<std::string, option_t *> model_options = {};
    // if (options.count("model_type")) model_options[model_type] =
    // to_int(options.at("model_type"));
    use_MX_ = false;
    integrator::set_model(model, options);

    solver_ = sim_create(config_, dims_, opts_);
}


void integrator::set_model(const casadi::Function &model, std::map<std::string, option_t *> options)
{
    if (options.count("model_type")) model_type_ = (model_t) to_int(options.at("model_type"));

    if (options.count("use_MX")) use_MX_ = (to_int(options.at("use_MX")>0);

    if (!use_MX_){
        if (model_type_ == GNSF)
        {
            throw std::invalid_argument("Not supported model type.");
        }
        else if (model_type_ == IMPLICIT)
        {
            // module_["impl_ode_fun"] = generate_impl_ode_fun(model);
            // int status = sim_set_model_internal(config_, in_->model,
            //     "impl_ode_fun", (void *) module_["impl_ode_fun"].as_external_function());

            throw std::invalid_argument("Not supported model type.");
        }
        else  // explicit default
        {
            if (opts_->sens_forw)
            {
                module_["expl_vde_for"] = generate_forward_vde(model);
                int status =
                    sim_set_model_internal(config_, in_->model, "expl_vde_for",
                                        (void *) module_["expl_vde_for"].as_external_function());
            }
            else
            {
                module_["expl_ode_fun"] = generate_expl_ode_fun(model);
                int status =
                    sim_set_model_internal(config_, in_->model, "expl_ode_fun",
                                        (void *) module_["expl_ode_fun"].as_external_function());
            }
        }
    }
    else

    // todo: write generators for all types of functions
}


void integrator::set_step(const double step) { in_->T = step; }


std::vector<double> integrator::integrate(std::vector<double> x, std::vector<double> u)
{
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
    sim_in_set_x(config_, dims_, x.data(), in_);
    sim_in_set_u(config_, dims_, u.data(), in_);

    /*
    if ( options.count("sens_forw") )
        _opts->sens_forw = (to_int(options.at("sens_forw")) >= 1);

        // optional parameters
    sim_in_set_xdot( _config, _dims, double *xdot,  _in)
    sim_in_set_Sx(   _config, _dims, double *Sx,    _in)
    sim_in_set_Su(   _config, _dims, double *Su,    _in)
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
    // @ tobi, do we need the special free functions that just wrap free basically?!
    // edit: talked with giaf about this. we should try to just use sim_interface basicalle
    //      to decouple from core.
    sim_config_free(config_);
    sim_dims_free(dims_);
    sim_opts_free(opts_);
    sim_in_free(in_);
    sim_out_free(out_);
    sim_free(solver_);
}



}  // namespace acados
