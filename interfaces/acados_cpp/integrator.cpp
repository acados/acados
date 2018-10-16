
#include "acados_cpp/integrator.hpp"

#include <algorithm>

#include "acados/utils/types.h"
#include "acados_c/sim_interface.h"


namespace acados
{
using std::map;
using std::string;
using std::vector;

/* CONSTRUCTOR */
integrator::integrator(const casadi::Function &model_fun, std::map<std::string, option_t *> options = {})
{
    if (!is_valid_model(model)) throw std::invalid_argument("Model is invalid.");

    sim_solver_plan sim_plan;

    if ( options.count("integrator") )
    {
        if (to_string(options.at("integrator")) == "ERK")
            sim_plan.sim_solver = ERK;
        else
            throw std::invalid_argument("Invalid integrator.");
    }
    else
            sim_plan.sim_solver = ERK;

    _config = sim_config_create(sim_plan);

    void *dims = sim_dims_create(_config);

    // get dimensions from model
    casadi::SX x = model.sx_in(0);
    casadi::SX u = model.sx_in(1);
    int nx = x.size1();
    int nu = u.size1();

    // set dimensions
    _config->set_nx(dims, nx);
    _config->set_nu(dims, nu);

    // sim opts
    _opts = sim_opts_create(_config, dims);

    if ( options.count("sens_forw") )
        _opts->sens_forw = (to_int(options.at("sens_forw")) >= 1);

    if ( options.count("sens_adj") )
        _opts->sens_adj = (to_int(options.at("sens_adj")) >= 1);
    
    if ( options.count("sens_hess") )
        _opts->sens_hess = (to_int(options.at("sens_hess")) >= 1);

    if ( options.count("jac_reuse") )
        _opts->jac_reuse = (to_int(options.at("jac_reuse")) >= 1);  

    if ( options.count("sens_algebraic") )
        _opts->sens_algebraic = (to_int(options.at("sens_algebraic")) >= 1);

    if ( options.count("output_z") )
        _opts->output_z = (to_int(options.at("output_z")) >= 1);

    if ( options.count("newton_iter") )
        _opts->newton_iter = to_int(options.at("newton_iter"));

    if ( options.count("num_steps") )
        _opts->num_steps = to_int(options.at("num_steps"));

    if ( options.count("stages") )
        _opts->ns = to_int(options.at("stages"));

    // sim_in & sim_out in C
    _in  = sim_in_create(_config, dims);
    _out = sim_out_create(_config, dims);

    // TODO: generate and set model
    // use external_function_generic stuff


    _solver = sim_create(config, dims, opts);

}


std::map<std::string, option_t *> integrator::integrate( std::map<std::string, option_t *> in )
{
    // TODO: get all these doubles from the "in"
    sim_in_set_xdot( _config, _dims, double *xdot,  _in)
    sim_in_set_x(    _config, _dims, double *x,     _in)
    sim_in_set_u(    _config, _dims, double *u,     _in)
    sim_in_set_Sx(   _config, _dims, double *Sx,    _in)
    sim_in_set_Su(   _config, _dims, double *Su,    _in)

    // cast in/out?!
    
    acados_return = sim_solve(sim_solver, _in, _out);

    // new getters
    sim_out_get_xn( _config, _dims, _out, double *xn)
    sim_out_get_Sxn(_config, _dims, _out, double *Sxn)
    sim_out_get_Sun(_config, _dims, _out, double *Sun)

}


/* DESTRUCTOR */
integrator::~integrator()
{
    // @ tobi, do we need the special free functions that just wrap free basically?!
    // edit: talked with giaf about this. we should try to just use sim_interface basicalle
    //      to decouple from core.
    sim_config_free(_config);
    sim_dims_free(_dims);
    sim_opts_free(_opts);
    sim_in_free(_in);
    sim_out_free(_out);
    sim_free(_solver);
}

// Todo: modify this, such that only expression is taken.
static bool is_valid_model(const casadi::Function &model)
{
    if (model.n_in() != 2)
        throw std::invalid_argument("An ODE model should have 2 inputs: states and controls.");
    if (model.n_out() != 1)
        throw std::runtime_error("An ODE model should have 1 output: the right hand side");

    casadi::SX x = model.sx_in(0);
    casadi::SX u = model.sx_in(1);
    int_t nx = x.size1();
    std::vector<casadi::SX> input{x, u};

    casadi::SX rhs = casadi::SX::vertcat(model(input));
    if (rhs.size1() != nx)
        throw std::runtime_error("Length of right hand size should equal number of states");
    return true;
}




}  // namespace acados
