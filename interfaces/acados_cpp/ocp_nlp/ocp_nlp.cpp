
#include "acados_cpp/ocp_nlp/ocp_nlp.hpp"

#include <string>

#include "acados/ocp_nlp/ocp_nlp_dynamics_cont.h"
#include "acados/ocp_nlp/ocp_nlp_cost_ls.h"
#include "acados/ocp_nlp/ocp_nlp_sqp.h"
#include "acados_cpp/ocp_dimensions.hpp"

#include "casadi/casadi.hpp"
#include "casadi/mem.h"

#if (defined _WIN32 || defined _WIN64 || defined __MINGW32__ || defined __MINGW64__)
#include <windows.h>
std::string compiler {"mex"};
#else
#include <dlfcn.h>
std::string compiler {"cc"};
#endif
int global_library_counter = 0;

namespace acados {

using std::vector;
using std::string;

ocp_nlp::ocp_nlp(std::vector<int> nx, std::vector<int> nu, std::vector<int> ng, std::vector<int> nh,
                 std::vector<int> ns)
                 : N(nx.size()-1), needs_initializing_(true) {
    
    // Number of controls on last stage should be zero;
    if (!nu.empty()) nu.back() = 0;

    std::map<std::string, std::vector<int>> dims {
        {"nx", nx}, {"nu", nu}, {"nbx", nx}, {"nbu", nu}, {"ng", ng}, {"nh", nh}, {"ns", ns}
    };

    if (!are_valid_ocp_dimensions(dims, {"nx", "nu", "nbx", "nbu", "nb", "ng", "nh", "ns"}))
        throw std::invalid_argument("Invalid dimensions map.");

    dimensions_["nx"] = nx;
    dimensions_["nu"] = nu;
    dimensions_["nbx"] = nx;
    dimensions_["nbu"] = nu;
    dimensions_["ng"] = ng;
    dimensions_["nh"] = nh;
    dimensions_["ns"] = ns;

    std::vector<int> ny(N+1, 0);
    std::transform(nx.begin(), nx.end(), nu.begin(), ny.begin(), std::plus<int>());
    dimensions_["ny"] = ny;

    int bytes = ocp_nlp_solver_config_calculate_size(N);
	void *config_mem = calloc(1, bytes);
	config_.reset(ocp_nlp_solver_config_assign(N, config_mem));

    config_->N = N;
    ocp_nlp_sqp_config_initialize_default(config_.get());
    config_->qp_solver = ocp_qp_config_create({PARTIAL_CONDENSING_HPIPM});

    for (int i = 0; i <= N; ++i) {
        ocp_nlp_cost_ls_config_initialize_default(config_->cost[i]);
        ocp_nlp_constraints_config_initialize_default(config_->constraints[i]);
    }

    for (int stage = 0; stage <= N; ++stage) {
        cached_bounds["lbx"].push_back(vector<double>(nx[stage], -INFINITY));
        cached_bounds["ubx"].push_back(vector<double>(nx[stage], +INFINITY));
        cached_bounds["lbu"].push_back(vector<double>(nu[stage], -INFINITY));
        cached_bounds["ubu"].push_back(vector<double>(nu[stage], +INFINITY));
    }

}

ocp_nlp::ocp_nlp(int N, int nx, int nu, int ng, int nh, int ns)
                 : ocp_nlp(std::vector<int>(N+1, nx), std::vector<int>(N+1, nu),
                           std::vector<int>(N+1, ng), std::vector<int>(N+1, nh),
                           std::vector<int>(N+1, ns)) {}

void ocp_nlp::initialize_solver(std::string solver_name, std::map<std::string, option_t *> options) {

    dims_.reset(ocp_nlp_dims_create(config_.get()));
    ocp_nlp_dims_initialize(config_.get(), dimensions_["nx"].data(), dimensions_["nu"].data(),
                            dimensions_["ny"].data(), dimensions_["nbx"].data(), dimensions_["nbu"].data(),
                            dimensions_["ng"].data(), dimensions_["nh"].data(), std::vector<int>(N+1, 0).data(),
                            dimensions_["ns"].data(), dims_.get());

    nlp_.reset(ocp_nlp_in_create(config_.get(), dims_.get()));

    void *fcn_handle = load_function("expl_vde_for", dynamics_handle["expl_vde_for"]);

    for (int stage = 0; stage < N; ++stage)
        nlp_set_model_in_stage(config_.get(), nlp_.get(), stage, "expl_vde_for", fcn_handle);

    squeeze_dimensions(cached_bounds);
}

// ocp_nlp_solution ocp_nlp::solve() {

//     // fill_in_bounds();

//     // config_->evaluate();

// }

void ocp_nlp::set_field(string field, vector<double> v) {
    for (int stage = 0; stage <= N; ++stage)
        set_field(field, stage, v);
}

void ocp_nlp::set_field(string field, int stage, vector<double> v) {

    if (field == "lbx" || field == "ubx") {

        if (v.size() != dimensions_["nx"].at(stage))
            throw std::invalid_argument("Expected size " + std::to_string(dimensions_["nx"].at(stage))
                                        + " but got " + std::to_string(v.size()) + " instead.");

        cached_bounds[field].at(stage) = v;

    } else if (field == "lbu" || field == "ubu") {

        if (v.size() != dimensions_["nu"].at(stage))
            throw std::invalid_argument("Expected size " + std::to_string(dimensions_["nu"].at(stage))
                                        + " but got " + std::to_string(v.size()) + " instead.");
        
        cached_bounds[field].at(stage) = v;

    } else {
        throw std::invalid_argument("OCP NLP does not contain field '" + field + "'.");
    }

}

static bool is_valid_model(const casadi::Function& model) {
    if (model.n_in() != 2 && model.n_in() != 3)
        throw std::invalid_argument("An ODE model should have 2 inputs: states and controls.");
    if (model.n_out() != 1)
        throw std::runtime_error("An ODE model should have 1 output: the right hand side");
    casadi::SX x = model.sx_in(0);
    casadi::SX u = model.sx_in(1);
    int_t nx = x.size1();
    std::vector<casadi::SX> input {x, u};

    casadi::SX rhs = casadi::SX::vertcat(model(input));
    if (rhs.size1() != nx)
        throw std::runtime_error("Length of right hand size should equal number of states");
    return true;
}

void ocp_nlp::set_dynamics(const casadi::Function& model, std::map<std::string, option_t *> options) {

    if (!is_valid_model(model))
        throw std::invalid_argument("Model is invalid.");

    auto functions = create_explicit_ode_functions(model);

    for (auto elem : functions) {
        elem.second.generate(elem.first + ".c", {{"with_header", true}, {"with_export", false}});
        dynamics_handle[elem.first] = compile_and_load(elem.first);
    }

    sim_solver_plan sim_plan;
    if (to_string(options.at("integrator")) == "rk4")
        sim_plan.sim_solver = ERK;
    else
        throw std::invalid_argument("Invalid integrator.");
    
    for (int i = 0; i < N; ++i) {
        ocp_nlp_dynamics_cont_config_initialize_default(config_->dynamics[i]);
        config_->dynamics[i]->sim_solver = sim_config_create(sim_plan);
    }
};

void *load_function(std::string function_name, void *handle) {
    #if (defined _WIN32 || defined _WIN64 || defined __MINGW32__ || defined __MINGW64__)
        return GetProcAddress((HMODULE) handle, function_name.c_str());
    #else
        return dlsym(handle, function_name.c_str());
    #endif
}

static std::string load_error_message() {
    #if (defined _WIN32 || defined _WIN64 || defined __MINGW32__ || defined __MINGW64__)

    // Retrieve the system error message for the last-error code
    LPVOID lpMsgBuf;
    DWORD dw = GetLastError();

    FormatMessage(
        FORMAT_MESSAGE_ALLOCATE_BUFFER |
        FORMAT_MESSAGE_FROM_SYSTEM |
        FORMAT_MESSAGE_IGNORE_INSERTS,
        NULL,
        dw,
        MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
        (LPTSTR) &lpMsgBuf,
        0, NULL);

    return std::string((LPTSTR) lpMsgBuf);

    #else

    return std::string(dlerror());

    #endif

}

std::map<std::string, casadi::Function> ocp_nlp::create_explicit_ode_functions(const casadi::Function& model) {
    casadi::SX x = model.sx_in(0);
    casadi::SX u = model.sx_in(1);
    int_t nx = x.size1();
    int_t nu = u.size1();
    casadi::SX rhs = casadi::SX::vertcat(model(std::vector<casadi::SX>({x, u})));

    casadi::SX Sx = casadi::SX::sym("Sx", nx, nx);
    casadi::SX Su = casadi::SX::sym("Su", nx, nu);

    casadi::SX vde_x = casadi::SX::jtimes(rhs, x, Sx);
    casadi::SX vde_u = casadi::SX::jacobian(rhs, u) + casadi::SX::jtimes(rhs, x, Su);

    casadi::Function vde_fun("expl_vde_for", {x, Sx, Su, u}, {rhs, vde_x, vde_u});

    casadi::SX jac_x = casadi::SX::zeros(nx, nx) + casadi::SX::jacobian(rhs, x);
    casadi::Function jac_fun("expl_ode_jac" + model.name(), {x, Sx, Su, u}, {rhs, jac_x});

    return {{"expl_vde_for", vde_fun}, {"expl_ode_jac", jac_fun}};
}

void *compile_and_load(std::string name) {
    void *handle;
    std::string library_name = name + std::to_string(global_library_counter++) + std::string(".so");
    std::string path_to_library = std::string("./") + library_name;
    char command[MAX_STR_LEN];
    snprintf(command, sizeof(command), "%s -fPIC -shared -g %s.c -o %s", compiler.c_str(), name.c_str(),
        library_name.c_str());
    int compilation_failed = system(command);
    if (compilation_failed)
        throw std::runtime_error("Something went wrong when compiling the model.");
    #if (defined _WIN32 || defined _WIN64 || defined __MINGW32__ || defined __MINGW64__)
    handle = LoadLibrary(path_to_library.c_str());
    #else
    handle = dlopen(path_to_library.c_str(), RTLD_LAZY);
    #endif
    if (handle == NULL)
        throw std::runtime_error("Loading of " + path_to_library + " failed. Error message: "
                                 + load_error_message());
    return handle;
}

void ocp_nlp::squeeze_dimensions(std::map<std::string, std::vector<std::vector<double>>> bounds) {
    ocp::squeeze_dimensions(bounds);
}

void ocp_nlp::change_bound_dimensions(std::vector<int> nbx, std::vector<int> nbu) {

    if (nbx.size() != (num_stages()+1) || nbu.size() != (num_stages()+1))
        throw std::invalid_argument("Expected " + std::to_string(num_stages()) + " bounds.");

	ocp_nlp_constraints_dims **constraint_dims = (ocp_nlp_constraints_dims **) dims_->constraints;

    for (int i = 0; i <= N; ++i) {
        constraint_dims[i]->nbx = nbx[i];
        constraint_dims[i]->nbu = nbu[i];
    }

    dimensions_["nbx"] = nbx;
    dimensions_["nbu"] = nbu;
}

bool ocp_nlp::needs_initializing() {
    return needs_initializing_;
}

void ocp_nlp::needs_initializing(bool flag) {
    needs_initializing_ = flag;
}

void ocp_nlp::set_bounds_indices(std::string bound, int stage, std::vector<int> idxb) {

	ocp_nlp_constraints_model **constraints = (ocp_nlp_constraints_model **) nlp_->constraints;
	ocp_nlp_constraints_dims **constraint_dims = (ocp_nlp_constraints_dims **) dims_->constraints;

    if (bound == "x") {
        if (idxb.size() != constraint_dims[stage]->nbx)
            throw std::invalid_argument("Expected " + std::to_string(constraint_dims[stage]->nbx)
                                        + " bound indices but got " + std::to_string(idxb.size()) + ".");
        for (int i = 0; i < constraint_dims[stage]->nbx; ++i)
            constraints[stage]->idxb[constraint_dims[stage]->nbu + i] = idxb[i];
    } else if (bound == "u") {
        if (idxb.size() != constraint_dims[stage]->nbu)
            throw std::invalid_argument("Expected " + std::to_string(constraint_dims[stage]->nbu)
                                        + " bound indices but got " + std::to_string(idxb.size()) + ".");
        for (int i = 0; i < constraint_dims[stage]->nbu; ++i)
            constraints[stage]->idxb[i] = idxb[i];
    } else {
        throw std::invalid_argument("Can only set bounds for 'x' and 'u'.");
    }

}


int ocp_nlp::num_stages() {
    return N;
}

}  // namespace acados
