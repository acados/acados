
#include "acados_cpp/ocp_nlp/ocp_nlp.hpp"

#include <algorithm>
#include <functional>
#include <string>

#include "acados/ocp_nlp/ocp_nlp_cost_ls.h"
#include "acados/ocp_nlp/ocp_nlp_dynamics_cont.h"
#include "acados/ocp_nlp/ocp_nlp_sqp.h"

#include "acados_c/external_function_interface.h"
#include "acados_cpp/ocp_bounds.hpp"
#include "acados_cpp/ocp_dimensions.hpp"
#include "acados_cpp/utils.hpp"

#include "blasfeo/include/blasfeo_d_aux.h"

#include "casadi/casadi.hpp"
#include "casadi/mem.h"

#if (defined _WIN32 || defined _WIN64 || defined __MINGW32__ || defined __MINGW64__)
#include <windows.h>
#else
#include <dlfcn.h>
#endif
int global_library_counter = 0;

namespace acados
{
using std::map;
using std::vector;
using std::string;

ocp_nlp::ocp_nlp(std::vector<int> nx, std::vector<int> nu, std::vector<int> ng, std::vector<int> nh,
                 std::vector<int> ns)
    : N(nx.size() - 1), needs_initializing_(true)
{
    // Number of controls on last stage should be zero;
    if (!nu.empty()) nu.back() = 0;

    std::map<std::string, std::vector<int>> dims{{"nx", nx}, {"nu", nu}, {"nbx", nx}, {"nbu", nu},
                                                 {"ng", ng}, {"nh", nh}, {"ns", ns}};

    if (!are_valid_ocp_dimensions(dims, {"nx", "nu", "nbx", "nbu", "nb", "ng", "nh", "ns"}))
        throw std::invalid_argument("Invalid dimensions map.");

    dimensions_["nx"] = nx;
    dimensions_["nu"] = nu;
    dimensions_["nbx"] = nx;
    dimensions_["nbu"] = nu;
    dimensions_["ng"] = ng;
    dimensions_["nh"] = nh;
    dimensions_["ns"] = ns;

    std::vector<int> ny(N + 1, 0);
    std::transform(nx.begin(), nx.end(), nu.begin(), ny.begin(), std::plus<int>());
    dimensions_["ny"] = ny;

    step_ = std::vector<double>(N);

    int bytes = ocp_nlp_solver_config_calculate_size(N);
    void *config_mem = calloc(1, bytes);
    config_.reset(ocp_nlp_solver_config_assign(N, config_mem));

    config_->N = N;
    ocp_nlp_sqp_config_initialize_default(config_.get());
    config_->qp_solver = ocp_qp_config_create({PARTIAL_CONDENSING_HPIPM});

    for (int i = 0; i <= N; ++i)
    {
        ocp_nlp_cost_ls_config_initialize_default(config_->cost[i]);
        ocp_nlp_constraints_config_initialize_default(config_->constraints[i]);
    }

    nlp_.reset(new ocp_nlp_in);
    nlp_->constraints = (void **) calloc(N + 1, sizeof(void *));
    nlp_->cost = (void **) calloc(N + 1, sizeof(void *));
    nlp_->dynamics = (void **) calloc(N, sizeof(void *));

    dims_.reset(new ocp_nlp_dims);
    dims_->N = N;
    dims_->nv = (int *) calloc(N + 1, sizeof(int));
    dims_->nx = (int *) calloc(N + 1, sizeof(int));
    dims_->nu = (int *) calloc(N + 1, sizeof(int));
    dims_->ni = (int *) calloc(N + 1, sizeof(int));
    dims_->constraints = (void **) calloc(N + 1, sizeof(void *));
    dims_->cost = (void **) calloc(N + 1, sizeof(void *));
    dims_->dynamics = (void **) calloc(N, sizeof(void *));

    for (int i = 0; i <= N; ++i)
    {
        dims_->nv[i] = dimensions_["nx"][i] + dimensions_["nu"][i] + dimensions_["ns"][i];
        dims_->nx[i] = dimensions_["nx"][i];
        dims_->nu[i] = dimensions_["nu"][i];
        dims_->ni[i] = dimensions_["nbx"][i] + dimensions_["nbu"][i] + dimensions_["ng"][i] +
                       dimensions_["nh"][i];
    }

    for (int i = 0; i <= N; ++i)
    {
        dims_->constraints[i] = malloc(sizeof(ocp_nlp_constraints_dims));
        ocp_nlp_constraints_dims *con_dims = (ocp_nlp_constraints_dims *) dims_->constraints[i];
        con_dims->nx = dimensions_["nx"][i];
        con_dims->nu = dimensions_["nu"][i];
        con_dims->nb = dimensions_["nbx"][i] + dimensions_["nbu"][i];
        con_dims->nbx = dimensions_["nbx"][i];
        con_dims->nbu = dimensions_["nbu"][i];
        con_dims->ng = dimensions_["ng"][i];
        con_dims->nh = dimensions_["nh"][i];
        con_dims->np = 0;
        con_dims->ns = dimensions_["ns"][i];
        int num_bytes = ocp_nlp_constraints_model_calculate_size(config_->constraints[i], con_dims);
        void *raw_memory = calloc(1, num_bytes);
        nlp_->constraints[i] =
            ocp_nlp_constraints_model_assign(config_->constraints[i], con_dims, raw_memory);

        dims_->cost[i] = malloc(sizeof(ocp_nlp_cost_ls_dims));
        ocp_nlp_cost_ls_dims *cost_dims = (ocp_nlp_cost_ls_dims *) dims_->cost[i];
        cost_dims->nx = dimensions_["nx"][i];
        cost_dims->nu = dimensions_["nu"][i];
        cost_dims->ns = dimensions_["ns"][i];
        cost_dims->ny = dimensions_["ny"][i];
        num_bytes = ocp_nlp_cost_ls_model_calculate_size(config_->cost[i], cost_dims);
        raw_memory = calloc(1, num_bytes);
        nlp_->cost[i] = ocp_nlp_cost_ls_model_assign(config_->cost[i], cost_dims, raw_memory);
    }

    int num_bytes = ocp_qp_dims_calculate_size(N);
    void *raw_memory = calloc(1, num_bytes);
    dims_->qp_solver = ocp_qp_dims_assign(N, raw_memory);

    nlp_->dims = dims_.get();

    for (int stage = 0; stage <= N; ++stage)
    {
        cached_bounds["lbx"].push_back(vector<double>(nx[stage], -INFINITY));
        cached_bounds["ubx"].push_back(vector<double>(nx[stage], +INFINITY));
        cached_bounds["lbu"].push_back(vector<double>(nu[stage], -INFINITY));
        cached_bounds["ubu"].push_back(vector<double>(nu[stage], +INFINITY));
    }
}

ocp_nlp::ocp_nlp(int N, int nx, int nu, int ng, int nh, int ns)
    : ocp_nlp(std::vector<int>(N + 1, nx), std::vector<int>(N + 1, nu), std::vector<int>(N + 1, ng),
              std::vector<int>(N + 1, nh), std::vector<int>(N + 1, ns))
{
}

void ocp_nlp::initialize_solver(std::string solver_name, std::map<std::string, option_t *> options)
{
    squeeze_dimensions(cached_bounds);

    ocp_nlp_dims_initialize(config_.get(), dimensions_["nx"].data(), dimensions_["nu"].data(),
                            dimensions_["ny"].data(), dimensions_["nbx"].data(),
                            dimensions_["nbu"].data(), dimensions_["ng"].data(),
                            dimensions_["nh"].data(), std::vector<int>(N + 1, 0).data(),
                            dimensions_["ns"].data(), nlp_->dims);

    nlp_->Ts = step_.data();

    vector<double> xref(2, 0);
    vector<double> uref(1, 0);

    for (int i = 0; i <= N; ++i)
    {
        ocp_nlp_cost_ls_model *stage_cost_ls = (ocp_nlp_cost_ls_model *) nlp_->cost[i];
        // Cyt
        blasfeo_dgese(dimensions_["nu"][i] + dimensions_["nx"][i], dimensions_["ny"][i], 0.0,
                      &stage_cost_ls->Cyt, 0, 0);
        for (int j = 0; j < dimensions_["nu"][i]; j++)
            BLASFEO_DMATEL(&stage_cost_ls->Cyt, j, dimensions_["nx"][i] + j) = 1.0;
        for (int j = 0; j < dimensions_["nx"][i]; j++)
            BLASFEO_DMATEL(&stage_cost_ls->Cyt, dimensions_["nu"][i] + j, j) = 1.0;

        // W
        blasfeo_dgese(dimensions_["ny"][i], dimensions_["ny"][i], 0.0, &stage_cost_ls->W, 0, 0);
        for (int j = 0; j < dimensions_["nx"][i]; j++)
            BLASFEO_DMATEL(&stage_cost_ls->W, j, j) = 100.0;
        for (int j = 0; j < dimensions_["nu"][i]; j++)
            BLASFEO_DMATEL(&stage_cost_ls->W, dimensions_["nx"][i] + j, dimensions_["nx"][i] + j) =
                1.0;

        // y_ref
        blasfeo_pack_dvec(dimensions_["nx"][i], xref.data(), &stage_cost_ls->y_ref, 0);
        blasfeo_pack_dvec(dimensions_["nu"][i], uref.data(), &stage_cost_ls->y_ref,
                          dimensions_["nx"][i]);
    }

    solver_options_.reset(ocp_nlp_opts_create(config_.get(), dims_.get()));
}

ocp_nlp_solution ocp_nlp::solve()
{
    fill_bounds(cached_bounds);

    solver_.reset(ocp_nlp_create(config_.get(), dims_.get(), solver_options_.get()));

    auto result = std::shared_ptr<ocp_nlp_out>(ocp_nlp_out_create(config_.get(), dims_.get()));

    ocp_nlp_solve(solver_.get(), nlp_.get(), result.get());

    return ocp_nlp_solution(result, dims_);
}

void ocp_nlp::set_field(string field, vector<double> v)
{
    for (int stage = 0; stage <= N; ++stage) set_field(field, stage, v);
}

void ocp_nlp::set_field(string field, int stage, vector<double> v)
{
    if (field == "lbx" || field == "ubx")
    {
        if (v.size() != (size_t) dimensions_["nx"].at(stage))
            throw std::invalid_argument("Expected size " +
                                        std::to_string(dimensions_["nx"].at(stage)) + " but got " +
                                        std::to_string(v.size()) + " instead.");

        cached_bounds[field].at(stage) = v;
    }
    else if (field == "lbu" || field == "ubu")
    {
        if (v.size() != (size_t) dimensions_["nu"].at(stage))
            throw std::invalid_argument("Expected size " +
                                        std::to_string(dimensions_["nu"].at(stage)) + " but got " +
                                        std::to_string(v.size()) + " instead.");

        cached_bounds[field].at(stage) = v;
    }
    else
    {
        throw std::invalid_argument("OCP NLP does not contain field '" + field + "'.");
    }
}

static bool is_valid_model(const casadi::Function &model)
{
    if (model.n_in() != 2 && model.n_in() != 3)
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

void ocp_nlp::set_dynamics(const casadi::Function &model, std::map<std::string, option_t *> options)
{
    if (!is_valid_model(model)) throw std::invalid_argument("Model is invalid.");

    if (!options.count("step")) throw std::invalid_argument("Expected 'step' as an option.");

    std::fill(step_.begin(), step_.end(), to_double(options["step"]));

    sim_solver_plan sim_plan;
    if (to_string(options.at("integrator")) == "rk4")
        sim_plan.sim_solver = ERK;
    else
        throw std::invalid_argument("Invalid integrator.");

    for (int i = 0; i < N; ++i)
    {
        ocp_nlp_dynamics_cont_config_initialize_default(config_->dynamics[i]);
        config_->dynamics[i]->sim_solver = sim_config_create(sim_plan);
    }

    for (int i = 0; i < N; ++i)
    {
        int num_bytes = ocp_nlp_dynamics_cont_dims_calculate_size(config_->dynamics[i]);
        void *raw_memory = calloc(1, num_bytes);
        dims_->dynamics[i] = ocp_nlp_dynamics_cont_dims_assign(config_->dynamics[i], raw_memory);
        ocp_nlp_dynamics_cont_dims_initialize(config_->dynamics[i], dims_->dynamics[i],
                                              dimensions_["nx"][i], dimensions_["nu"][i],
                                              dimensions_["nx"][i + 1], dimensions_["nu"][i + 1]);

        num_bytes =
            ocp_nlp_dynamics_cont_model_calculate_size(config_->dynamics[i], dims_->dynamics[i]);
        raw_memory = calloc(1, num_bytes);
        nlp_->dynamics[i] = ocp_nlp_dynamics_cont_model_assign(config_->dynamics[i],
                                                               dims_->dynamics[i], raw_memory);
    }

    auto functions = create_explicit_ode_functions(model);

    for (auto elem : functions)
    {
        elem.second.generate(
            elem.first + ".c",
            {{"with_header", true}, {"with_export", false}, {"casadi_int", "int"}});
        dynamics_handle_[elem.first] = compile_and_load(elem.first);
    }

    forw_vde_.casadi_fun = reinterpret_cast<casadi_eval_t>(
        load_function("expl_vde_for", dynamics_handle_["expl_vde_for"]));
    forw_vde_.casadi_n_in = reinterpret_cast<casadi_getint_t>(
        load_function("expl_vde_for_n_in", dynamics_handle_["expl_vde_for"]));
    forw_vde_.casadi_n_out = reinterpret_cast<casadi_getint_t>(
        load_function("expl_vde_for_n_out", dynamics_handle_["expl_vde_for"]));
    forw_vde_.casadi_sparsity_in = reinterpret_cast<casadi_sparsity_t>(
        load_function("expl_vde_for_sparsity_in", dynamics_handle_["expl_vde_for"]));
    forw_vde_.casadi_sparsity_out = reinterpret_cast<casadi_sparsity_t>(
        load_function("expl_vde_for_sparsity_out", dynamics_handle_["expl_vde_for"]));
    forw_vde_.casadi_work = reinterpret_cast<casadi_work_t>(
        load_function("expl_vde_for_work", dynamics_handle_["expl_vde_for"]));

    external_function_casadi_create_array(1, &forw_vde_);

    for (int stage = 0; stage < N; ++stage)
        nlp_set_model_in_stage(config_.get(), nlp_.get(), stage, "expl_vde_for", &forw_vde_);
};

void *load_function(std::string function_name, void *handle)
{
#if (defined _WIN32 || defined _WIN64 || defined __MINGW32__ || defined __MINGW64__)
    return GetProcAddress((HMODULE) handle, function_name.c_str());
#else
    return dlsym(handle, function_name.c_str());
#endif
}

std::map<std::string, casadi::Function> ocp_nlp::create_explicit_ode_functions(
    const casadi::Function &model)
{
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

void *compile_and_load(std::string name)
{
#if (defined _WIN32 || defined _WIN64 || defined __MINGW32__ || defined __MINGW64__)
    std::string compiler{"mex"};
#else
    std::string compiler{"cc"};
#endif
    void *handle;
    std::string library_name = name + std::to_string(global_library_counter++) + std::string(".so");
    std::string path_to_library = std::string("./") + library_name;
    char command[MAX_STR_LEN];
    snprintf(command, sizeof(command), "%s -fPIC -shared -g %s.c -o %s", compiler.c_str(),
             name.c_str(), library_name.c_str());
    int compilation_failed = system(command);
    if (compilation_failed)
        throw std::runtime_error("Something went wrong when compiling the model.");
#if (defined _WIN32 || defined _WIN64 || defined __MINGW32__ || defined __MINGW64__)
    handle = LoadLibrary(path_to_library.c_str());
#else
    handle = dlopen(path_to_library.c_str(), RTLD_LAZY);
#endif
    if (handle == NULL)
        throw std::runtime_error("Loading of " + path_to_library +
                                 " failed. Error message: " + load_error_message());
    return handle;
}

void ocp_nlp::set_bound(std::string bound, int stage, std::vector<double> new_bound)
{
    ocp_nlp_constraints_model **constraints = (ocp_nlp_constraints_model **) nlp_->constraints;
    ocp_nlp_constraints_dims **constraint_dims = (ocp_nlp_constraints_dims **) dims_->constraints;

    int nbx = constraint_dims[stage]->nbx;
    int nbu = constraint_dims[stage]->nbu;
    int ni = dims_->ni[stage];

    if (bound == "lbx")
    {
        blasfeo_pack_dvec(nbx, new_bound.data(), &constraints[stage]->d, nbu);
    }
    else if (bound == "ubx")
    {
        blasfeo_pack_dvec(nbx, new_bound.data(), &constraints[stage]->d, ni + nbu);
    }
    else if (bound == "lbu")
    {
        blasfeo_pack_dvec(nbu, new_bound.data(), &constraints[stage]->d, 0);
    }
    else if (bound == "ubu")
    {
        blasfeo_pack_dvec(nbu, new_bound.data(), &constraints[stage]->d, ni);
    }
    else
    {
        throw std::invalid_argument("Expected one of {'lbx', 'ubx', 'lbu', 'ubu'} but got " +
                                    bound + ".");
    }
}

void ocp_nlp::change_bound_dimensions(std::vector<int> nbx, std::vector<int> nbu)
{
    if (nbx.size() != (size_t)(num_stages() + 1) || nbu.size() != (size_t)(num_stages() + 1))
        throw std::invalid_argument("Expected " + std::to_string(num_stages()) + " bounds.");

    ocp_nlp_constraints_dims **constraint_dims = (ocp_nlp_constraints_dims **) dims_->constraints;

    for (int i = 0; i <= N; ++i)
    {
        constraint_dims[i]->nbx = nbx[i];
        constraint_dims[i]->nbu = nbu[i];
    }

    dimensions_["nbx"] = nbx;
    dimensions_["nbu"] = nbu;
}

bool ocp_nlp::needs_initializing() { return needs_initializing_; }

void ocp_nlp::needs_initializing(bool flag) { needs_initializing_ = flag; }

vector<int> ocp_nlp::get_bound_indices(std::string bound, int stage)
{
    ocp_nlp_constraints_model **constraints = (ocp_nlp_constraints_model **) nlp_->constraints;
    ocp_nlp_constraints_dims **constraint_dims = (ocp_nlp_constraints_dims **) dims_->constraints;

    vector<int> idxb;

    int nbx = constraint_dims[stage]->nbx;
    int nbu = constraint_dims[stage]->nbu;
    int nu = dims_->nu[stage];

    if (bound == "lbx" || bound == "ubx")
    {
        for (int i = 0; i < nbx; ++i) idxb.push_back(constraints[stage]->idxb[nbu + i] - nu);
    }
    else if (bound == "lbu" || bound == "ubu")
    {
        for (int i = 0; i < nbu; ++i) idxb.push_back(constraints[stage]->idxb[i]);
    }
    else
    {
        throw std::invalid_argument("Expected 'lbx', 'ubx', 'lbu' or 'ubu' but got '" + bound +
                                    "'.");
    }

    return idxb;
}

void ocp_nlp::set_bound_indices(std::string bound, int stage, std::vector<int> idxb)
{
    ocp_nlp_constraints_model **constraints = (ocp_nlp_constraints_model **) nlp_->constraints;
    ocp_nlp_constraints_dims **constraint_dims = (ocp_nlp_constraints_dims **) dims_->constraints;

    int nbx = constraint_dims[stage]->nbx;
    int nbu = constraint_dims[stage]->nbu;

    if (bound == "x")
    {
        if (idxb.size() != (size_t) nbx)
            throw std::invalid_argument("Expected " + std::to_string(constraint_dims[stage]->nbx) +
                                        " bound indices but got " + std::to_string(idxb.size()) +
                                        ".");
        for (int i = 0; i < nbx; ++i)
            constraints[stage]->idxb[nbu + i] = dims_->nu[stage] + idxb[i];
    }
    else if (bound == "u")
    {
        if (idxb.size() != (size_t) nbu)
            throw std::invalid_argument("Expected " + std::to_string(constraint_dims[stage]->nbu) +
                                        " bound indices but got " + std::to_string(idxb.size()) +
                                        ".");
        for (int i = 0; i < nbu; ++i) constraints[stage]->idxb[i] = idxb[i];
    }
    else
    {
        throw std::invalid_argument("Can only get bound indices for 'x' and 'u'.");
    }
}

int ocp_nlp::num_stages() { return N; }

ocp_nlp::~ocp_nlp()
{
    for (int i = 0; i < N; ++i)
    {
        free(dims_->constraints[i]);
        free(dims_->cost[i]);
        free(dims_->dynamics[i]);

        free(nlp_->constraints[i]);
        free(nlp_->cost[i]);
        free(nlp_->dynamics[i]);
    }

    free(dims_->constraints[N]);
    free(dims_->cost[N]);

    free(nlp_->constraints[N]);
    free(nlp_->cost[N]);

    free(dims_->nv);
    free(dims_->nx);
    free(dims_->nu);
    free(dims_->ni);

    free(dims_->constraints);
    free(dims_->cost);
    free(dims_->dynamics);
    free(dims_->qp_solver);

    free(nlp_->constraints);
    free(nlp_->cost);
    free(nlp_->dynamics);

    free(config_->qp_solver);

    for (auto &elem : dynamics_handle_) dlclose(elem.second);
}

}  // namespace acados
