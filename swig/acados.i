/*
 *    This file is part of acados.
 *
 *    acados is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    acados is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with acados; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

%module acados

#if defined(SWIGMATLAB)
typedef mxArray LangObject;
#define NONE NULL
#elif defined(SWIGPYTHON)
typedef PyObject LangObject;
#define NONE Py_None
%{
#define SWIG_FILE_WITH_INIT
%}
#endif

%include "exception.i"

%exception {
  try {
    $action
  } catch (const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}

%include "conversions.i"
%include "ocp_typemaps.i"

%{
// TODO(dimitris): support compilation with visual studio
#if (defined _WIN32 || defined _WIN64 || defined __MINGW32__ || defined __MINGW64__)
#include <windows.h>
char compiler[16] = "gcc";
#else
#include <dlfcn.h>
char compiler[16] = "cc";
#endif
// #include <xmmintrin.h>  // for floating point exceptions
// _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);

#include <cstdlib>
#include <string>

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_condensing_hpipm.h"
#include "acados/ocp_qp/ocp_qp_condensing_qpoases.h"
#include "acados/ocp_qp/ocp_qp_hpipm.h"
#include "acados/ocp_qp/ocp_qp_hpmpc.h"
#ifdef OOQP
#include "acados/ocp_qp/ocp_qp_ooqp.h"
#endif
#include "acados/ocp_nlp/allocate_ocp_nlp.h"
#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/ocp_nlp/ocp_nlp_sm_common.h"
#include "acados/ocp_nlp/ocp_nlp_sm_gn.h"
#include "acados/ocp_nlp/ocp_nlp_sqp.h"
#include "acados/ocp_qp/ocp_qp_qpdunes.h"
#include "acados/sim/allocate_sim.h"
#include "acados/sim/sim_casadi_wrapper.h"
#include "acados/sim/sim_erk_integrator.h"
#include "acados/sim/sim_lifted_irk_integrator.h"
#include "acados/sim/sim_rk_common.h"
#include "acados/utils/casadi_wrapper.h"
#include "acados/utils/print.h"
#include "acados/utils/types.h"

#include "casadi/casadi.hpp"

%}

%{

int global_library_counter = 1;

// static bool is_valid_sim_dimensions_map(const LangObject *input) {
//     if (!is_map(input))
//         return false;
//     LangObject *nx = from(input, "nx");
//     LangObject *nu = from(input, "nu");
//     if (!is_integer(nx) || !is_integer(nu))
//         return false;
//     return true;
// }

std::string load_error_message() {
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

static int(*load_dims_function(void *handle, std::string name))(int_t *, int_t *,
                                                                 int_t *, int_t *) {
    #if (defined _WIN32 || defined _WIN64 || defined __MINGW32__ || defined __MINGW64__)
        return (int_t(*)(int_t*, int_t*, int_t*, int_t*))
                   GetProcAddress((HMODULE)handle, name.c_str());
    #else
        return (int_t(*)(int_t*, int_t*, int_t*, int_t*)) dlsym(handle, name.c_str());
    #endif
}

static const int_t *(*load_sparsity_function(void *handle, std::string name))(int_t) {
    #if (defined _WIN32 || defined _WIN64 || defined __MINGW32__ || defined __MINGW64__)
        return (const int_t *(*)(int_t)) GetProcAddress((HMODULE)handle, name.c_str());
    #else
        return (const int_t *(*)(int_t)) dlsym(handle, name.c_str());
    #endif
}

static void validate_model(casadi::Function& model) {
    if (model.n_in() != 2)
        throw std::runtime_error("An ODE model should have 2 inputs: state and controls");
    if (model.n_out() != 1)
        throw std::runtime_error("An ODE model should have 1 output: the right hand side");
    casadi::SX x = model.sx_in(0);
    casadi::SX u = model.sx_in(1);
    int_t nx = x.size1();
    const std::vector<casadi::SX> input = {x, u};
    casadi::SX rhs = casadi::SX::vertcat(model(input));
    if (rhs.size1() != nx)
        throw std::runtime_error("Length of right hand size should equal number of states");
}

static std::string generate_vde_function(casadi::Function& model) {
    validate_model(model);
    casadi::SX x = model.sx_in(0);
    casadi::SX u = model.sx_in(1);
    int_t nx = x.size1();
    int_t nu = u.size1();
    const std::vector<casadi::SX> states_controls = {x, u};
    casadi::SX rhs = casadi::SX::vertcat(model(states_controls));
    casadi::SX Sx = casadi::SX::sym("Sx", nx, nx);
    casadi::SX Su = casadi::SX::sym("Su", nx, nu);
    casadi::SX vde_x = casadi::SX::jtimes(rhs, x, Sx);
    casadi::SX vde_u = casadi::SX::jacobian(rhs, u) + casadi::SX::jtimes(rhs, x, Su);
    const std::vector<casadi::SX> input = {x, Sx, Su, u};
    const std::vector<casadi::SX> output = {rhs, vde_x, vde_u};
    std::string full_name = std::string("vde_") + model.name();
    std::string generated_file = full_name + std::string(".c");
    casadi::Function vde = casadi::Function(full_name, input, output);
    casadi::Dict opts;
    opts["with_header"] = casadi::GenericType(true);
    vde.generate(generated_file, opts);
    return full_name;
}

static std::string generate_jac_function(casadi::Function& model) {
    validate_model(model);
    casadi::SX x = model.sx_in(0);
    casadi::SX u = model.sx_in(1);
    const std::vector<casadi::SX> states_controls = {x, u};
    casadi::SX rhs = casadi::SX::vertcat(model(states_controls));
    int_t nx = x.size1();
    casadi::SX jac_x = casadi::SX::zeros(nx, nx) + casadi::SX::jacobian(rhs, x);
    const std::vector<casadi::SX> input = {x, u};
    const std::vector<casadi::SX> output = {rhs, jac_x};
    std::string full_name = std::string("jac_") + model.name();
    std::string generated_file = full_name + std::string(".c");
    casadi::Function jac = casadi::Function(full_name, input, output);
    casadi::Dict opts;
    opts["with_header"] = casadi::GenericType(true);
    jac.generate(generated_file, opts);
    return full_name;
}

enum generation_mode {
    GENERATE_VDE,
    GENERATE_JAC
};

static casadi_function_t compile_and_load(std::string name, void **handle) {
    std::string library_name = name + std::to_string(global_library_counter++) + std::string(".so");
    std::string path_to_library = library_name;
    char command[MAX_STR_LEN];
    snprintf(command, sizeof(command), "%s -fPIC -shared -g %s.c -o %s", compiler, name.c_str(),
        library_name.c_str());
    int compilation_failed = system(command);
    if (compilation_failed)
        throw std::runtime_error("Something went wrong when compiling the model.");
    #if (defined _WIN32 || defined _WIN64 || defined __MINGW32__ || defined __MINGW64__)
    *handle = LoadLibrary(path_to_library.c_str());
    #else
    *handle = dlopen(path_to_library.c_str(), RTLD_LAZY);
    #endif
    if (*handle == NULL)
        throw std::runtime_error("Loading of " + path_to_library + " failed. Error message: "
                                 + load_error_message());
    typedef int (*casadi_function_t)(const double** arg, double** res, int* iw, double* w, int mem);
    #if (defined _WIN32 || defined _WIN64 || defined __MINGW32__ || defined __MINGW64__)
    return (casadi_function_t) GetProcAddress((HMODULE)*handle, name.c_str());
    #else
    return (casadi_function_t) dlsym(*handle, name.c_str());
    #endif
}

void set_model(sim_in *sim, casadi::Function& f, double step, enum generation_mode mode) {
    std::string vde_name = generate_vde_function(f);
    void *vde_handle = malloc(sizeof(void *));
    sim->vde = compile_and_load(vde_name, &vde_handle);
    sim->jac = NULL;
    if (mode == GENERATE_JAC) {
        std::string jac_name = generate_jac_function(f);
        void *jac_handle = malloc(sizeof(void *));
        sim->jac = compile_and_load(jac_name, &jac_handle);
    }
    sim->forward_vde_wrapper = &vde_fun;
    sim->jacobian_wrapper = &jac_fun;
    sim->step = step;
}

LangObject *sim_output(const sim_in *in, const sim_out *out) {
    int_t x_dims[2] = {in->nx, 1};
    LangObject *x_final = new_matrix(x_dims, out->xn);
    int_t S_dims[2] = {in->num_forw_sens, in->nx+in->nu};
    LangObject *S_forward = new_matrix(S_dims, out->S_forw);
    return new_sim_output_tuple(x_final, S_forward);
}

%}

%ignore sim_in;
%ignore sim_out;
%include "acados/sim/sim_common.h"

%extend sim_solver {
    sim_solver(const char *solver_name, casadi::Function& model, LangObject *options = NONE) {

        const char *fieldnames[2] = {"time_step", "order"};
        real_t time_step = 0.1;
        int_t order = 4;
        if (options != NONE) {
            if (has(options, fieldnames[0]))
                time_step = real_from(options, fieldnames[0]);
            if (has(options, fieldnames[1]))
                order = int_from(options, fieldnames[1]);
        }
        sim_solver *solver = (sim_solver *) malloc(sizeof(sim_solver));
        validate_model(model);
        int_t nx = model.sx_in(0).size1();
        int_t nu = model.sx_in(1).size1();
        sim_in *input = (sim_in *) malloc(sizeof(sim_in));
        input->nx = nx;
        input->nu = nu;
        input->sens_forw = true;
        input->sens_adj = false;
        input->sens_hess = false;
        input->num_forw_sens = nx + nu;
        input->num_steps = 10;
        allocate_sim_in(input, order);

        sim_out *output = (sim_out *) malloc(sizeof(sim_out));
        allocate_sim_out(input, output);

        void *args = NULL;
        void *memory = NULL;
        int_t workspace_size;
        void *workspace = NULL;
        if (!strcmp("explicit runge-kutta", solver_name) || !strcmp("erk", solver_name) ||
                                                            !strcmp("rk", solver_name)) {
            solver->fun = sim_erk;
            args = (void *) malloc(sizeof(sim_rk_opts));
            sim_erk_create_arguments(args, order);
            workspace_size = sim_erk_calculate_workspace_size(input, args);
            workspace = malloc(workspace_size);
            set_model(input, model, time_step / input->num_steps, GENERATE_VDE);
        } else if (!strcmp("implicit runge-kutta", solver_name) || !strcmp("irk", solver_name)
                   || !strcmp("in", solver_name) || !strcmp("inis", solver_name)) {
            solver->fun = sim_lifted_irk;
            args = (void *) malloc(sizeof(sim_rk_opts));
            sim_irk_create_arguments(args, 2, "Gauss");
            if (!strcmp("implicit runge-kutta", solver_name)
                || !strcmp("irk", solver_name))
                sim_irk_create_Newton_scheme(args, 2, "Gauss", exact);
            if (!strcmp("in", solver_name))
                sim_irk_create_Newton_scheme(args, 2, "Gauss", simplified_in);
            if (!strcmp("inis", solver_name))
                sim_irk_create_Newton_scheme(args, 2, "Gauss", simplified_inis);
            memory = (sim_lifted_irk_memory *) malloc(sizeof(sim_lifted_irk_memory));
            sim_lifted_irk_create_memory(input, args, (sim_lifted_irk_memory *) memory);
            workspace_size = sim_lifted_irk_calculate_workspace_size(input, args);
            workspace = (void *) malloc(workspace_size);
            set_model(input, model, time_step / input->num_steps, GENERATE_JAC);
        } else {
            throw std::invalid_argument("Integrator name not known!");
        }
        solver->in = input;
        solver->out = output;
        solver->args = args;
        solver->mem = memory;
        solver->work = workspace;
        return solver;
    }

    LangObject *evaluate(LangObject *initial_state, LangObject *control) {
        fill_array_from(initial_state, $self->in->x, $self->in->nx);
        fill_array_from(control, $self->in->u, $self->in->nu);
        int_t return_code = $self->fun($self->in, $self->out, $self->args, $self->mem, $self->work);
        if (return_code != 0) {
            throw std::runtime_error("Integrator failed!");
        }
        return sim_output($self->in, $self->out);
    }
}

%{

static bool is_valid_ocp_dimensions_map(const LangObject *input) {
    if (!is_map(input))
        return false;
    int_t N = int_from(input, "N");
    LangObject *nx = from(input, "nx");
    if (!is_integer(nx) && !is_sequence(nx, N+1)) {
        return false;
    }
    LangObject *nu = from(input, "nu");
    if (!is_integer(nu) && !is_sequence(nu, N)) {
        return false;
    }
    return true;
}

// static bool qp_dimensions_equal(const ocp_qp_in *qp1, const ocp_qp_in *qp2) {
//     if (qp1->N != qp2->N)
//         return false;
//     int_t N = qp1->N;
//     for (int_t i = 0; i < N; i++) {
//         if (qp1->nx[i] != qp2->nx[i])
//             return false;
//         else if (qp1->nu[i] != qp2->nu[i])
//             return false;
//         else if (qp1->nb[i] != qp2->nb[i])
//             return false;
//         else if (qp1->nc[i] != qp2->nc[i])
//             return false;
//     }
//     if (qp1->nx[N] != qp2->nx[N])
//         return false;
//     else if (qp1->nb[N] != qp2->nb[N])
//         return false;
//     else if (qp1->nc[N] != qp2->nc[N])
//         return false;
//     return true;
// }

LangObject *ocp_qp_output(const ocp_qp_in *in, const ocp_qp_out *out) {
    const real_t **states_copy, **controls_copy;
    states_copy = (const real_t **) malloc((in->N+1) * sizeof(real_t *));
    controls_copy = (const real_t **) malloc((in->N+1) * sizeof(real_t *));
    for (int_t i = 0; i <= in->N; i++) {
        states_copy[i] = (const real_t *) calloc(in->nx[i], sizeof(real_t));
        memcpy((void *)states_copy[i], (void *)out->x[i], in->nx[i] * sizeof(real_t));
        controls_copy[i] = (const real_t *) calloc(in->nu[i], sizeof(real_t));
        memcpy((void *)controls_copy[i], (void *)out->u[i], in->nu[i] * sizeof(real_t));
    }

    LangObject *x_star = new_sequence_from(states_copy, in->N+1, in->nx);
    LangObject *u_star = new_sequence_from(controls_copy, in->N+1, in->nu);
    return new_ocp_output_tuple(x_star, u_star);
}

%}

%ignore ocp_qp_in_calculate_size;
%ignore assign_ocp_qp_in;
%ignore create_ocp_qp_in;
%ignore ocp_qp_out_calculate_size;
%ignore assign_ocp_qp_out;
%ignore create_ocp_qp_out;
%ignore ocp_qp_out;
%ignore create_ocp_qp_solver;
%rename(ocp_qp) ocp_qp_in;
%include "acados/ocp_qp/ocp_qp_common.h"

%extend ocp_qp_in {
#if defined(SWIGMATLAB)
    %matlabcode %{
    function self = subsasgn(self, s, v)
      if numel(s) == 1 && strcmp(s.type, '.')
        self.(s.subs)(v)
      elseif numel(s) == 2 && strcmp(s(1, 1).type, '.') && strcmp(s(1, 2).type, '{}')
        for cell_no = 1:numel(s(1, 2).subs)
          index_group = s(1, 2).subs{cell_no};
          for index = index_group
            cell_array = self.(s(1, 1).subs)();
            cell_array{index} = v;
            self.(s(1, 1).subs)(cell_array)
          end
        end
      else
        self = builtin('subsasgn', self, s, v);
      end
    end
    %}
#endif
    ocp_qp_in(LangObject *input_map) {
        if (!is_valid_ocp_dimensions_map(input_map)) {
            char err_msg[MAX_STR_LEN];
            snprintf(err_msg, sizeof(err_msg), "Input must be a valid OCP %s that specifies at "
                "least N, nx, nu", LANG_MAP_NAME);
            throw std::invalid_argument(err_msg);
        }
        int_t N = int_from(input_map, "N");
        int_t nx[N+1], nu[N+1], nb[N+1], nc[N+1];
        fill_array_from(input_map, "nx", nx, N+1);
        fill_array_from(input_map, "nu", nu, N+1);
        fill_array_from(input_map, "nb", nb, N+1);
        fill_array_from(input_map, "nc", nc, N+1);
        nu[N] = 0;
        // Default behavior is that initial state is fixed
        if (!has(input_map, "nb")) {
            nb[0] = nx[0];
        }
        ocp_qp_in *qp_in = create_ocp_qp_in(N, nx, nu, nb, nc);
        // Initial state is fixed
        if (!has(input_map, "nb")) {
            int idxb[nb[0]];
            for (int_t i = 0; i < nb[0]; i++)
                idxb[i] = i;
            memcpy((void *) qp_in->idxb[0], idxb, sizeof(idxb));
        }
        return qp_in;
    }
}

%extend ocp_qp_solver {
    ocp_qp_solver(const char *solver_name, ocp_qp_in *qp_in, LangObject *options = NONE) {
        ocp_qp_solver *solver = create_ocp_qp_solver(qp_in, solver_name, NULL);
        return solver;
    }

    LangObject *evaluate() {
        int_t return_code = $self->fun($self->qp_in, $self->qp_out, $self->args, \
            $self->mem, $self->work);
        if (return_code != 0) {
            throw std::runtime_error("qp solver failed!");
        }
        return ocp_qp_output($self->qp_in, $self->qp_out);
    }

    LangObject *evaluate(LangObject *x0) {
        fill_array_from(x0, (real_t **) $self->qp_in->lb, 1, $self->qp_in->nx);
        fill_array_from(x0, (real_t **) $self->qp_in->ub, 1, $self->qp_in->nx);
        int_t return_code = $self->fun($self->qp_in, $self->qp_out, $self->args, \
            $self->mem, $self->work);
        if (return_code != 0)
            throw std::runtime_error("qp solver failed!");
        return ocp_qp_output($self->qp_in, $self->qp_out);
    }
}

%typemap(in) real_t ** ls_cost_matrix {
    $1 = arg1->W;
    int_t W_dimensions[arg1->N+1];
    for (int_t i = 0; i < arg1->N+1; i++) {
        W_dimensions[i] = arg1->fun[i]->ny;
    }
    fill_array_from($input, $1, arg1->N+1, W_dimensions, W_dimensions);
}

%typemap(out) real_t ** ls_cost_matrix {
    int_t W_dimensions[arg1->N+1];
    for (int_t i = 0; i < arg1->N+1; i++) {
        W_dimensions[i] = arg1->fun[i]->ny;
    }
    $result = new_sequence_from((const real_t **) $1, arg1->N+1, W_dimensions, W_dimensions);
}

%typemap(in) real_t **ls_cost_ref {
    $1 = arg1->y_ref;
    int_t y_dimensions[arg1->N + 1];
    int_t ones[arg1->N + 1];
    for (int_t i = 0; i < arg1->N + 1; i++) {
        y_dimensions[i] = arg1->fun[i]->ny;
        ones[i] = 1;
    }
    fill_array_from($input, $1, arg1->N + 1, y_dimensions, ones);
}

%typemap(out) real_t **ls_cost_ref {
    int_t y_dimensions[arg1->N + 1];
    int_t ones[arg1->N + 1];
    for (int_t i = 0; i < arg1->N + 1; i++) {
        y_dimensions[i] = arg1->fun[i]->ny;
        ones[i] = 1;
    }
    $result = new_sequence_from((const real_t **)$1, arg1->N + 1, y_dimensions, ones);
}

%{
LangObject *ocp_nlp_output(const ocp_nlp_in *in, const ocp_nlp_out *out) {
    LangObject *x_star = new_sequence_from((const real_t **) out->x, in->N+1, in->nx);
    LangObject *u_star = new_sequence_from((const real_t **) out->u, in->N+1, in->nu);
    return new_ocp_output_tuple(x_star, u_star);
}

void ocp_nlp_ls_cost_ls_cost_matrix_set(ocp_nlp_ls_cost *ls_cost, real_t **matrix) {
    ls_cost->W = matrix;
}

real_t **ocp_nlp_ls_cost_ls_cost_matrix_get(ocp_nlp_ls_cost *ls_cost) {
    return ls_cost->W;
}

void ocp_nlp_ls_cost_ls_cost_ref_set(ocp_nlp_ls_cost *ls_cost, real_t **y_ref) {
    ls_cost->y_ref = y_ref;
}

real_t **ocp_nlp_ls_cost_ls_cost_ref_get(ocp_nlp_ls_cost *ls_cost) {
    return ls_cost->y_ref;
}
%}

%ignore ocp_nlp_memory;
%ignore ocp_nlp_out;
%ignore ocp_nlp_calculate_memory_size;
%ignore ocp_nlp_assign_memory;
%ignore ocp_nlp_create_memory;
%ignore ocp_nlp_destroy;
%rename(ocp_nlp) ocp_nlp_in;
%include "acados/ocp_nlp/ocp_nlp_common.h"

%ignore ocp_nlp_sm;
%ignore ocp_nlp_sm_in;
%ignore ocp_nlp_sm_out;
%include "acados/ocp_nlp/ocp_nlp_sm_common.h"

%ignore ocp_nlp_sm_gn_args;
%ignore ocp_nlp_sm_gn_memory;
%ignore ocp_nlp_sm_gn_workspace;
%ignore ocp_nlp_sm_gn_create_arguments;
%ignore ocp_nlp_sm_gn_calculate_memory_size;
%ignore ocp_nlp_sm_gn_assign_memory;
%ignore ocp_nlp_sm_gn_create_memory;
%ignore ocp_nlp_sm_gn_calculate_workspace_size;
%ignore ocp_nlp_sm_gn_assign_workspace;
%ignore ocp_nlp_sm_gn_create_workspace;
%ignore ocp_nlp_sm_gn_initialize;
%ignore ocp_nlp_sm_gn_destroy;
%include "acados/ocp_nlp/ocp_nlp_sm_gn.h"

%extend ocp_nlp_function {
    ocp_nlp_function(casadi::Function& cas_fun, LangObject *options = NONE) {
        casadi::SX x = cas_fun.sx_in(0);
        casadi::SX u = cas_fun.sx_in(1);
        casadi::SX y = cas_fun.sx_out(0);
        int_t nx = x.size1();
        int_t nu = u.size1();
        int_t ny = y.size1();

        // Derive sensitivities
        const std::vector<casadi::SX> input = {x, u};
        casadi::SX xu = casadi::SX::vertcat(input);
        casadi::SX rhs = casadi::SX::vertcat(cas_fun(input));
        casadi::SX jac = casadi::SX::jacobian(rhs, xu);
        casadi::SX hess = casadi::SX::jacobian(jac, xu);

        // Generate code
        const std::vector<casadi::SX> input_vector = {x, u};
        const std::vector<casadi::SX> output_vector = {rhs, jac, hess};
        std::string full_name = cas_fun.name();
        std::string generated_file = full_name + std::string(".c");
        casadi::Function extended_function =
            casadi::Function(full_name, input_vector, output_vector);
        extended_function.generate(generated_file);

        void *handle = malloc(sizeof(void *));
        casadi_function_t generated_function = compile_and_load(full_name, &handle);

        // Create nlp_function object
        ocp_nlp_function *nlp_function = (ocp_nlp_function *)malloc(sizeof(ocp_nlp_function));

        // Initialize dimensions
        nlp_function->nx = nx;
        nlp_function->nu = nu;
        nlp_function->np = 0;
        nlp_function->ny = ny;

        // Initialize casadi wrapper input
        nlp_function->in = (casadi_wrapper_in *)malloc(sizeof(casadi_wrapper_in));
        nlp_function->in->compute_jac = true;
        nlp_function->in->compute_hess = true;

        // Initialize casadi wrapper output
        nlp_function->out = (casadi_wrapper_out *)malloc(sizeof(casadi_wrapper_out));

        // Write casadi wrapper arguments
        nlp_function->args = casadi_wrapper_create_arguments();
        nlp_function->args->fun = generated_function;
        nlp_function->args->dims = load_dims_function(handle, full_name + "_work");
        nlp_function->args->sparsity = load_sparsity_function(handle, full_name + "_sparsity_out");
        casadi_wrapper_initialize(nlp_function->in, nlp_function->args,
                                  &nlp_function->work);

        return nlp_function;
    }

    // TODO(nielsvd): write destructor

    LangObject *evaluate(LangObject *x, LangObject *u) {
        int_t nx = $self->nx;
        int_t nu = $self->nu;
        int_t ny = $self->ny;

        real_t *x_in = (real_t *)malloc(nx * sizeof(real_t));
        real_t *u_in = (real_t *)malloc(nu * sizeof(real_t));
        real_t *p_in = NULL;

        $self->in->x = x_in;
        $self->in->u = u_in;
        $self->in->p = p_in;

        real_t *y_out = (real_t *)malloc(ny * sizeof(real_t));
        real_t *jac_y_out = (real_t *)malloc(ny * (nx + nu) * sizeof(real_t));
        real_t *hess_y_out = (real_t *)malloc(ny * (nx + nu) * (nx + nu) * sizeof(real_t));

        $self->out->y = y_out;
        $self->out->jac_y = jac_y_out;
        $self->out->hess_y = hess_y_out;

        fill_array_from(x, x_in, nx);
        fill_array_from(u, u_in, nu);

        int_t return_code = casadi_wrapper($self->in, $self->out, $self->args, $self->work);
        if (return_code != 0) {
            throw std::runtime_error("Function evaluation failed!");
        }

        int_t y_dims[2] = {ny, 1};
        LangObject *y = new_matrix(y_dims, y_out);
        int_t jac_dims[2] = {ny, nx+nu};
        LangObject *jac_y = new_matrix(jac_dims, jac_y_out);
        int_t hess_dims[2] = {ny * (nx + nu), nx + nu};
        LangObject *hess_y = new_matrix(hess_dims, hess_y_out);

        free(x_in);
        free(u_in);

        return new_ocp_nlp_function_output_tuple(y, jac_y, hess_y);
    }
}

%extend ocp_nlp_ls_cost{
#if defined(SWIGMATLAB)
    %matlabcode %{
    function self = subsasgn(self, s, v)
      if numel(s) == 1 && strcmp(s.type, '.')
        self.(s.subs)(v)
      elseif numel(s) == 2 && strcmp(s(1, 1).type, '.') && strcmp(s(1, 2).type, '{}')
        for cell_no = 1:numel(s(1, 2).subs)
          index_group = s(1, 2).subs{cell_no};
          for index = index_group
            cell_array = self.(s(1, 1).subs)();
            cell_array{index} = v;
            self.(s(1, 1).subs)(cell_array)
          end
        end
      else
        self = builtin('subsasgn', self, s, v);
      end
    end
    %}
#endif
    real_t **ls_cost_matrix;
    real_t **ls_cost_ref;
    ocp_nlp_ls_cost(LangObject *N, LangObject* stage_costs, LangObject *options = NONE){
        int_t NN = int_from(N);
        if (!is_sequence(stage_costs, NN + 1)) {
            throw std::runtime_error(
                "The array of stage cost functions must be of length N+1!");
        }

        ocp_nlp_ls_cost *ls_cost = (ocp_nlp_ls_cost *) malloc(sizeof(ocp_nlp_ls_cost));
        ls_cost->N = NN;

        // Read NLS cost output functions
        ls_cost->fun = (ocp_nlp_function **) malloc((NN+1)*sizeof(ocp_nlp_function *));
        for (int_t i = 0; i <= NN; i++) {
            LangObject* stage_cost_lo = from(stage_costs, i);
            ocp_nlp_function *stage_cost;
            SWIG_ConvertPtr(stage_cost_lo, (void **)&stage_cost, SWIGTYPE_p_ocp_nlp_function, 0);
            // Initialize LS cost
            ls_cost->fun[i] = (ocp_nlp_function *) malloc(sizeof(ocp_nlp_function));
            ls_cost->fun[i]->nx = stage_cost->nx;
            ls_cost->fun[i]->nu = stage_cost->nu;
            ls_cost->fun[i]->np = stage_cost->np;
            ls_cost->fun[i]->ny = stage_cost->ny;
            ls_cost->fun[i]->in = (casadi_wrapper_in *)malloc(sizeof(casadi_wrapper_in));
            ls_cost->fun[i]->in->compute_jac = true;
            ls_cost->fun[i]->in->compute_hess = false;
            ls_cost->fun[i]->out = (casadi_wrapper_out *)malloc(sizeof(casadi_wrapper_out));
            ls_cost->fun[i]->args = casadi_wrapper_create_arguments();
            ls_cost->fun[i]->args->fun = stage_cost->args->fun;
            ls_cost->fun[i]->args->dims = stage_cost->args->dims;
            ls_cost->fun[i]->args->sparsity = stage_cost->args->sparsity;
            casadi_wrapper_initialize(ls_cost->fun[i]->in, ls_cost->fun[i]->args,
                                      &ls_cost->fun[i]->work);
        }

        // Prepare memory for cost and reference
        ls_cost->W = (real_t **) malloc((NN+1)*sizeof(real_t *));
        ls_cost->y_ref = (real_t **) malloc((NN+1)*sizeof(real_t *));

        for (int_t i = 0; i <= NN; i++) {
            int_t ny = ls_cost->fun[i]->ny;
            ls_cost->W[i] = (real_t *) malloc(ny*ny*sizeof(real_t));
            ls_cost->y_ref[i] = (real_t *) malloc(ny*sizeof(real_t));
            for (int_t j = 0; j < ny; j++) {
                ls_cost->y_ref[i][j] = 0;
                for (int_t k = 0; k < ny; k++) ls_cost->W[i][j*ny+k] = 0;
            }
        }

        return ls_cost;
    }
}

%extend ocp_nlp_in {
#if defined(SWIGMATLAB)
    %matlabcode %{
    function self = subsasgn(self, s, v)
      if numel(s) == 1 && strcmp(s.type, '.')
        self.(s.subs)(v)
      elseif numel(s) == 2 && strcmp(s(1, 1).type, '.') && strcmp(s(1, 2).type, '{}')
        for cell_no = 1:numel(s(1, 2).subs)
          index_group = s(1, 2).subs{cell_no};
          for index = index_group
            cell_array = self.(s(1, 1).subs)();
            cell_array{index} = v;
            self.(s(1, 1).subs)(cell_array)
          end
        end
      else
        self = builtin('subsasgn', self, s, v);
      end
    end
    %}
#endif
    ocp_nlp_in(LangObject *input_map) {
        ocp_nlp_in *nlp_in = (ocp_nlp_in *) malloc(sizeof(ocp_nlp_in));
        if (!is_valid_ocp_dimensions_map(input_map)) {
            char err_msg[MAX_STR_LEN];
            snprintf(err_msg, sizeof(err_msg), "Input must be a valid OCP %s that specifies at "
                "least N, nx, nu", LANG_MAP_NAME);
            throw std::invalid_argument(err_msg);
        }
        int_t N = int_from(input_map, "N");
        int_t nx[N+1], nu[N+1], nb[N+1], ng[N+1];
        fill_array_from(input_map, "nx", nx, N+1);
        fill_array_from(input_map, "nu", nu, N+1);
        fill_array_from(input_map, "nb", nb, N+1);
        fill_array_from(input_map, "ng", ng, N+1);
        nu[N] = 0;
        // Default behavior is that initial state is fixed
        if (!has(input_map, "nb")) {
            nb[0] = nx[0];
        }
        allocate_ocp_nlp_in(N, nx, nu, nb, ng, 0, nlp_in);
        if (!has(input_map, "nb")) {
            int idxb[nb[0]];
            for (int_t i = 0; i < nb[0]; i++)
                idxb[i] = i;
            memcpy((void *) nlp_in->idxb[0], idxb, sizeof(idxb));
        }

        // TODO(nielsvd): cost, path_constraints
        return nlp_in;
    }

    // TODO(nielsvd): write destructor

    void set_model(casadi::Function& f, double step) {
        std::string model_name = generate_vde_function(f);
        void *handle = malloc(sizeof(void *));
        casadi_function_t eval = compile_and_load(model_name, &handle);
        sim_solver **simulators = (sim_solver **)$self->sim;
        for (int_t i = 0; i < $self->N; i++) {
            simulators[i]->in->vde = eval;
            simulators[i]->in->forward_vde_wrapper = &vde_fun;
            simulators[i]->in->step = step;
        }
    }

    void set_cost(ocp_nlp_ls_cost *ls_cost) {
        // Make compatible with general cost
        $self->cost = (void *) ls_cost;
    }

    void set_path_constraints(LangObject *path_constraints) {
        int_t NN = $self->N;
        $self->path_constraints = (void **) malloc((NN + 1) * sizeof(ocp_nlp_function *));
        ocp_nlp_function **pathcons =(ocp_nlp_function **)$self->path_constraints;
        for (int_t i = 0; i <= NN; i++) {
            LangObject *path_constraint_lo = from(path_constraints, i);
            ocp_nlp_function *path_constraint;
            SWIG_ConvertPtr(path_constraint_lo, (void **)&path_constraint,
                            SWIGTYPE_p_ocp_nlp_function, 0);
            // Initialize path constraint
            pathcons[i] = (ocp_nlp_function *) malloc(sizeof(ocp_nlp_function));
            pathcons[i]->nx = path_constraint->nx;
            pathcons[i]->nu = path_constraint->nu;
            pathcons[i]->np = path_constraint->np;
            pathcons[i]->ny = path_constraint->ny;
            pathcons[i]->in =
                (casadi_wrapper_in *)malloc(sizeof(casadi_wrapper_in));
            pathcons[i]->in->compute_jac = true;
            pathcons[i]->in->compute_hess = false;
            pathcons[i]->out =
                (casadi_wrapper_out *)malloc(sizeof(casadi_wrapper_out));
            pathcons[i]->args = casadi_wrapper_create_arguments();
            pathcons[i]->args->fun = path_constraint->args->fun;
            pathcons[i]->args->dims = path_constraint->args->dims;
            pathcons[i]->args->sparsity = path_constraint->args->sparsity;
            casadi_wrapper_initialize(pathcons[i]->in, pathcons[i]->args,
                                      &pathcons[i]->work);
        }
    }
}

%extend ocp_nlp_solver {
    ocp_nlp_solver(const char *solver_name, ocp_nlp_in *nlp_in, LangObject *options = NONE) {
        const char *fieldnames[4] = {"qp_solver", "sensitivity_method",
                                     "integrator_steps", "SQP_steps"};
        const char *qp_solver = "qpdunes";
        const char *sensitivity_method = "gauss-newton";
        int_t integrator_steps = 1;
        int_t sqp_steps = 1;
        if (options == NONE) {
#if defined(SWIGMATLAB)
            const mwSize dims[1] = {(const mwSize) 1};
            options = mxCreateStructArray(1, dims, 4, fieldnames);
            mxArray *qp_solver_string = mxCreateString(qp_solver);
            mxArray *sensitivity_method_string = mxCreateString(sensitivity_method);
            mxArray *integrator_steps_array = mxCreateDoubleScalar(integrator_steps);
            mxArray *sqp_steps_array = mxCreateDoubleScalar(sqp_steps);
            mxSetField(options, 0, fieldnames[0], qp_solver_string);
            mxSetField(options, 0, fieldnames[1], sensitivity_method_string);
            mxSetField(options, 0, fieldnames[2], integrator_steps_array);
            mxSetField(options, 0, fieldnames[3], sqp_steps_array);
#elif defined(SWIGPYTHON)
            options = PyDict_New();
            PyDict_SetItem(options,
                           PyString_FromString(fieldnames[0]),
                           PyString_FromString(qp_solver));
            PyDict_SetItem(options,
                           PyString_FromString(fieldnames[1]),
                           PyString_FromString(sensitivity_method));
            PyDict_SetItem(options,
                           PyString_FromString(fieldnames[2]),
                           PyLong_FromLong((long) integrator_steps));
            PyDict_SetItem(options,
                           PyString_FromString(fieldnames[3]),
                           PyLong_FromLong((long) sqp_steps));
#endif
        }
        if (!is_map(options)) {
            char err_msg[MAX_STR_LEN];
            snprintf(err_msg, sizeof(err_msg), "Please pass options as a %s", LANG_MAP_NAME);
            throw std::runtime_error(err_msg);
        }
        if (has(options, fieldnames[0]))
            qp_solver = char_from(options, fieldnames[0]);
        if (has(options, fieldnames[1]))
            sensitivity_method = char_from(options, fieldnames[1]);
        if (has(options, fieldnames[2]))
            integrator_steps = int_from(options, fieldnames[2]);
        if (has(options, fieldnames[3]))
            sqp_steps = int_from(options, fieldnames[3]);
        ocp_nlp_solver *solver = (ocp_nlp_solver *) malloc(sizeof(ocp_nlp_solver));
        void *args = NULL;
        void *mem = NULL;
        void *workspace = NULL;
        int_t N = nlp_in->N;

        if (!strcmp("sqp", solver_name)) {
            solver->fun = &ocp_nlp_sqp;

            args = (ocp_nlp_sqp_args *) malloc(sizeof(ocp_nlp_sqp_args));

            // Select QP solver based on user input
            ((ocp_nlp_sqp_args *)args)->qp_solver = (ocp_qp_solver *) malloc(sizeof(ocp_qp_solver));
            ocp_qp_solver *qpsol = ((ocp_nlp_sqp_args *)args)->qp_solver;
            qpsol->qp_in = create_ocp_qp_in(
                nlp_in->N, nlp_in->nx, nlp_in->nu, nlp_in->nb, nlp_in->ng);
            qpsol->qp_out = create_ocp_qp_out(
                nlp_in->N, nlp_in->nx, nlp_in->nu, nlp_in->nb, nlp_in->ng);
            // TODO(nielsvd): lines below should go
            int_t **idxb = (int_t **)qpsol->qp_in->idxb;
            for (int_t i = 0; i <= N; i++)
                for (int_t j = 0; j < nlp_in->nb[i]; j++)
                    idxb[i][j] = nlp_in->idxb[i][j];
            if (!strcmp(qp_solver, "qpdunes")) {
                qpsol->args = (void *)ocp_qp_qpdunes_create_arguments(QPDUNES_NONLINEAR_MPC);
                qpsol->fun = &ocp_qp_qpdunes;
                qpsol->initialize = &ocp_qp_qpdunes_initialize;
                qpsol->destroy = &ocp_qp_qpdunes_destroy;
#ifdef OOQP
            } else if (!strcmp(qp_solver, "ooqp")) {
                qpsol->args = (void *)ocp_qp_ooqp_create_arguments();
                qpsol->fun = &ocp_qp_ooqp;
                qpsol->initialize = &ocp_qp_ooqp_initialize;
                qpsol->destroy = &ocp_qp_ooqp_destroy;
#endif
            } else if (!strcmp(qp_solver, "condensing_qpoases")) {
                qpsol->args =
                    (void *)ocp_qp_condensing_qpoases_create_arguments(
                        qpsol->qp_in);
                qpsol->fun = &ocp_qp_condensing_qpoases;
                qpsol->initialize = &ocp_qp_condensing_qpoases_initialize;
                qpsol->destroy = &ocp_qp_condensing_qpoases_destroy;
            } else if (!strcmp(qp_solver, "hpmpc")) {
                qpsol->args = (void *)ocp_qp_hpmpc_create_arguments(
                    qpsol->qp_in, HPMPC_DEFAULT_ARGUMENTS);
                qpsol->fun = &ocp_qp_hpmpc;
                qpsol->initialize = &ocp_qp_hpmpc_initialize;
                qpsol->destroy = &ocp_qp_hpmpc_destroy;
            } else if (!strcmp(qp_solver, "condensing_hpipm")) {
                qpsol->args =
                    ocp_qp_condensing_hpipm_create_arguments(qpsol->qp_in);
                qpsol->fun = &ocp_qp_condensing_hpipm;
                qpsol->initialize = &ocp_qp_condensing_hpipm_initialize;
                qpsol->destroy = &ocp_qp_condensing_hpipm_destroy;
            } else if (!strcmp(qp_solver, "hpipm")) {
                qpsol->args = ocp_qp_hpipm_create_arguments(qpsol->qp_in);
                qpsol->fun = &ocp_qp_hpipm;
                qpsol->initialize = &ocp_qp_hpipm_initialize;
                qpsol->destroy = &ocp_qp_hpipm_destroy;
            } else {
                throw std::invalid_argument("Chosen QP solver not available!");
            }

            // Select sensitivity method based on user input
            ((ocp_nlp_sqp_args *)args)->sensitivity_method =
                (ocp_nlp_sm *)malloc(sizeof(ocp_nlp_sm));
            ocp_nlp_sm *sm = ((ocp_nlp_sqp_args *)args)->sensitivity_method;
            if (!strcmp(sensitivity_method, "gauss-newton")) {
                sm->fun = &ocp_nlp_sm_gn;
                sm->initialize = &ocp_nlp_sm_gn_initialize;
                sm->destroy = &ocp_nlp_sm_gn_destroy;
            } else {
                throw std::invalid_argument(
                    "Chosen sensitivity method not available!");
            }

            ((ocp_nlp_sqp_args *) args)->maxIter = sqp_steps;
            sim_solver** simulators = (sim_solver**) nlp_in->sim;
            for (int_t i = 0; i < N; i++) {
                simulators[i]->in->nx = nlp_in->nx[i];
                simulators[i]->in->nu = nlp_in->nu[i];
                simulators[i]->in->sens_forw = true;
                simulators[i]->in->sens_adj = false;
                simulators[i]->in->sens_hess = false;
                simulators[i]->in->num_forw_sens =
                    nlp_in->nx[i] + nlp_in->nu[i];
                simulators[i]->in->num_steps = integrator_steps;
                simulators[i]->args = (void *)malloc(sizeof(sim_rk_opts));
                sim_erk_create_arguments(simulators[i]->args, 4);
                int_t erk_workspace_size = sim_erk_calculate_workspace_size(
                    simulators[i]->in, simulators[i]->args);
                simulators[i]->work = (void *)malloc(erk_workspace_size);
                simulators[i]->fun = &sim_erk;
            }

            ocp_nlp_sqp_initialize(nlp_in, args, &mem, &workspace);

            // Set initial guess to zero
            ocp_nlp_sqp_memory *sqp_mem = (ocp_nlp_sqp_memory *) mem;
            real_t **x = (real_t **) sqp_mem->common->x;
            real_t **u = (real_t **) sqp_mem->common->u;
            real_t **pi = (real_t **) sqp_mem->common->pi;
            real_t **lam = (real_t **) sqp_mem->common->lam;
            for (int_t i = 0; i <= N; i++) {
                for (int_t j = 0; j < nlp_in->nx[i]; j++) x[i][j] = 0.0;
                for (int_t j = 0; j < nlp_in->nu[i]; j++) u[i][j] = 0.0;
                for (int_t j = 0; j < 2 * nlp_in->nb[i] + 2 * nlp_in->ng[i]; j++) {
                    lam[i][j] = 0.0;
                }
            }
            for (int_t i = 0; i < N; i++) {
                for (int_t j = 0; j < nlp_in->nx[i+1]; j++) pi[i][j] = 0.0;
            }
        } else {
            throw std::invalid_argument("Solver name not known!");
        }

        solver->nlp_in = nlp_in;
        ocp_nlp_out *nlp_out = (ocp_nlp_out *) malloc(sizeof(ocp_nlp_out));
        allocate_ocp_nlp_out(nlp_in, nlp_out);
        solver->nlp_out = nlp_out;
        solver->args = args;
        solver->mem = mem;
        solver->work = workspace;
        return solver;
    }

    // TODO(nielsvd): write destructor

    LangObject *evaluate() {
        int_t fail = $self->fun($self->nlp_in, $self->nlp_out, $self->args, $self->mem,
                                $self->work);
        if (fail)
            throw std::runtime_error("nlp solver failed with error code " + std::to_string(fail));
        return ocp_nlp_output($self->nlp_in, $self->nlp_out);
    }

    LangObject *evaluate(LangObject *x0) {
        fill_array_from(x0, (real_t **)$self->nlp_in->lb, 1, $self->nlp_in->nx);
        fill_array_from(x0, (real_t **)$self->nlp_in->ub, 1, $self->nlp_in->nx);
        int_t fail = $self->fun($self->nlp_in, $self->nlp_out, $self->args,
                                $self->mem, $self->work);

        if (fail)
            throw std::runtime_error("nlp solver failed with error code " + std::to_string(fail));
        return ocp_nlp_output($self->nlp_in, $self->nlp_out);
    }
}
