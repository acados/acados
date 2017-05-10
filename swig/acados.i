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
#elif defined(SWIGPYTHON)
typedef PyObject LangObject;
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
#include <dlfcn.h>
// #include <xmmintrin.h>  // for floating point exceptions

#include <cstdlib>
#include <string>

#include "acados/ocp_qp/allocate_ocp_qp.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_condensing_qpoases.h"
#include "acados/ocp_qp/ocp_qp_ooqp.h"
#include "acados/ocp_qp/ocp_qp_qpdunes.h"
#include "acados/ocp_nlp/allocate_ocp_nlp.h"
#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/ocp_nlp/ocp_nlp_gn_sqp.h"
#include "acados/sim/model_wrapper.h"
#include "acados/sim/sim_erk_integrator.h"
#include "acados/sim/sim_rk_common.h"
#include "acados/utils/types.h"

%}

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

static bool qp_dimensions_equal(const ocp_qp_in *qp1, const ocp_qp_in *qp2) {
    if (qp1->N != qp2->N)
        return false;
    int_t N = qp1->N;
    for (int_t i = 0; i < N; i++) {
        if (qp1->nx[i] != qp2->nx[i])
            return false;
        else if (qp1->nu[i] != qp2->nu[i])
            return false;
        else if (qp1->nb[i] != qp2->nb[i])
            return false;
        else if (qp1->nc[i] != qp2->nc[i])
            return false;
    }
    if (qp1->nx[N] != qp2->nx[N])
        return false;
    else if (qp1->nb[N] != qp2->nb[N])
        return false;
    else if (qp1->nc[N] != qp2->nc[N])
        return false;
    return true;
}
%}

%ignore ocp_qp_out;
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
        ocp_qp_in *qp_in = (ocp_qp_in *) malloc(sizeof(ocp_qp_in));
        if (!is_valid_ocp_dimensions_map(input_map)) {
            char err_msg[256];
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
        allocate_ocp_qp_in(N, nx, nu, nb, nc, qp_in);
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
    ocp_qp_solver(const char *solver_name, ocp_qp_in *qp_in) {
        ocp_qp_solver *solver = (ocp_qp_solver *) malloc(sizeof(ocp_qp_solver));
        void *args = NULL;
        void *mem = NULL;
        int_t workspace_size;
        void *workspace = NULL;
        if (!strcmp(solver_name, "condensing_qpoases")) {
            solver->fun = ocp_qp_condensing_qpoases;
            args = (ocp_qp_condensing_qpoases_args *) \
                malloc(sizeof(ocp_qp_condensing_qpoases_args));
#if defined(OOQP)
        } else if (!strcmp(solver_name, "ooqp")) {
            solver->fun = ocp_qp_ooqp;
            args = (ocp_qp_ooqp_args *) malloc(sizeof(ocp_qp_ooqp_args));
            mem = (ocp_qp_ooqp_memory *) malloc(sizeof(ocp_qp_ooqp_memory));
            ocp_qp_ooqp_create_memory(qp_in, args, mem);
            workspace_size = ocp_qp_ooqp_calculate_workspace_size(qp_in, args);
            workspace = (void *) malloc(workspace_size);
#endif
        } else if (!strcmp(solver_name, "qpdunes")) {
            solver->fun = ocp_qp_qpdunes;
            args = (ocp_qp_qpdunes_args *) malloc(sizeof(ocp_qp_qpdunes_args));
            ocp_qp_qpdunes_create_arguments(args, QPDUNES_DEFAULT_ARGUMENTS);
            mem = (ocp_qp_qpdunes_memory *) malloc(sizeof(ocp_qp_qpdunes_memory));
            ocp_qp_qpdunes_create_memory(qp_in, args, mem);
            workspace_size = ocp_qp_qpdunes_calculate_workspace_size(qp_in, args);
            workspace = (void *) malloc(workspace_size);
        }  else {
            throw std::invalid_argument("Solver name not known!");
            return NULL;
        }
        solver->qp_in = qp_in;
        ocp_qp_out *qp_out = (ocp_qp_out *) malloc(sizeof(ocp_qp_out));
        allocate_ocp_qp_out(qp_in, qp_out);
        solver->qp_out = qp_out;
        solver->args = args;
        solver->mem = mem;
        solver->work = workspace;
        return solver;
    }

    LangObject *solve() {
        int_t return_code = $self->fun($self->qp_in, $self->qp_out, $self->args, \
            $self->mem, $self->work);
        if (return_code != 0) {
            throw std::runtime_error("qp solver failed!");
        }
        return new_sequence_from($self->qp_out->u[0], $self->qp_in->nu[0]);
    }

    LangObject *solve(ocp_qp_in *qp_in) {
        if (!qp_dimensions_equal(qp_in, $self->qp_in)) {
            throw std::invalid_argument("Not allowed to change dimensions of variables "
                "between calls to solver");
        }
        $self->qp_in = qp_in;
        int_t return_code = $self->fun($self->qp_in, $self->qp_out, $self->args, \
            $self->mem, $self->work);
        if (return_code != 0) {
            throw std::runtime_error("qp solver failed!");
        }
        return new_sequence_from($self->qp_out->u[0], $self->qp_in->nu[0]);
    }
}

%typemap(in) real_t ** ls_cost_matrix {
    $1 = ((ocp_nlp_ls_cost *) arg1->cost)->W;
    int_t W_dimensions[arg1->N+1];
    for (int_t i = 0; i < arg1->N+1; i++) {
        W_dimensions[i] = arg1->nx[i] + arg1->nu[i];
    }
    fill_array_from($input, $1, arg1->N+1, W_dimensions, W_dimensions);
}

%typemap(out) real_t ** ls_cost_matrix {
    int_t W_dimensions[arg1->N+1];
    for (int_t i = 0; i < arg1->N+1; i++) {
        W_dimensions[i] = arg1->nx[i] + arg1->nu[i];
    }
    $result = new_sequence_from((const real_t **) $1, arg1->N+1, W_dimensions, W_dimensions);
}

%ignore ocp_nlp_function;
%ignore ocp_nlp_ls_cost;
%ignore ocp_nlp_stage_cost;
%ignore ocp_nlp_args;
%ignore ocp_nlp_memory;
%ignore ocp_nlp_work;
%ignore ocp_nlp_out;
%ignore ocp_nlp_calculate_workspace_size;
%ignore ocp_nlp_cast_workspace;
%ignore ocp_nlp_create_memory;
%ignore ocp_nlp_free_memory;
%rename(ocp_nlp) ocp_nlp_in;
%include "acados/ocp_nlp/ocp_nlp_common.h"

%{
void ocp_nlp_in_ls_cost_matrix_set(ocp_nlp_in *nlp, real_t **matrix) {
    ((ocp_nlp_ls_cost *) nlp->cost)->W = matrix;
}

real_t **ocp_nlp_in_ls_cost_matrix_get(ocp_nlp_in *nlp) {
    return ((ocp_nlp_ls_cost *) nlp->cost)->W;
}
%}

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
    real_t **ls_cost_matrix;
    ocp_nlp_in(LangObject *input_map) {
        ocp_nlp_in *nlp_in = (ocp_nlp_in *) malloc(sizeof(ocp_nlp_in));
        if (!is_valid_ocp_dimensions_map(input_map)) {
            char err_msg[256];
            snprintf(err_msg, sizeof(err_msg), "Input must be a valid OCP %s that specifies at "
                "least N, nx, nu", LANG_MAP_NAME);
            throw std::invalid_argument(err_msg);
        }
        int_t N = int_from(input_map, "N");
        int_t nx[N+1], nu[N+1], nb[N+1], nc[N+1], ng[N+1];
        fill_array_from(input_map, "nx", nx, N+1);
        fill_array_from(input_map, "nu", nu, N+1);
        fill_array_from(input_map, "nb", nb, N+1);
        fill_array_from(input_map, "nc", nc, N+1);
        fill_array_from(input_map, "ng", ng, N+1);
        nu[N] = 0;
        // Default behavior is that initial state is fixed
        if (!has(input_map, "nb")) {
            nb[0] = nx[0];
        }
        allocate_ocp_nlp_in(N, nx, nu, nb, nc, ng, nlp_in);
        if (!has(input_map, "nb")) {
            int idxb[nb[0]];
            for (int_t i = 0; i < nb[0]; i++)
                idxb[i] = i;
            memcpy((void *) nlp_in->idxb[0], idxb, sizeof(idxb));
        }
        return nlp_in;
    }

    void set_model(char *model_name) {
        char library_name[256], path_to_library[256];
        snprintf(library_name, sizeof(library_name), "%s.so", model_name);
        snprintf(path_to_library, sizeof(path_to_library), "./%s", library_name);
        char command[256];
        snprintf(command, sizeof(command), "cc -fPIC -shared -g %s.c -o %s", \
            model_name, library_name);
        int compilation_failed = system(command);
        if (compilation_failed)
            throw std::runtime_error("Something went wrong when compiling the model.");
        void *handle;
        handle = dlopen(path_to_library, RTLD_LAZY);
        if (handle == 0) {
            char err_msg[256];
            snprintf(err_msg, sizeof(err_msg), \
                "Something went wrong when loading the model. dlerror(): %s", dlerror());
            throw std::runtime_error(err_msg);
        }
        typedef int (*eval_t)(const double** arg, double** res, int* iw, double* w, int mem);
        eval_t eval = (eval_t)dlsym(handle, model_name);
        for (int_t i = 0; i < $self->N; i++) {
            $self->sim[i].in->vde = eval;
            $self->sim[i].in->VDE_forw = &vde_fun;
        }
        // dlclose(handle);
    }
}

%extend ocp_nlp_solver {
    ocp_nlp_solver(char *solver_name, ocp_nlp_in *nlp_in) {
        ocp_nlp_solver *solver = (ocp_nlp_solver *) malloc(sizeof(ocp_nlp_solver));
        void *args = NULL;
        void *mem = NULL;
        int_t workspace_size;
        void *workspace = NULL;
        if (!strcmp("gauss-newton-sqp", solver_name)) {
            solver->fun = ocp_nlp_gn_sqp;
            args = (ocp_nlp_gn_sqp_args *) malloc(sizeof(ocp_nlp_gn_sqp_args));
          ((ocp_nlp_gn_sqp_args *) args)->common = (ocp_nlp_args *) malloc(sizeof(ocp_nlp_args));
            snprintf(((ocp_nlp_gn_sqp_args *) args)->qp_solver_name, \
                sizeof(((ocp_nlp_gn_sqp_args *) args)->qp_solver_name), "qpdunes");
            mem = (ocp_nlp_gn_sqp_memory *) malloc(sizeof(ocp_nlp_gn_sqp_memory));
            ((ocp_nlp_gn_sqp_memory *) mem)->common = \
                (ocp_nlp_memory *) malloc(sizeof(ocp_nlp_memory));
            ocp_nlp_gn_sqp_create_memory(nlp_in, args, mem);
          ((ocp_nlp_gn_sqp_args *) args)->common = (ocp_nlp_args *) malloc(sizeof(ocp_nlp_args));
            workspace_size = ocp_nlp_gn_sqp_calculate_workspace_size(nlp_in, args);
            workspace = (void *) malloc(workspace_size);
            int_t N = nlp_in->N;
            ((ocp_nlp_gn_sqp_args *) args)->common->maxIter = 1;
            nlp_in->freezeSens = false;
            for (int_t i = 0; i < N; i++) {
                nlp_in->sim[i].in->nx = nlp_in->nx[i];
                nlp_in->sim[i].in->nu = nlp_in->nu[i];
                nlp_in->sim[i].in->sens_forw = true;
                nlp_in->sim[i].in->sens_adj = false;
                nlp_in->sim[i].in->sens_hess = false;
                nlp_in->sim[i].in->nsens_forw = nlp_in->nx[i] + nlp_in->nu[i];
                nlp_in->sim[i].in->nSteps = 2;
                nlp_in->sim[i].in->step = 0.1;
                nlp_in->sim[i].args = (void*) malloc(sizeof(sim_RK_opts));
                sim_erk_create_arguments(nlp_in->sim[i].args, 4);
                int_t erk_workspace_size = sim_erk_calculate_workspace_size(nlp_in->sim[i].in, \
                    nlp_in->sim[i].args);
                nlp_in->sim[i].work = (void *) malloc(erk_workspace_size);
                nlp_in->sim[i].fun = &sim_erk;
            }
        } else {
            throw std::invalid_argument("Solver name not known!");
            return NULL;
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

    LangObject *solve() {
        // _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & (~_MM_MASK_INVALID));
        int_t return_code = $self->fun($self->nlp_in, $self->nlp_out, $self->args, \
            $self->mem, $self->work);
        if (return_code != 0) {
            throw std::runtime_error("nlp solver failed!");
        }
        return new_sequence_from((const real_t **) $self->nlp_out->x, $self->nlp_in->N+1,
            $self->nlp_in->nx);
    }

    LangObject *solve(LangObject *x0) {
        fill_array_from(x0, (real_t **) $self->nlp_in->lb, 1, $self->nlp_in->nx);
        fill_array_from(x0, (real_t **) $self->nlp_in->ub, 1, $self->nlp_in->nx);
        int_t return_code = $self->fun($self->nlp_in, $self->nlp_out, $self->args, \
            $self->mem, $self->work);
        if (return_code != 0) {
            throw std::runtime_error("nlp solver failed!");
        }
        return new_sequence_from((const real_t **) $self->nlp_out->x, $self->nlp_in->N+1,
            $self->nlp_in->nx);
    }
}
