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
%}

/* ------------------------
/        OCP 
/  ----------------------*/
%{

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados_c/ocp_qp.h"

bool is_valid_ocp_dimensions_map(const LangObject *input) {
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

LangObject *ocp_qp_output(const ocp_qp_in *in, const ocp_qp_out *out) {
    real_t **states_copy, **controls_copy;
    ocp_qp_dims *dims = in->dim;
    states_copy = (real_t **) malloc((dims->N+1) * sizeof(real_t *));
    controls_copy = (real_t **) malloc((dims->N+1) * sizeof(real_t *));
    for (int_t i = 0; i <= dims->N; i++) {
        states_copy[i] = (real_t *) calloc(dims->nx[i], sizeof(real_t));
        for (int_t j = 0; j < dims->nx[i]; j++)
            states_copy[i][j] = DVECEL_LIBSTR(out->ux, dims->nu[i] + j);
        controls_copy[i] = (real_t *) calloc(dims->nu[i], sizeof(real_t));
        for (int_t j = 0; j < dims->nu[i]; j++)
            controls_copy[i][j] = DVECEL_LIBSTR(out->ux, j);
    }

    LangObject *x_star = new_sequence_from(states_copy, dims->N+1, dims->nx);
    LangObject *u_star = new_sequence_from(controls_copy, dims->N+1, dims->nu);
    return new_ocp_output_tuple(x_star, u_star);
}

%}

%rename($ignore, %$isclass) ""; // Only ignore all classes
%rename("%s") ocp_qp_solver;
%include "acados_c/ocp_qp.h"

%{
void ocp_qp_solver_blabla_set(ocp_qp_solver *solver, ocp_qp_in *blabla) {
    solver->blabla = blabla;
}

void ocp_qp_solver_result_set(ocp_qp_solver *solver, ocp_qp_out *result) {
    solver->result = result;
}

ocp_qp_in *ocp_qp_solver_blabla_get(ocp_qp_solver *solver) {
    return solver->blabla;
}

ocp_qp_out *ocp_qp_solver_result_get(ocp_qp_solver *solver) {
    return solver->result;
}
%}

%extend ocp_qp_solver {

    ocp_qp_in *blabla = NULL;
    ocp_qp_out *result = NULL;

    ocp_qp_solver(ocp_qp_solver_t solver_name, ocp_qp_in *qp_in, LangObject *options = NONE) {
        
        ocp_qp_solver *solver = (ocp_qp_solver *) malloc(sizeof(ocp_qp_solver *));
        
        solver->blabla = qp_in;
        solver->result = create_ocp_qp_out(blabla->dim);
        
        ocp_qp_dims *dims = qp_in->dim;
        
        ocp_qp_solver_plan plan;
        plan.qp_solver = solver_name;
        
        void *args = ocp_qp_create_args(&plan, dims);
        // ocp_qp_solver *solver = ocp_qp_create(&plan, dims, args);
        return solver;
    }

    LangObject *evaluate() {
        int_t return_code = ocp_qp_solve($self, $self->blabla, $self->result);
        if (return_code != 0) {
            throw std::runtime_error("qp solver failed!");
        }
        return ocp_qp_output($self->blabla, $self->result);
    }
}
