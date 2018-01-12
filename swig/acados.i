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

%include "std_vector.i"
namespace std {
    %template(vectori) vector<int>;
};

%rename($ignore, %$isclass) ""; // Only ignore all classes

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

#include <iostream>

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados_c/ocp_qp.h"
#include "acados_cpp/ocp_qp.h"

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

bool qp_dimensions_equal(const ocp_qp_dims *qp1, const ocp_qp_dims *qp2) {
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
        else if (qp1->ng[i] != qp2->ng[i])
            return false;
    }
    if (qp1->nx[N] != qp2->nx[N])
        return false;
    else if (qp1->nb[N] != qp2->nb[N])
        return false;
    else if (qp1->ng[N] != qp2->ng[N])
        return false;
    return true;
}

ocp_qp_dims *map_to_ocp_qp_dims(const LangObject *map) {
    if (!is_valid_ocp_dimensions_map(map)) {
        std::string err_msg =
            std::string("Input must be a valid OCP %s that specifies at least N, nx, nu")
            + std::string(LANG_MAP_NAME);
        throw std::invalid_argument(err_msg);
    }

    int_t N = int_from(map, "N");
    ocp_qp_dims *qp_dims = create_ocp_qp_dims(N);
    fill_array_from(map, "nx", qp_dims->nx, N+1);
    fill_array_from(map, "nu", qp_dims->nu, N+1);
    fill_array_from(map, "nb", qp_dims->nb, N+1);
    fill_array_from(map, "nc", qp_dims->ng, N+1);
    qp_dims->nu[N] = 0;
    // Default behavior is that initial state is fixed
    if (!has(map, "nb")) {
        qp_dims->nb[0] = qp_dims->nx[0];
    }
    return qp_dims;
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

%rename("%s") OcpQp;
%include "acados_cpp/ocp_qp.h"

%extend OcpQp {

    void setQ(LangObject *input) {
        for (int i = 0; i <= $self->N; i++)
            $self->setQ(i, as_dpointer(input));
    }

    void setS(LangObject *input) {
        for (int i = 0; i <= $self->N; i++)
            $self->setS(i, as_dpointer(input));
    }

    void setR(LangObject *input) {
        for (int i = 0; i <= $self->N; i++)
            $self->setR(i, as_dpointer(input));
    }

    void setq(LangObject *input) {
        for (int i = 0; i <= $self->N; i++)
            $self->setq(i, as_dpointer(input));
    }

    void setr(LangObject *input) {
        for (int i = 0; i <= $self->N; i++)
            $self->setr(i, as_dpointer(input));
    }

    void setA(LangObject *input) {
        for (int i = 0; i < $self->N; i++)
            $self->setA(i, as_dpointer(input));
    }

    void setB(LangObject *input) {
        for (int i = 0; i < $self->N; i++)
            $self->setB(i, as_dpointer(input));
    }

    void setb(LangObject *input) {
        for (int i = 0; i < $self->N; i++)
            $self->setb(i, as_dpointer(input));
    }

    void setlb(LangObject *input) {
        for (int i = 0; i <= $self->N; i++)
            $self->setlb(i, as_dpointer(input));
    }

    void setub(LangObject *input) {
        for (int i = 0; i <= $self->N; i++)
            $self->setub(i, as_dpointer(input));
    }

    void setC(LangObject *input) {
        for (int i = 0; i <= $self->N; i++)
            $self->setC(i, as_dpointer(input));
    }

    void setD(LangObject *input) {
        for (int i = 0; i <= $self->N; i++)
            $self->setD(i, as_dpointer(input));
    }

    void setlg(LangObject *input) {
        for (int i = 0; i <= $self->N; i++)
            $self->setlg(i, as_dpointer(input));
    }

    void setug(LangObject *input) {
        for (int i = 0; i <= $self->N; i++)
            $self->setug(i, as_dpointer(input));
    }

    char *__str__() {
        static char tmp[1];
        std::cout << *($self);
        return tmp;
    }
}

%rename("%s") ocp_qp_solver;
%include "acados_c/ocp_qp.h"

%extend ocp_qp_solver {

    ocp_qp_solver(ocp_qp_solver_t solver_name, LangObject *dimensions, LangObject *options = NONE) {

        ocp_qp_dims *dims = map_to_ocp_qp_dims(dimensions);

        ocp_qp_solver_plan plan;
        plan.qp_solver = solver_name;
        
        void *args = ocp_qp_create_args(&plan, dims);
        ocp_qp_solver *solver = ocp_qp_create(&plan, dims, args);
        return solver;
    }

    LangObject *evaluate(ocp_qp_in *input) {
        ocp_qp_out *result = create_ocp_qp_out(input->dim);
        int_t return_code = ocp_qp_solve($self, input, result);
        if (return_code != 0)
            throw std::runtime_error("qp solver failed!");
        return ocp_qp_output(input, result);
    }

}
