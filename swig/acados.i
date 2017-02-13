%module acados
%{
#define SWIG_FILE_WITH_INIT
#define PyArray_SimpleNewFromDataF(nd, dims, typenum, data) \
        PyArray_New(&PyArray_Type, nd, dims, typenum, NULL, \
                    data, 0, NPY_ARRAY_FARRAY, NULL)
#include <typeinfo>
#include "acados/utils/types.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_condensing_qpoases.h"
#include "acados/utils/allocate_ocp_qp.h"
%}

%include "numpy.i"
%fragment("NumPy_Fragments");
%init %{
import_array();
%}

%include "acados/utils/types.h"

%{
static bool is_valid_integer(PyObject *input) {
    if (!PyInt_Check(input))
        return false;
    return true;
}

static bool is_sequence_with_length(PyObject *input, int_t expected_length) {
    if (!PySequence_Check(input))
        return false;
    int_t length_of_sequence = PySequence_Length(input);
    if (length_of_sequence != expected_length) {
        char err_msg[256];
        snprintf(err_msg, sizeof(err_msg), "Length of sequence must be %d", expected_length);
        SWIG_Error(SWIG_ValueError, err_msg);
    }
    return true;
}

static bool is_valid_2dim_array(PyObject * const input) {
    if (!PyArray_Check(input))
        return false;
    if (array_numdims(input) != 2)
        return false;
    return true;
}

static bool is_valid_1dim_array(PyObject * const input) {
    if (!PyArray_Check(input))
        return false;
    if (array_numdims(input) != 1)
        return false;
    return true;
}

static bool key_has_valid_integer_value(PyObject *dictionary, const char *key) {
    PyObject *value = PyDict_GetItemString(dictionary, key);
    if (value == NULL) {
        char err_msg[256];
        snprintf(err_msg, sizeof(err_msg), "Input dictionary must have an '%s' key "
            "with as value an integer number", key);
        SWIG_Error(SWIG_ValueError, err_msg);
    } else if (!is_valid_integer(value)) {
        char err_msg[256];
        snprintf(err_msg, sizeof(err_msg), "'%s' must be an integer number", key);
        SWIG_Error(SWIG_ValueError, err_msg);
    }
    return true;
}

static bool key_has_valid_integer_or_sequence_value(PyObject *dictionary,
    const char *key, int_t expected_length) {

    PyObject *value = PyDict_GetItemString(dictionary, key);
    if (value == NULL) {
        char err_msg[256];
        snprintf(err_msg, sizeof(err_msg), "Input dictionary must have an '%s' key "
            "with as value an integer number or a list of integer numbers", key);
        SWIG_Error(SWIG_ValueError, err_msg);
    } else if (!is_valid_integer(value) && !is_sequence_with_length(value, expected_length)) {
        char err_msg[256];
        snprintf(err_msg, sizeof(err_msg), "'%s' must be an integer number or a sequence with "
            "length %d", key, expected_length);
        SWIG_Error(SWIG_ValueError, err_msg);
    }
    return true;
}

static bool is_valid_qp_dictionary(PyObject * const input) {
    if (!PyDict_Check(input))
        return false;
    if (!key_has_valid_integer_value(input, "N"))
        return false;
    int_t N = (int_t) PyInt_AsLong(PyDict_GetItemString(input, "N"));
    if (!key_has_valid_integer_or_sequence_value(input, "nx", N+1))
        return false;
    if (!key_has_valid_integer_or_sequence_value(input, "nu", N))
        return false;

    return true;
}

static void fill_array_from_sequence(int_t * const array, const int_t length_of_array,
    PyObject * const sequence) {

    for (int_t i = 0; i < length_of_array; i++) {
        PyObject *item = PySequence_GetItem(sequence, i);
        if (!is_valid_integer(item)) {
            Py_XDECREF(item);
            SWIG_Error(SWIG_ValueError, "Sequence elements must be integer numbers");
        }
        array[i] = (int_t) PyInt_AsLong(item);
        Py_DECREF(item);
    }
}

static PyObject *convert_to_sequence(int_t *array, const int_t length) {
    PyObject *sequence = PyList_New(length);
    for (int_t i = 0; i < length; i++) {
        PyList_SetItem(sequence, i, PyInt_FromLong((long) array[i]));
    }
    return sequence;
}

static void convert_to_c_array(PyObject *input, int_t * const array, const int_t length_of_array) {
    if (is_valid_integer(input)) {
        int_t integer_number = (int_t) PyInt_AsLong(input);
        for (int_t i = 0; i < length_of_array; i++)
            array[i] = integer_number;
    } else if (is_sequence_with_length(input, length_of_array)) {
        fill_array_from_sequence(array, length_of_array, input);
    } else {
        char err_msg[256];
        snprintf(err_msg, sizeof(err_msg), "Expected integer number or sequence of "
            "length %d", length_of_array);
        SWIG_Error(SWIG_ValueError, err_msg);
    }
}

template<typename T>
int_t get_numpy_type() {
    if (typeid(T) == typeid(real_t))
        return NPY_DOUBLE;
    else if (typeid(T) == typeid(int_t))
        return NPY_INT;
    return NPY_NOTYPE;
}

template<typename T>
static PyObject *convert_to_sequence_of_2dim_arrays(T **c_array,
    const int_t length, const int_t *dimensions1, const int_t *dimensions2) {

    PyObject *pModule = PyImport_Import(PyString_FromString("list_of_arrays"));
    PyObject *pDict = PyModule_GetDict(pModule);
    PyObject *pClass = PyDict_GetItemString(pDict, "list_of_arrays");
    PyObject *sequence = NULL;
    if (PyCallable_Check(pClass)) {
        sequence = PyObject_CallObject(pClass, NULL);
    }
    for (int_t i = 0; i < length; i++) {
        npy_intp dims[2] = {dimensions1[i], dimensions2[i]};
        PyObject *py_array = PyArray_SimpleNewFromDataF(2, dims, get_numpy_type<T>(), \
            (void*) c_array[i]);
        PyList_Append(sequence, py_array);
    }
    return sequence;
}

template<typename T>
static PyObject *convert_to_sequence_of_1dim_arrays(T **c_array,
    const int_t length, const int_t *dimensions) {

    PyObject *pModule = PyImport_Import(PyString_FromString("list_of_arrays"));
    PyObject *pDict = PyModule_GetDict(pModule);
    PyObject *pClass = PyDict_GetItemString(pDict, "list_of_arrays");
    PyObject *sequence = NULL;
    if (PyCallable_Check(pClass)) {
        sequence = PyObject_CallObject(pClass, NULL);
    }
    for (int_t i = 0; i < length; i++) {
        npy_intp dims[1] = {dimensions[i]};
        PyObject *py_array = PyArray_SimpleNewFromDataF(1, dims, get_numpy_type<T>(), \
            (void*) c_array[i]);
        PyList_Append(sequence, py_array);
    }
    return sequence;
}

template<typename T>
static PyArrayObject *array_with_type(PyObject *input) {
    PyArrayObject *f_array = (PyArrayObject *) PyArray_FROM_OF(input, NPY_ARRAY_F_CONTIGUOUS);
    PyObject *obj = PyArray_Cast(f_array, get_numpy_type<T>());
    return reinterpret_cast<PyArrayObject *>(obj);
}

template<typename T>
static PyArrayObject *object_to_2dim_array(PyObject *input, int_t dim1, int_t dim2) {
    PyArrayObject *input_array = array_with_type<T>(input);
    if (PyArray_NDIM(input_array) != 2) {
        SWIG_Error(SWIG_ValueError, "Expected a 2D array as input");
    }
    npy_intp *dims = PyArray_DIMS(input_array);
    if (dims[0] != dim1 || dims[1] != dim2) {
        SWIG_Error(SWIG_ValueError, "Input array with wrong dimensions");
    }
    return input_array;
}

template<typename T>
static PyArrayObject *object_to_1dim_array(PyObject *input, int_t dim) {
    PyArrayObject *input_array = array_with_type<T>(input);
    if (PyArray_NDIM(input_array) != 1) {
        SWIG_Error(SWIG_ValueError, "Expected a 1D array as input");
    }
    npy_intp *dims = PyArray_DIMS(input_array);
    if (dims[0] != dim) {
        SWIG_Error(SWIG_ValueError, "Input array with wrong dimensions");
    }
    return input_array;
}

template<typename T>
static void convert_to_2dim_c_array(PyObject * const input, T ** const array,
    const int_t length_of_array, const int_t *dimensions1, const int_t *dimensions2) {

    if (is_valid_2dim_array(input)) {
        int_t dim1 = dimensions1[0];
        int_t dim2 = dimensions2[0];
        PyArrayObject *input_array = object_to_2dim_array<T>(input, dim1, dim2);
        for (int_t i = 0; i < length_of_array; i++) {
            memcpy((void *) array[i], (T *) array_data(input_array), dim1*dim2*sizeof(T));
        }
    } else if (is_sequence_with_length(input, length_of_array)) {
        for (int_t i = 0; i < length_of_array; i++) {
            PyObject *item = PySequence_GetItem(input, i);
            if (is_valid_2dim_array(item)) {
                int_t dim1 = dimensions1[i];
                int_t dim2 = dimensions2[i];
                PyArrayObject *input_array = object_to_2dim_array<T>(item, dim1, dim2);
                memcpy((void *) array[i], (T *) array_data(input_array), dim1*dim2*sizeof(T));
            }
        }
    } else {
        SWIG_Error(SWIG_ValueError, "Expected a 2-dimensional array as input");
    }
}

template<typename T>
static void convert_to_1dim_c_array(PyObject * const input, T ** const array,
    const int_t length_of_array, const int_t *dimensions) {

    if (is_valid_1dim_array(input)) {
        PyArrayObject *input_array = object_to_1dim_array<T>(input, dimensions[0]);
        for (int_t i = 0; i < length_of_array; i++) {
            memcpy((void *) array[i], (T *) array_data(input_array), dimensions[0]*sizeof(T));
        }
    } else if (is_sequence_with_length(input, length_of_array)) {
        for (int_t i = 0; i < length_of_array; i++) {
            PyObject *item = PySequence_GetItem(input, i);
            if (is_valid_1dim_array(item)) {
                PyArrayObject *input_array = object_to_1dim_array<T>(item, dimensions[i]);
                memcpy((void *) array[i], (T *) array_data(input_array), dimensions[i]*sizeof(T));
            }
        }
    } else {
        SWIG_Error(SWIG_ValueError, "Expected a 1-dimensional array as input");
    }
}

template<typename T>
static void read_array_from_dictionary(PyObject *dictionary, const char *key_name,
    T *array, int_t array_length) {

    if (PyDict_GetItemString(dictionary, key_name) == NULL) {
        memset(array, 0, array_length*sizeof(*array));
    } else {
        PyObject *sequence = PyDict_GetItemString(dictionary, key_name);
        convert_to_c_array(sequence, array, array_length);
    }
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

%typemap(in) int_t N {
    $1 = ($1_ltype) arg1->$1_name;
    SWIG_Error(SWIG_ValueError, "It's not allowed to change number of stages");
}

%typemap(in) const int_t * nx {
    SWIG_Error(SWIG_ValueError, "It's not allowed to change dimension of state vector");
}

%typemap(out) const int_t * nx {
    return convert_to_sequence($1, arg1->N+1);
}

%typemap(in) const int_t * nu {
    SWIG_Error(SWIG_ValueError, "It's not allowed to change dimension of vector of controls");
}

%typemap(out) const int_t * nu {
    return convert_to_sequence($1, arg1->N);
}

%typemap(in) const int_t * nb {
    SWIG_Error(SWIG_ValueError, "It's not allowed to change number of bounds");
}

%typemap(out) const int_t * nb {
    return convert_to_sequence($1, arg1->N+1);
}

%typemap(in) const int_t * nc {
    SWIG_Error(SWIG_ValueError, "It's not allowed to change number of polytopic constraints");
}

%typemap(out) const int_t * nc {
    return convert_to_sequence($1, arg1->N+1);
}

%typemap(in) const real_t ** A {
    $1 = ($1_ltype) arg1->$1_name;
    convert_to_2dim_c_array($input, $1, arg1->N, &(arg1->nx[1]), &(arg1->nx[0]));
}

%typemap(out) const real_t ** A {
    $result = convert_to_sequence_of_2dim_arrays($1, arg1->N, &arg1->nx[1], &arg1->nx[0]);
}

%typemap(in) const real_t ** B {
    $1 = ($1_ltype) arg1->$1_name;
    convert_to_2dim_c_array($input, $1, arg1->N, &(arg1->nx[1]), &(arg1->nu[0]));
}

%typemap(out) const real_t ** B {
    $result = convert_to_sequence_of_2dim_arrays($1, arg1->N, &arg1->nx[1], &arg1->nu[0]);
}

%typemap(in) const real_t ** b {
    $1 = ($1_ltype) arg1->$1_name;
    convert_to_1dim_c_array($input, $1, arg1->N, &(arg1->nx[1]));
}

%typemap(out) const real_t ** b {
    $result = convert_to_sequence_of_1dim_arrays($1, arg1->N, &(arg1->nx[1]));
}

%typemap(in) const real_t ** Q {
    $1 = ($1_ltype) arg1->$1_name;
    convert_to_2dim_c_array($input, $1, arg1->N+1, arg1->nx, arg1->nx);
}

%typemap(out) const real_t ** Q {
    $result = convert_to_sequence_of_2dim_arrays($1, arg1->N+1, arg1->nx, arg1->nx);
}

%typemap(in) const real_t ** S {
    $1 = ($1_ltype) arg1->$1_name;
    convert_to_2dim_c_array($input, $1, arg1->N, arg1->nu, arg1->nx);
}

%typemap(out) const real_t ** S {
    $result = convert_to_sequence_of_2dim_arrays($1, arg1->N, arg1->nu, arg1->nx);
}

%typemap(in) const real_t ** R {
    $1 = ($1_ltype) arg1->$1_name;
    convert_to_2dim_c_array($input, $1, arg1->N, arg1->nu, arg1->nu);
}

%typemap(out) const real_t ** R {
    $result = convert_to_sequence_of_2dim_arrays($1, arg1->N, arg1->nu, arg1->nu);
}

%typemap(in) const real_t ** q {
    $1 = ($1_ltype) arg1->$1_name;
    convert_to_1dim_c_array($input, $1, arg1->N+1, arg1->nx);
}

%typemap(out) const real_t ** q {
    $result = convert_to_sequence_of_1dim_arrays($1, arg1->N+1, arg1->nx);
}

%typemap(in) const real_t ** r {
    $1 = ($1_ltype) arg1->$1_name;
    convert_to_1dim_c_array($input, $1, arg1->N, arg1->nu);
}

%typemap(out) const real_t ** r {
    $result = convert_to_sequence_of_1dim_arrays($1, arg1->N, arg1->nu);
}

%typemap(in) const int_t ** idxb {
    $1 = ($1_ltype) arg1->$1_name;
    convert_to_1dim_c_array($input, $1, arg1->N+1, arg1->nb);
}

%typemap(out) const int_t ** idxb {
    $result = convert_to_sequence_of_1dim_arrays($1, arg1->N+1, arg1->nb);
}

%typemap(in) const real_t ** lb {
    $1 = ($1_ltype) arg1->$1_name;
    convert_to_1dim_c_array($input, $1, arg1->N+1, arg1->nb);
}

%typemap(out) const real_t ** lb {
    $result = convert_to_sequence_of_1dim_arrays($1, arg1->N+1, arg1->nb);
}

%typemap(in) const real_t ** ub {
    $1 = ($1_ltype) arg1->$1_name;
    convert_to_1dim_c_array($input, $1, arg1->N+1, arg1->nb);
}

%typemap(out) const real_t ** ub {
    $result = convert_to_sequence_of_1dim_arrays($1, arg1->N+1, arg1->nb);
}

%typemap(in) const real_t ** Cx {
    $1 = ($1_ltype) arg1->$1_name;
    convert_to_2dim_c_array($input, $1, arg1->N+1, arg1->nc, arg1->nx);
}

%typemap(out) const real_t ** Cx {
    $result = convert_to_sequence_of_2dim_arrays($1, arg1->N+1, arg1->nc, arg1->nx);
}

%typemap(in) const real_t ** Cu {
    $1 = ($1_ltype) arg1->$1_name;
    convert_to_2dim_c_array($input, $1, arg1->N, arg1->nc, arg1->nu);
}

%typemap(out) const real_t ** Cu {
    $result = convert_to_sequence_of_2dim_arrays($1, arg1->N, arg1->nc, arg1->nu);
}

%typemap(in) const real_t ** lc {
    $1 = ($1_ltype) arg1->$1_name;
    convert_to_1dim_c_array($input, $1, arg1->N+1, arg1->nc);
}

%typemap(out) const real_t ** lc {
    $result = convert_to_sequence_of_1dim_arrays($1, arg1->N+1, arg1->nc);
}

%typemap(in) const real_t ** uc {
    $1 = ($1_ltype) arg1->$1_name;
    convert_to_1dim_c_array($input, $1, arg1->N+1, arg1->nc);
}

%typemap(out) const real_t ** uc {
    $result = convert_to_sequence_of_1dim_arrays($1, arg1->N+1, arg1->nc);
}

%include "acados/ocp_qp/ocp_qp_common.h"

%extend ocp_qp_in {
    ocp_qp_in(PyObject *dictionary) {
        ocp_qp_in *qp_in = (ocp_qp_in *) malloc(sizeof(ocp_qp_in));
        if (!is_valid_qp_dictionary(dictionary)) {
            SWIG_Error(SWIG_ValueError, "Input must be a dictionary");
        }
        int_t N = (int_t) PyInt_AsLong(PyDict_GetItemString(dictionary, "N"));
        int_t nx[N+1], nu[N], nb[N+1], nc[N+1];
        read_array_from_dictionary(dictionary, "nx", nx, N+1);
        read_array_from_dictionary(dictionary, "nu", nu, N);
        read_array_from_dictionary(dictionary, "nb", nb, N+1);
        read_array_from_dictionary(dictionary, "nc", nc, N+1);
        allocate_ocp_qp_in(N, nx, nu, nb, nc, qp_in);
        return qp_in;
    }
}

%extend ocp_qp_solver {
    ocp_qp_solver(const char *solver_name, ocp_qp_in *qp_in) {
        ocp_qp_solver *solver = (ocp_qp_solver *) malloc(sizeof(ocp_qp_solver));
        if (strcmp(solver_name, "condensing_qpoases")) {
            SWIG_Error(SWIG_ValueError, "Solver name not known!");
            return NULL;
        }
        solver->fun = ocp_qp_condensing_qpoases;
        solver->qp_in = qp_in;
        ocp_qp_out *qp_out = (ocp_qp_out *) malloc(sizeof(ocp_qp_out));
        allocate_ocp_qp_out(qp_in, qp_out);
        solver->qp_out = qp_out;
        ocp_qp_condensing_qpoases_args *args = \
            (ocp_qp_condensing_qpoases_args *) malloc(sizeof(ocp_qp_condensing_qpoases_args));
        solver->mem = args;
        real_t *qpoases_work = NULL;
        solver->work = qpoases_work;
        initialise_qpoases(qp_in);
        return solver;
    }
    PyObject *solve() {
        int_t return_code = $self->fun($self->qp_in, $self->qp_out, $self->mem, $self->work);
        if (return_code != 0) {
            SWIG_Error(SWIG_RuntimeError, "qp solver failed!");
        }

        return convert_to_sequence_of_1dim_arrays($self->qp_out->x, \
            $self->qp_in->N+1, $self->qp_in->nx);
    }
    PyObject *solve(ocp_qp_in *qp_in) {
        if(!qp_dimensions_equal(qp_in, $self->qp_in)) {
            SWIG_Error(SWIG_ValueError, "Not allowed to change dimensions of variables "
                "between calls to solver");
        }
        $self->qp_in = qp_in;
        int_t return_code = $self->fun($self->qp_in, $self->qp_out, $self->mem, $self->work);
        if (return_code != 0) {
            SWIG_Error(SWIG_RuntimeError, "qp solver failed!");
        }
        return convert_to_sequence_of_1dim_arrays($self->qp_out->x, \
            $self->qp_in->N+1, $self->qp_in->nx);
    }
}

%include "acados/ocp_qp/ocp_qp_condensing_qpoases.h"
