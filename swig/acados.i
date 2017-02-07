%module acados
%{
#define SWIG_FILE_WITH_INIT
#define PyArray_SimpleNewFromDataF(nd, dims, typenum, data) \
        PyArray_New(&PyArray_Type, nd, dims, typenum, NULL, \
                    data, 0, NPY_ARRAY_FARRAY, NULL)
#include <stdexcept>
#include "acados/utils/types.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/allocate_ocp_qp.h"
%}

%include "std_string.i"
%include "exception.i"

%exception {
    try {
        $action
    } catch (std::invalid_argument& e) {
        SWIG_exception(SWIG_ValueError, const_cast<char*>(e.what()));
    } catch (...) {
        SWIG_exception(SWIG_UnknownError, "Unknown exception");
    }
}

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
        throw std::invalid_argument(err_msg);
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
        throw std::invalid_argument(err_msg);
    } else if (!is_valid_integer(value)) {
        char err_msg[256];
        snprintf(err_msg, sizeof(err_msg), "'%s' must be an integer number", key);
        throw std::invalid_argument(err_msg);
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
        throw std::invalid_argument(err_msg);
    } else if (!is_valid_integer(value) && !is_sequence_with_length(value, expected_length)) {
        char err_msg[256];
        snprintf(err_msg, sizeof(err_msg), "'%s' must be an integer number or a sequence with "
            "length %d", key, expected_length);
        throw std::invalid_argument(err_msg);
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
            throw std::invalid_argument("Sequence elements must be integer numbers");
        }
        array[i] = (int_t) PyInt_AsLong(item);
        Py_DECREF(item);
    }
}

static PyObject *convert_to_sequence(const int_t *array, const int_t length) {
    PyObject *sequence = PyList_New(length);
    for (int_t i = 0; i < length; i++) {
        PyList_SetItem(sequence, i, PyInt_FromLong((long) array[i]));
    }
    return sequence;
}


static PyObject *convert_to_sequence_of_2dim_arrays(real_t **c_array,
    const int_t length, const int_t *dimensions1, const int_t *dimensions2) {

    PyObject *sequence = PyList_New(length);
    for (int_t i = 0; i < length; i++) {
        npy_intp dims[2] = {dimensions1[i], dimensions2[i]};
        PyObject *py_array = PyArray_SimpleNewFromDataF(2, dims, NPY_DOUBLE, (void*) c_array[i]);
        PyList_SetItem(sequence, i, py_array);
    }
    return sequence;
}

static PyObject *convert_to_sequence_of_1dim_arrays(real_t **c_array,
    const int_t length, const int_t *dimensions) {

    PyObject *sequence = PyList_New(length);
    for (int_t i = 0; i < length; i++) {
        npy_intp dims[1] = {dimensions[i]};
        PyArrayObject *py_array = (PyArrayObject *) PyArray_SimpleNewFromDataF(1, \
            dims, NPY_DOUBLE, (void*) c_array[i]);
        PyList_SetItem(sequence, i, (PyObject *) py_array);
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
        throw std::invalid_argument(err_msg);
    }
}

static PyArrayObject *object_to_2dim_array(PyObject *input, int_t dim1, int_t dim2) {
    PyObject *obj = PyArray_FROM_OTF(input, NPY_DOUBLE, NPY_ARRAY_F_CONTIGUOUS);
    PyArrayObject *input_array = reinterpret_cast<PyArrayObject *>(obj);
    if (PyArray_NDIM(input_array) != 2) {
        throw std::invalid_argument("Expected a 2D array as input");
    }
    npy_intp *dims = PyArray_DIMS(input_array);
    if (dims[0] != dim1 || dims[1] != dim2) {
        throw std::invalid_argument("Input array with wrong dimensions");
    }
    return input_array;
}

static PyArrayObject *object_to_1dim_array(PyObject *input, int_t dim) {
    PyObject *obj = PyArray_FROM_OTF(input, NPY_DOUBLE, NPY_ARRAY_F_CONTIGUOUS);
    PyArrayObject *input_array = reinterpret_cast<PyArrayObject *>(obj);
    if (PyArray_NDIM(input_array) != 1) {
        throw std::invalid_argument("Expected a 1D array as input");
    }
    npy_intp *dims = PyArray_DIMS(input_array);
    if (dims[0] != dim) {
        throw std::invalid_argument("Input array with wrong dimensions");
    }
    return input_array;
}

static void convert_to_2dim_c_array(PyObject * const input, real_t ** const array,
    const int_t length_of_array, const int_t *dimensions1, const int_t *dimensions2) {

    if (is_valid_2dim_array(input)) {
        int_t dim1 = dimensions1[0];
        int_t dim2 = dimensions2[0];
        PyArrayObject *input_array = object_to_2dim_array(input, dim1, dim2);
        for (int_t i = 0; i < length_of_array; i++) {
            memcpy((void *) array[i], (real_t *) array_data(input_array), \
                dim1*dim2*sizeof(real_t));
        }
    } else if (is_sequence_with_length(input, length_of_array)) {
        for (int_t i = 0; i < length_of_array; i++) {
            PyObject *item = PySequence_GetItem(input, i);
            if (is_valid_2dim_array(item)) {
                int_t dim1 = dimensions1[i];
                int_t dim2 = dimensions2[i];
                PyArrayObject *input_array = object_to_2dim_array(item, dim1, dim2);
                memcpy((void *) array[i], (real_t *) array_data(input_array), \
                    dim1*dim2*sizeof(real_t));
            }
        }
    } else {
        throw std::invalid_argument("Expected a 2-dimensional array as input");
    }
}

static void convert_to_1dim_c_array(PyObject * const input, real_t ** const array,
    const int_t length_of_array, const int_t *dimensions) {

    if (is_valid_1dim_array(input)) {
        PyArrayObject *input_array = object_to_1dim_array(input, dimensions[0]);
        for (int_t i = 0; i < length_of_array; i++) {
            memcpy((void *) array[i], (real_t *) array_data(input_array), \
                dimensions[0]*sizeof(real_t));
        }
    } else if (is_sequence_with_length(input, length_of_array)) {
        for (int_t i = 0; i < length_of_array; i++) {
            PyObject *item = PySequence_GetItem(input, i);
            if (is_valid_1dim_array(item)) {
                PyArrayObject *input_array = object_to_1dim_array(item, dimensions[i]);
                memcpy((void *) array[i], (real_t *) array_data(input_array), \
                    dimensions[i]*sizeof(real_t));
            }
        }
    } else {
        throw std::invalid_argument("Expected a 1-dimensional array as input");
    }
}

static void read_int_array_from_dictionary(PyObject *dictionary, const char *key_name,
    int_t *array, int_t array_length) {

    if (PyDict_GetItemString(dictionary, key_name) == NULL) {
        memset(array, 0, array_length*sizeof(*array));
    } else {
        PyObject *sequence = PyDict_GetItemString(dictionary, key_name);
        convert_to_c_array(sequence, array, array_length);
    }
}

static bool strings_equal(const char *string1, const char *string2) {
    if (strcmp(string1, string2) == 0)
        return true;
    else
        return false;
}

%}

%typemap(in) const real_t ** {
    std::string name("$symname");
    if (name.find("ocp_qp_in") == std::string::npos)
        SWIG_exception(SWIG_ValueError, "This SWIG typemap is meant for ocp_qp_in, not $symname");
    $1 = ($1_ltype) arg1->$1_name;
    try {
        if (strings_equal("$1_name", "A")) {
            convert_to_2dim_c_array($input, $1, arg1->N, &(arg1->nx[1]), &(arg1->nx[0]));
        } else if (strings_equal("$1_name", "B")) {
            convert_to_2dim_c_array($input, $1, arg1->N, &(arg1->nx[1]), &(arg1->nu[0]));
        } else if (strings_equal("$1_name", "b")) {
            convert_to_1dim_c_array($input, $1, arg1->N, &(arg1->nx[1]));
        } else if (strings_equal("$1_name", "Q")) {
            convert_to_2dim_c_array($input, $1, arg1->N+1, arg1->nx, arg1->nx);
        } else if (strings_equal("$1_name", "S")) {
            convert_to_2dim_c_array($input, $1, arg1->N, arg1->nu, arg1->nx);
        } else if (strings_equal("$1_name", "R")) {
            convert_to_2dim_c_array($input, $1, arg1->N, arg1->nu, arg1->nu);
        } else if (strings_equal("$1_name", "q")) {
            convert_to_1dim_c_array($input, $1, arg1->N+1, arg1->nx);
        } else if (strings_equal("$1_name", "r")) {
            convert_to_1dim_c_array($input, $1, arg1->N, arg1->nu);
        } else if (strings_equal("$1_name", "lb")) {
            convert_to_1dim_c_array($input, $1, arg1->N+1, arg1->nb);
        } else if (strings_equal("$1_name", "ub")) {
            convert_to_1dim_c_array($input, $1, arg1->N+1, arg1->nb);
        } else if (strings_equal("$1_name", "Cx")) {
            convert_to_2dim_c_array($input, $1, arg1->N+1, arg1->nc, arg1->nx);
        } else if (strings_equal("$1_name", "Cu")) {
            convert_to_2dim_c_array($input, $1, arg1->N, arg1->nc, arg1->nu);
        } else if (strings_equal("$1_name", "lc")) {
            convert_to_1dim_c_array($input, $1, arg1->N+1, arg1->nc);
        } else if (strings_equal("$1_name", "uc")) {
            convert_to_1dim_c_array($input, $1, arg1->N+1, arg1->nc);
        } else {
            throw std::invalid_argument("There is no field called $1_name");
        }
    } catch (std::invalid_argument &e) {
        SWIG_exception(SWIG_ValueError, e.what());
    } catch (...) {
        SWIG_exception(SWIG_UnknownError, "Unknown exception");
    }
}

%typemap(out) const real_t ** {
    std::string name("$symname");
    if (name.find("ocp_qp_in") == std::string::npos)
        SWIG_exception(SWIG_ValueError, "This SWIG typemap is meant for ocp_qp_in, not $symname");
    try {
        if (strings_equal("$1_name", "A")) {
            $result = convert_to_sequence_of_2dim_arrays($1, arg1->N, &arg1->nx[1], &arg1->nx[0]);
        } else if (strings_equal("$1_name", "B")) {
            $result = convert_to_sequence_of_2dim_arrays($1, arg1->N, &arg1->nx[1], &arg1->nu[0]);
        } else if (strings_equal("$1_name", "b")) {
            $result = convert_to_sequence_of_1dim_arrays($1, arg1->N, &(arg1->nx[1]));
        } else if (strings_equal("$1_name", "Q")) {
            $result = convert_to_sequence_of_2dim_arrays($1, arg1->N+1, arg1->nx, arg1->nx);
        } else if (strings_equal("$1_name", "S")) {
            $result = convert_to_sequence_of_2dim_arrays($1, arg1->N, arg1->nu, arg1->nx);
        } else if (strings_equal("$1_name", "R")) {
            $result = convert_to_sequence_of_2dim_arrays($1, arg1->N, arg1->nu, arg1->nu);
        } else if (strings_equal("$1_name", "q")) {
            $result = convert_to_sequence_of_1dim_arrays($1, arg1->N+1, arg1->nx);
        } else if (strings_equal("$1_name", "r")) {
            $result = convert_to_sequence_of_1dim_arrays($1, arg1->N, arg1->nu);
        } else if (strings_equal("$1_name", "lb")) {
            $result = convert_to_sequence_of_1dim_arrays($1, arg1->N+1, arg1->nb);
        } else if (strings_equal("$1_name", "ub")) {
            $result = convert_to_sequence_of_1dim_arrays($1, arg1->N+1, arg1->nb);
        } else if (strings_equal("$1_name", "Cx")) {
            $result = convert_to_sequence_of_2dim_arrays($1, arg1->N+1, arg1->nc, arg1->nx);
        } else if (strings_equal("$1_name", "Cu")) {
            $result = convert_to_sequence_of_2dim_arrays($1, arg1->N, arg1->nc, arg1->nu);
        } else if (strings_equal("$1_name", "lc")) {
            $result = convert_to_sequence_of_1dim_arrays($1, arg1->N+1, arg1->nc);
        } else if (strings_equal("$1_name", "uc")) {
            $result = convert_to_sequence_of_1dim_arrays($1, arg1->N+1, arg1->nc);
        } else {
            throw std::invalid_argument("There is no field called $1_name");
        }
    } catch (std::invalid_argument &e) {
        SWIG_exception(SWIG_ValueError, e.what());
    } catch (...) {
        SWIG_exception(SWIG_UnknownError, "Unknown exception");
    }
}

%include "acados/ocp_qp/ocp_qp_common.h"

%extend ocp_qp_in {
    ocp_qp_in(PyObject *dictionary) {
        if (!is_valid_qp_dictionary(dictionary)) {
            throw std::invalid_argument("Input must be a dictionary");
        }
        ocp_qp_in *qp = (ocp_qp_in *) malloc(sizeof(ocp_qp_in));
        int_t N = (int_t) PyInt_AsLong(PyDict_GetItemString(dictionary, "N"));
        int_t nx[N+1], nu[N], nb[N+1], nc[N+1];
        read_int_array_from_dictionary(dictionary, "nx", nx, N+1);
        read_int_array_from_dictionary(dictionary, "nu", nu, N);
        read_int_array_from_dictionary(dictionary, "nb", nb, N+1);
        read_int_array_from_dictionary(dictionary, "nc", nc, N+1);
        allocate_ocp_qp_in(N, nx, nu, nb, nc, qp);
        return qp;
    }
    PyObject *get_nx() {
        return convert_to_sequence($self->nx, $self->N+1);
    }
    void set_nx(PyObject *input) {
        convert_to_c_array(input, (int_t *) $self->nx, $self->N+1);
    }
    PyObject *get_nu() {
        return convert_to_sequence($self->nu, $self->N);
    }
    void set_nu(PyObject *input) {
        convert_to_c_array(input, (int_t *) $self->nu, $self->N);
    }
    PyObject *get_nb() {
        return convert_to_sequence($self->nb, $self->N+1);
    }
    void set_nb(PyObject *input) {
        convert_to_c_array(input, (int_t *) $self->nb, $self->N+1);
    }
    PyObject *get_nc() {
        return convert_to_sequence($self->nc, $self->N+1);
    }
    void set_nc(PyObject *input) {
        convert_to_c_array(input, (int_t *) $self->nc, $self->N+1);
    }
    %pythoncode %{
    __swig_getmethods__["nx"] = get_nx
    __swig_setmethods__["nx"] = set_nx
    __swig_getmethods__["nu"] = get_nu
    __swig_setmethods__["nu"] = set_nu
    __swig_getmethods__["nb"] = get_nb
    __swig_setmethods__["nb"] = set_nb
    __swig_getmethods__["nc"] = get_nc
    __swig_setmethods__["nc"] = set_nc
    if _newclass: nx = property(get_nx, set_nx)
    if _newclass: nu = property(get_nu, set_nu)
    if _newclass: nb = property(get_nb, set_nb)
    if _newclass: nc = property(get_nc, set_nc)
%}
}

%extend qp_solver {
    qp_solver() {
        qp_solver *solver = (qp_solver *) malloc(sizeof(qp_solver));
        return solver;
    }
}
