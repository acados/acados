%module acados
%{
#define SWIG_FILE_WITH_INIT
#define PyArray_SimpleNewFromDataF(nd, dims, typenum, data) \
        PyArray_New(&PyArray_Type, nd, dims, typenum, NULL, \
                    data, 0, NPY_ARRAY_FARRAY, NULL)
#include <stdexcept>
#include "acados/utils/types.h"
#include "acados/ocp_qp/ocp_qp_common.h"
%}

%include "exception.i"

%exception {
    try {
        $action
    } catch (std::invalid_argument& e) {
        SWIG_exception(SWIG_ValueError, const_cast<char*>(e.what()));
    }
}

%include "numpy.i"
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


static PyObject *convert_to_sequence_of_arrays(const real_t **c_array,
    const int_t length, const int_t *sizes) {

    PyObject *sequence = PyList_New(length);
    for (int_t i = 0; i < length; i++) {
        npy_intp dims[2] = {sizes[i+1], sizes[i]};
        PyObject *py_array = PyArray_SimpleNewFromDataF(2, dims, NPY_DOUBLE, (void*) c_array[i]);
        PyList_SetItem(sequence, i, py_array);
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
        throw std::invalid_argument("Expected a sequence or an integer number");
    }
}

static PyArrayObject *array_from_object(PyObject *input, int_t dim1, int_t dim2) {
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

static void convert_to_2dim_c_array(PyObject * const input, real_t ** const array,
    const int_t length_of_array, const int_t *sizes) {

    if (is_valid_2dim_array(input)) {
        int_t dim = sizes[0];
        PyArrayObject *input_array = array_from_object(input, dim, dim);
        npy_intp element[2] = {0, 0};
        for (int_t i = 0; i < length_of_array; i++) {
            array[i] = (real_t *) PyArray_GetPtr(input_array, element);
        }
    } else if (is_sequence_with_length(input, length_of_array)) {
        for (int_t i = 0; i < length_of_array; i++) {
            PyObject *item = PySequence_GetItem(input, i);
            npy_intp element[2] = {0, 0};
            if (is_valid_2dim_array(item)) {
                PyArrayObject *input_array = array_from_object(item, sizes[i], sizes[i+1]);
                array[i] = (real_t *) PyArray_GetPtr(input_array, element);
            }
        }
    } else {
        throw std::invalid_argument("Expected an array as input");
    }
}

%}

%include "acados/ocp_qp/ocp_qp_common.h"

%extend ocp_qp_in {
    ocp_qp_in(int_t num_stages) {
        ocp_qp_in *qp = (ocp_qp_in *) malloc(sizeof(ocp_qp_in));
        qp->N = num_stages;
        qp->nx = (int_t *) calloc(num_stages+1, sizeof(int_t));
        return qp;
    }
    ocp_qp_in(PyObject *dictionary) {
        if (!PyDict_Check(dictionary)) {
            throw std::invalid_argument("Input must be a dictionary");
        }
        ocp_qp_in *qp = (ocp_qp_in *) malloc(sizeof(ocp_qp_in));
        PyObject *num_stages = PyDict_GetItemString(dictionary, "N");
        if (!is_valid_integer(num_stages)) {
            throw std::invalid_argument("N must be an integer number");
        }
        qp->N = (int_t) PyInt_AsLong(num_stages);
        qp->nx = (int_t *) calloc(qp->N+1, sizeof(int_t));
        PyObject *nx_list = PyDict_GetItemString(dictionary, "nx");
        convert_to_c_array(nx_list, (int_t *) qp->nx, qp->N+1);
        qp->A = (const real_t **) calloc(qp->N, sizeof(*qp->A));
        for (int_t i = 0; i < qp->N; i++) {
            real_t *temp = (real_t *) calloc(qp->nx[i]*qp->nx[i+1], sizeof(real_t));
            temp[1] = 1;
            temp[2] = 2;
            qp->A[i] = temp;
        }
        return qp;
    }
    PyObject *get_nx() {
        return convert_to_sequence($self->nx, $self->N+1);
    }
    void set_nx(PyObject *input) {
        convert_to_c_array(input, (int_t *) $self->nx, $self->N+1);
    }
    PyObject *get_A() {
        return convert_to_sequence_of_arrays($self->A, $self->N, $self->nx);
    }
    void set_A(PyObject *input) {
        convert_to_2dim_c_array(input, (real_t **) $self->A, $self->N, $self->nx);
    }
    %pythoncode %{
        __swig_getmethods__["nx"] = get_nx
        __swig_setmethods__["nx"] = set_nx
        __swig_getmethods__["A"] = get_A
        __swig_setmethods__["A"] = set_A
        if _newclass: nx = property(get_nx, set_nx)
        if _newclass: A = property(get_A, set_A)
    %}
}
