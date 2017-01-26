%module acados
%{
#define SWIG_FILE_WITH_INIT
#define PyArray_SimpleNewFromDataF(nd, dims, typenum, data) \
        PyArray_New(&PyArray_Type, nd, dims, typenum, NULL, \
                    data, 0, NPY_ARRAY_FARRAY, NULL)
#include "acados/utils/types.h"
#include "acados/ocp_qp/ocp_qp_common.h"
%}

%include "numpy.i"
%init %{
import_array();
%}

%include "acados/utils/types.h"

%{
static PyObject *convert_to_sequence(const int_t *array, const int_t length) {
    PyObject *sequence = PyList_New(length);
    for (int_t i = 0; i < length; i++) {
        PyList_SetItem(sequence, i, PyInt_FromLong((long) array[i]));
    }
    return sequence;
}
%}

%{
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
%}

%{
static int_t convert_to_c_array(PyObject *input, int_t * const array, const int_t length_of_array) {
    if (PyNumber_Check(input)) {
        int_t integer_number = (int_t) PyInt_AsLong(input);
        for (int_t i = 0; i < length_of_array; i++)
            array[i] = integer_number;
    } else if (PySequence_Check(input)) {
        int_t length_of_sequence = PySequence_Length(input);
        if (length_of_sequence != length_of_array) {
            char err_msg[256];
            snprintf(err_msg, sizeof(err_msg), "Length of sequence must be %d", length_of_array);
            SWIG_Error(SWIG_ValueError, err_msg);
            return 0;
        }
        for (int_t i = 0; i < length_of_array; i++) {
            PyObject *o = PySequence_GetItem(input, i);
            if (!PyInt_Check(o)) {
                Py_XDECREF(o);
                SWIG_Error(SWIG_ValueError, "Sequence elements must be integer numbers");
                return 0;
            }
            array[i] = (int_t) PyInt_AsLong(o);
            Py_DECREF(o);
        }
    } else {
        SWIG_Error(SWIG_ValueError, "Expected a sequence or an integer number");
        return 0;
    }
    return 1;
}
%}

%{
static int_t convert_to_2dim_c_array(PyObject *input, real_t ** const array, const int_t length_of_array, const int_t *sizes) {
    if (PyArray_Check(input)) {
        PyArrayObject *input_array = reinterpret_cast<PyArrayObject *>(PyArray_FROM_OTF(input, NPY_DOUBLE, NPY_ARRAY_F_CONTIGUOUS));
        if (PyArray_NDIM(input_array) != 2) {
            SWIG_Error(SWIG_ValueError, "Expected a 2D array as input");
            return 0;
        }
        npy_intp *dims = PyArray_DIMS(input_array);
        int_t size_of_all_dimensions = sizes[0];
        if (dims[0] != size_of_all_dimensions || dims[1] != size_of_all_dimensions) {
            SWIG_Error(SWIG_ValueError, "Input array with wrong dimensions");
            return 0;
        }
        printf("PyArray type: %d\n", PyArray_TYPE(input_array));
        npy_intp size = PyArray_SIZE(input_array);
        printf("PyArray size: %lu\n", size);
        for (int_t i = 0; i < length_of_array; i++) {
            npy_intp dims[2] = {0, 0};
            array[i] = (real_t *) PyArray_GetPtr(input_array, dims);
        }
    } else {
        SWIG_Error(SWIG_ValueError, "Expected an array as input");
        return 0;
    }
    return 1;
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
            SWIG_Error(SWIG_ValueError, "Input must be a dictionary");
            return NULL;
        }
        ocp_qp_in *qp = (ocp_qp_in *) malloc(sizeof(ocp_qp_in));
        PyObject *num_stages = PyDict_GetItemString(dictionary, "N");
        if (!PyInt_Check(num_stages)) {
            SWIG_Error(SWIG_ValueError, "N must be an integer number");
            return NULL;
        }
        qp->N = (int_t) PyInt_AsLong(num_stages);
        qp->nx = (int_t *) calloc(qp->N+1, sizeof(int_t));
        PyObject *nx_list = PyDict_GetItemString(dictionary, "nx");
        if (!convert_to_c_array(nx_list, (int_t *) qp->nx, qp->N+1)) {
            SWIG_Error(SWIG_ValueError, "nx is not a valid input");
            return NULL;
        }
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
