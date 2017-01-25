%module acados
%{
#define SWIG_FILE_WITH_INIT
/* Put header files here or function declarations like below */
#include "acados/utils/types.h"
#include "acados/ocp_qp/ocp_qp_common.h"
%}

%include "typemaps.i"
%include "numpy.i"
%init %{
import_array();
%}

%include "acados/utils/types.h"

%{
static PyObject *convert_to_sequence(const int_t *array, const int_t length_of_array) {
    PyObject *sequence = PyList_New(length_of_array);
    for (int_t i = 0; i < length_of_array; i++) {
        PyList_SetItem(sequence, i, PyInt_FromLong((long) array[i]));
    }
    return sequence;
}
%}

%{
static int_t convert_to_array(PyObject *input, int_t * const array, const int length_of_array) {
    if (PyNumber_Check(input)) {
        int_t integer_number = (int_t) PyInt_AsLong(input);
        for (int_t i = 0; i < length_of_array; i++)
            array[i] = integer_number;
    }
    else if (PySequence_Check(input)) {
        int_t length_of_sequence = PySequence_Length(input);
        if (length_of_sequence != length_of_array) {
            char err_msg[256];
            sprintf(err_msg, "Length of sequence must be %d", length_of_array);
            SWIG_Error(SWIG_ValueError, err_msg);
            return 1;
        }
        for (int_t i = 0; i < length_of_array; i++) {
            PyObject *o = PySequence_GetItem(input, i);
            if (!PyInt_Check(o)) {
                Py_XDECREF(o);
                SWIG_Error(SWIG_ValueError, "Sequence elements must be integer numbers");
                return 1;
            }
            array[i] = (int_t) PyInt_AsLong(o);
            Py_DECREF(o);
        }
    } else {
        SWIG_Error(SWIG_ValueError, "Expected a sequence or an integer number");
        return 1;
    }
    return 0;
}
%}

%include "acados/ocp_qp/ocp_qp_common.h"

%extend ocp_qp_in {
    ocp_qp_in(int_t N) {
        ocp_qp_in *qp_in;
        qp_in = (ocp_qp_in *) malloc(sizeof(ocp_qp_in));
        qp_in->N = N;
        qp_in->nx = (int_t *) calloc(N+1, sizeof(int_t));
        qp_in->nu = (int_t *) calloc(N, sizeof(int_t));
        qp_in->nb = (int_t *) calloc(N+1, sizeof(int_t));
        qp_in->nc = (int_t *) calloc(N+1, sizeof(int_t));
        return qp_in;
    }
    PyObject *get_nx() {
        return convert_to_sequence($self->nx, $self->N+1);
    }
    void set_nx(PyObject *input) {
        convert_to_array(input, (int_t *) $self->nx, $self->N+1);
    }
    %pythoncode %{
        __swig_getmethods__["nx"] = get_nx
        __swig_setmethods__["nx"] = set_nx
        if _newclass: nx = property(get_nx, set_nx)
    %}
}
