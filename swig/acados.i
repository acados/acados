%module acados
%{
#define SWIG_FILE_WITH_INIT
/* Put header files here or function declarations like below */
#include "swig/test.h"
%}

%include "typemaps.i"
%include "numpy.i"
%init %{
import_array();
%}

%typemap (in) int * {
    if (PyNumber_Check($input)) {
        printf("Reading one number\n");
        $1 = (int *) malloc(sizeof(int));
        $1[0] = (int) PyInt_AsLong($input);
    }
    else if (PySequence_Check($input)) {
        int length = PySequence_Length($input);
        printf("Reading sequence of length %d\n", length);
        $1 = (int *) malloc(length * sizeof(int));
        int i;
        for (i = 0; i < length; i++) {
            PyObject *o = PySequence_GetItem($input, i);
            if (PyNumber_Check(o)) {
                $1[i] = (float) PyFloat_AsDouble(o);
            } else {
                PyErr_SetString(PyExc_ValueError,"Sequence elements must be numbers");
                SWIG_fail;
            }
        }
    } else {
        PyErr_SetString(PyExc_ValueError,"Expected a sequence or a number");
        SWIG_fail;
    }
}

%typemap (out) int * {
    $result = PyInt_FromLong(*$1);
}


%include "swig/test.h"
