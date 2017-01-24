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
    if (!PySequence_Check($input)) {
        PyErr_SetString(PyExc_ValueError,"Expected a sequence");
        SWIG_fail;
    }
    $1 = (int *) malloc(sizeof(int));
    $1[0] = 55;
    printf("%p\n", $1);
}

%typemap (out) int * {
    $result = PyInt_FromLong(*$1);
}

%include "swig/test.h"
