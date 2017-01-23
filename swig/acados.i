%module acados
%{
#define SWIG_FILE_WITH_INIT
/* Put header files here or function declarations like below */
#include "swig/test.h"
%}

%include "numpy.i"
%init %{
import_array();
%}

%typemap (argout) my_struct * {
  /* some_int */
  {
    PyObject *member = PyInt_FromLong($1->some_int);
    if (!member) SWIG_fail;
    $result = member;
  }
  /* some_int_vector */
  {
    npy_intp dims[1] = { $1->some_int };
    PyObject *array = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, (void*)($1->some_int_vector));
    if (!array) SWIG_fail;
    $result = SWIG_Python_AppendOutput($result, array);
  }
}

%include "swig/test.h"
