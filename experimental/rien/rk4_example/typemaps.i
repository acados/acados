%{
 /* Put header files here or function declarations like below */
#include "common_header.h"
%}

/* FLOAT TYPEMAPS */
%typemap(memberin) float [ANY] {
  int i;
  for (i = 0; i < $1_dim0; i++) {
      $1[i] = $input[i];
  }
}
%typemap(in) float [ANY] {
  int i;
  if (!PySequence_Check($input)) {
    PyErr_SetString(PyExc_ValueError,"Expected a sequence");
    return NULL;
  }
  if (PySequence_Length($input) != $1_dim0) {
    PyErr_SetString(PyExc_ValueError,"Size mismatch. Expected $1_dim0 elements");
    return NULL;
  }
  $1 = (float *) malloc($1_dim0*sizeof(float));
  for (i = 0; i < $1_dim0; i++) {
    PyObject *o = PySequence_GetItem($input,i);
    if (PyNumber_Check(o)) {
      $1[i] = (float) PyFloat_AsDouble(o);
    } else {
      PyErr_SetString(PyExc_ValueError,"Sequence elements must be numbers");      
      free($1);
      return NULL;
    }
  }
}
%typemap(freearg) float [ANY] {
   if ($1) free($1);
}
%typemap(out) float [ANY] {
  int i;
  $result = PyList_New($1_dim0);
  for (i = 0; i < $1_dim0; i++) {
    PyObject *o = PyFloat_FromDouble((double) $1[i]);
    PyList_SetItem($result,i,o);
  }
}


/* DOUBLE TYPEMAPS */
%typemap(memberin) double [ANY] {
  int i;
  for (i = 0; i < $1_dim0; i++) {
      $1[i] = $input[i];
  }
}
%typemap(in) double [ANY] {
  int i;
  if (!PySequence_Check($input)) {
    PyErr_SetString(PyExc_ValueError,"Expected a sequence");
    return NULL;
  }
  if (PySequence_Length($input) != $1_dim0) {
    PyErr_SetString(PyExc_ValueError,"Size mismatch. Expected $1_dim0 elements");
    return NULL;
  }
  $1 = (double *) malloc($1_dim0*sizeof(double));
  for (i = 0; i < $1_dim0; i++) {
    PyObject *o = PySequence_GetItem($input,i);
    if (PyNumber_Check(o)) {
      $1[i] = (double) PyFloat_AsDouble(o);
    } else {
      PyErr_SetString(PyExc_ValueError,"Sequence elements must be numbers");      
      free($1);
      return NULL;
    }
  }
}
%typemap(freearg) double [ANY] {
   if ($1) free($1);
}
%typemap(out) double [ANY] {
  int i;
  $result = PyList_New($1_dim0);
  for (i = 0; i < $1_dim0; i++) {
    PyObject *o = PyFloat_FromDouble((double) $1[i]);
    PyList_SetItem($result,i,o);
  }
}

