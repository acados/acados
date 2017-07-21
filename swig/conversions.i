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

#if defined(SWIGPYTHON)
%include "numpy.i"
%fragment("NumPy_Fragments");
%init %{
import_array();
%}
#endif

#if defined(SWIGPYTHON)
%pythoncode %{
from numpy import copy, copyto

class sequence_of_arrays(list):
    def __init__(self):
        super().__init__()
    def __init__(self, iterable):
        super().__init__(iterable)
    def __getitem__(self, k):
        return copy(super().__getitem__(k))
    def __setitem__(self, k, v):
        try:
            indices_to_set = range(*k.indices(len(self)))
        except AttributeError:
            # k is probably an integer
            try:
                indices_to_set = [int(k)]
            except TypeError:
                # k is probably tuple
                indices_to_set = k
        for index in indices_to_set:
            copyto(super().__getitem__(index), v)

%}

%{
// Global variable for Python module
PyTypeObject tuple_type;
PyObject *pModule = NULL;
%}
#endif

%{
#include "swig/conversions.h"

#include <algorithm>
#include <stdexcept>
#include <typeinfo>

#if defined(SWIGMATLAB)
template<typename T>
mxClassID get_numeric_type() {
    if (typeid(T) == typeid(real_t))
        return mxDOUBLE_CLASS;
    else if (typeid(T) == typeid(int_t))
        return mxDOUBLE_CLASS;
    throw std::invalid_argument("Matrix can only have integer or floating point entries");
    return mxUNKNOWN_CLASS;
}
#elif defined(SWIGPYTHON)
template<typename T>
int get_numeric_type() {
    if (typeid(T) == typeid(int_t))
        return NPY_INT32;
    else if (typeid(T) == typeid(long))
        return NPY_INT64;
    else if (typeid(T) == typeid(real_t))
        return NPY_DOUBLE;
    throw std::invalid_argument("Matrix can only have integer or floating point entries");
    return 0;
}
#endif

bool is_integer(const LangObject *input) {
#if defined(SWIGMATLAB)
    if (!mxIsScalar(input) || !mxIsNumeric(input))
        return false;
    return true;
#elif defined(SWIGPYTHON)
    if (!PyLong_Check((PyObject *) input))
        return false;
    return true;
#endif
}

int_t int_from(const LangObject *scalar) {
#if defined(SWIGMATLAB)
    return (int_t) mxGetScalar(scalar);
#elif defined(SWIGPYTHON)
    return (int_t) PyLong_AsLong((PyObject *) scalar);
#endif
}

LangObject *to_scalar(int_t scalar) {
#if defined(SWIGMATLAB)
    return mxCreateDoubleScalar((double) scalar);
#elif defined(SWIGPYTHON)
    return PyLong_FromLong((long) scalar);
#endif
}

bool is_matrix(const LangObject *input) {
#if defined(SWIGMATLAB)
    if (!mxIsNumeric(input))
        return false;
    mwSize nb_dims = mxGetNumberOfDimensions(input);
    if (nb_dims != 2)
        return false;
    return true;
#elif defined(SWIGPYTHON)
    if (!PyArray_Check(input))
        return false;
    int nb_dims = PyArray_NDIM((PyArrayObject *) input);
    if (nb_dims < 1 || nb_dims > 2)
        return false;
    return true;
#endif
}

bool is_matrix(const LangObject *input, const int_t nb_rows, const int_t nb_columns) {
    if (!is_matrix(input))
        return false;
#if defined(SWIGMATLAB)
    const mwSize *dims = mxGetDimensions(input);
    int_t input_rows = dims[0], input_cols = dims[1];
    if (input_rows != nb_rows || input_cols != nb_columns)
        return false;
    return true;
#elif defined(SWIGPYTHON)
    int nb_dims = PyArray_NDIM((PyArrayObject *) input);
    npy_intp *dims = PyArray_DIMS((PyArrayObject *) input);
    if (dims[0] != nb_rows)
        return false;
    if (nb_dims == 1) {
        if (nb_columns != 1)
            return false;
    } else {
        if (dims[1] != nb_columns)
            return false;
    }
    return true;
#endif
}

template<typename T>
LangObject *new_matrix(const int_t *dims, const T *data) {
    int_t nb_rows = dims[0];
    int_t nb_cols = dims[1];
#if defined(SWIGMATLAB)
    mxArray *matrix = mxCreateNumericMatrix(nb_rows, nb_cols, get_numeric_type<T>(), mxREAL);
    T *new_array = (T *) mxCalloc(nb_rows*nb_cols, sizeof(T));
    for (int_t i = 0; i < nb_rows*nb_cols; i++)
        new_array[i] = data[i];
    mxSetData(matrix, new_array);
    return matrix;
#elif defined(SWIGPYTHON)
    PyObject *matrix = NULL;
    if (nb_cols == 1) {
        npy_intp npy_dims[1] = {nb_rows};
        matrix = PyArray_NewFromDataF(1, npy_dims, get_numeric_type<T>(), (void *) data);
    } else {
        npy_intp npy_dims[2] = {nb_rows, nb_cols};
        matrix = PyArray_NewFromDataF(2, npy_dims, get_numeric_type<T>(), (void *) data);
    }
    if (matrix == NULL)
        throw std::runtime_error("Something went wrong while copying array");
    return matrix;
#endif
}

bool is_sequence(const LangObject *object) {
#if defined(SWIGMATLAB)
    if (!mxIsCell(object))
        return false;
#elif defined(SWIGPYTHON)
    if (!PyList_Check((PyObject *) object))
        return false;
#endif
    return true;
}

bool is_sequence(const LangObject *input, int_t expected_length) {
    if (!is_sequence(input))
        return false;
#if defined(SWIGMATLAB)
    int_t length_of_sequence = mxGetNumberOfElements(input);
#elif defined(SWIGPYTHON)
    int_t length_of_sequence = PyList_Size((PyObject *) input);
#endif
    if (length_of_sequence != expected_length)
        return false;
    return true;
}

LangObject *from(const LangObject *sequence, int_t index) {
#if defined(SWIGMATLAB)
    return mxGetCell(sequence, index);
#elif defined(SWIGPYTHON)
    return PyList_GetItem((PyObject *) sequence, index);
#endif
}

LangObject *new_sequence(const int_t length) {
#if defined(SWIGMATLAB)
    const mwSize dims[1] = {(const mwSize) length};
    return mxCreateCellArray(1, dims);
#elif defined(SWIGPYTHON)
    return PyList_New(length);
#endif
}

template <typename T>
LangObject *new_sequence_from(const T *array, const int_t length) {
    LangObject *sequence = new_sequence(length);
    for (int_t index = 0; index < length; index++) {
        if (typeid(T) == typeid(int_t))
            write_int_to(sequence, index, array[index]);
        else if (typeid(T) == typeid(real_t))
            write_real_to(sequence, index, array[index]);
    }
    return sequence;
}

void fill_int_array_from(const LangObject *sequence, int_t *array, const int_t length) {
    for (int_t index = 0; index < length; index++) {
        LangObject *item = from(sequence, index);
        if (!is_integer(item)) {
            char err_msg[256];
            snprintf(err_msg, sizeof(err_msg), "Input %s elements must be scalars",
                LANG_SEQUENCE_NAME);
            throw std::invalid_argument(err_msg);
        }
        array[index] = int_from(item);
    }
}

bool is_map(const LangObject *object) {
#if defined(SWIGMATLAB)
    if (!mxIsStruct(object))
        return false;
#elif defined(SWIGPYTHON)
    if (!PyDict_Check(object))
        return false;
#endif
    return true;
}

bool has(const LangObject *map, const char *key) {
#if defined(SWIGMATLAB)
    if (mxGetField(map, 0, key) == NULL)
        return false;
#elif defined(SWIGPYTHON)
    if (PyDict_GetItemString((PyObject *) map, key) == NULL)
        return false;
#endif
    return true;
}

LangObject *from(const LangObject *map, const char *key) {
    if (!has(map, key)) {
        char err_msg[256];
        snprintf(err_msg, sizeof(err_msg), "Input %s has no key %s", LANG_MAP_NAME, key);
        throw std::invalid_argument(err_msg);
    }
#if defined(SWIGMATLAB)
    return mxGetField(map, 0, key);
#elif defined(SWIGPYTHON)
    return PyDict_GetItemString((PyObject *) map, key);
#endif
}

int_t int_from(const LangObject *map, const char *key) {
#if defined(SWIGMATLAB)
    mxArray *value_ptr = mxGetField(map, 0, key);
    return (int_t) mxGetScalar(value_ptr);
#elif defined(SWIGPYTHON)
    return (int_t) PyLong_AsLong(PyDict_GetItemString((PyObject *) map, key));
#endif
}

void fill_array_from(const LangObject *input, int_t *array, const int_t length) {
    if (is_integer(input)) {
        int_t number = int_from(input);
        for (int_t i = 0; i < length; i++)
            array[i] = number;
    } else if (is_sequence(input, length)) {
        fill_int_array_from(input, array, length);
    } else {
        char err_msg[256];
        snprintf(err_msg, sizeof(err_msg), \
            "Expected scalar or %s of length %d", LANG_SEQUENCE_NAME, length);
        throw std::invalid_argument(err_msg);
    }
}

void fill_array_from(const LangObject *map, const char *key, int_t *array, int_t array_length) {
    if (!has(map, key)) {
        memset(array, 0, array_length*sizeof(*array));
    } else {
        LangObject *item = from(map, key);
        fill_array_from(item, array, array_length);
    }
}

void to(LangObject *sequence, const int_t index, LangObject *item) {
#if defined(SWIGMATLAB)
    mxSetCell(sequence, index, item);
#elif defined(SWIGPYTHON)
    PyList_SetItem(sequence, index, item);
#endif
}

void write_int_to(LangObject *sequence, const int_t index, const int_t number) {
#if defined(SWIGMATLAB)
    mxArray *scalar = mxCreateDoubleScalar(number);
    to(sequence, index, scalar);
#elif defined(SWIGPYTHON)
    to(sequence, index, PyLong_FromLong((long) number));
#endif
}

void write_real_to(LangObject *sequence, const int_t index, const real_t number) {
#if defined(SWIGMATLAB)
    mxArray *scalar = mxCreateDoubleScalar(number);
    to(sequence, index, scalar);
#elif defined(SWIGPYTHON)
    to(sequence, index, PyFloat_FromDouble((double) number));
#endif
}

LangObject *new_sequence_of_arrays(const int_t length) {
#if defined(SWIGMATLAB)
    mxArray *sequence = new_sequence(length);
#elif defined(SWIGPYTHON)
    // Try loading Python module into global variable
    if (pModule == NULL)
        pModule = PyImport_Import(PyString_FromString("acados"));
    // Check if loading was succesful
    if (pModule == NULL)
        SWIG_Error(SWIG_RuntimeError, "Something went wrong when importing Python module");
    PyObject *pDict = PyModule_GetDict(pModule);
    PyObject *pClass = PyDict_GetItemString(pDict, "sequence_of_arrays");
    PyObject *sequence = NULL;
    if (PyCallable_Check(pClass)) {
        PyObject *args = PyTuple_New(1);
        PyObject *list = PyList_New(length);
        for (int_t index = 0; index < length; index++)
            PyList_SetItem(list, index, PyLong_FromLong((long) index));  // fill list with dummies
        PyTuple_SetItem(args, 0, list);
        sequence = PyObject_CallObject(pClass, args);
        if (sequence == NULL)
            PyErr_Print();
    } else {
        char err_msg[256];
        snprintf(err_msg, sizeof(err_msg), "Something went wrong during construction of %s "
            "with length %d", LANG_SEQUENCE_NAME, length);
        SWIG_Error(SWIG_RuntimeError, err_msg);
    }
#endif
    if (sequence == NULL) {
        char err_msg[256];
        snprintf(err_msg, sizeof(err_msg), "Something went wrong during construction of %s "
            "with length %d", LANG_SEQUENCE_NAME, length);
        SWIG_Error(SWIG_RuntimeError, err_msg);
    }
    return sequence;
}

template<typename T>
LangObject *new_sequence_from(const T **data, const int_t length,
    const int_t *nb_rows, const int_t *nb_columns) {

    LangObject *sequence = new_sequence_of_arrays(length);
    for (int_t index = 0; index < length; index++) {
        int_t dims[2] = {nb_rows[index], nb_columns[index]};
        LangObject *item = new_matrix<T>(dims, data[index]);
        to(sequence, index, item);
    }
    return sequence;
}

template<typename T>
LangObject *new_sequence_from(const T **data, const int_t length,
    const int_t *nb_elems) {

    int_t nb_columns[length];
    for (int_t i = 0; i < length; i++)
        nb_columns[i] = 1;
    return new_sequence_from(data, length, nb_elems, nb_columns);
}

bool dimensions_match(const LangObject *matrix, const int_t *nb_rows, const int_t *nb_cols,
    const int_t length) {

    int_t rows = nb_rows[0];
    int_t cols = nb_cols[0];
    for (int_t i = 1; i < length; i++) {
        if (nb_rows[i] != rows || nb_cols[i] != cols) {
            throw std::invalid_argument("If just given one matrix, dimensions for all stages "
                "must be equal");
            return false;
        }
    }
    if (!is_matrix(matrix, rows, cols)) {
        throw std::invalid_argument("Input matrix has wrong dimensions");
        return false;
    }
    return true;
}

template<typename T>
void copy_from(const LangObject *matrix, T *data, const int_t nb_elems) {
#if defined(SWIGMATLAB)
    if (!mxIsDouble(matrix))
        throw std::invalid_argument("Only matrices with double precision numbers allowed");
    double *matrix_data = (double *) mxGetData(matrix);
    std::copy(matrix_data, matrix_data + nb_elems, data);
#elif defined(SWIGPYTHON)
    if (PyArray_TYPE((PyArrayObject *) matrix) == get_numeric_type<int_t>()) {
        int_t *matrix_data = (int_t *) PyArray_DATA((PyArrayObject *) matrix);
        std::copy(matrix_data, matrix_data + nb_elems, data);
    } else if (PyArray_TYPE((PyArrayObject *) matrix) == get_numeric_type<long>()) {
        long *matrix_data = (long *) PyArray_DATA((PyArrayObject *) matrix);
        std::copy(matrix_data, matrix_data + nb_elems, data);
    } else if (PyArray_TYPE((PyArrayObject *) matrix) == get_numeric_type<real_t>()) {
        real_t *matrix_data = (real_t *) PyArray_DATA((PyArrayObject *) matrix);
        std::copy(matrix_data, matrix_data + nb_elems, data);
    } else {
        throw std::invalid_argument("Only matrices with integer numbers or double "
            "precision numbers allowed");
    }
#endif
}

template<typename T>
void fill_array_from(const LangObject *input, T **array,
    const int_t length, const int_t *nb_rows, const int_t *nb_columns) {

    if (is_matrix(input) && dimensions_match(input, nb_rows, nb_columns, length)) {
        int_t nb_elems = nb_rows[0]*nb_columns[0];
        for (int_t index = 0; index < length; index++) {
            copy_from(input, array[index], nb_elems);
        }
    } else if (is_sequence(input, length)) {
        for (int_t index = 0; index < length; index++) {
            LangObject *item = from(input, index);
            if (is_matrix(item, nb_rows[index], nb_columns[index]))
                copy_from(item, array[index], nb_rows[index]*nb_columns[index]);
        }
    } else {
        char err_msg[256];
        snprintf(err_msg, sizeof(err_msg),
            "Expected %s or %s as input", LANG_SEQUENCE_NAME, LANG_MATRIX_NAME);
        throw std::invalid_argument(err_msg);
    }
}

template<typename T>
void fill_array_from(const LangObject *input, T **array, const int_t length,
    const int_t *nb_elems) {

    int_t nb_columns[length];
    for (int_t i = 0; i < length; i++)
        nb_columns[i] = 1;
    fill_array_from(input, array, length, nb_elems, nb_columns);
}

// TODO(roversch): This can probably be merged with the new_sequence_from functions.
LangObject *new_output_list_from(const LangObject **input, const int_t length) {
    LangObject *output_list = new_sequence(length);
    for (int_t index = 0; index < length; index++) {
        to(output_list, index, (LangObject *) input[index]);
    }
    return output_list;
}

LangObject *sequence_concatenate(const LangObject *seq1, const LangObject *seq2) {
#if defined(SWIGMATLAB)
    if (mxGetNumberOfDimensions(seq1) != 1 || mxGetNumberOfDimensions(seq2) != 1)
        throw std::invalid_argument("Can only concatenate 1-D cell arrays");
    int_t length_seq1 = mxGetNumberOfElements(seq1);
    int_t length_seq2 = mxGetNumberOfElements(seq2);
    int_t total_length = length_seq1 + length_seq2;
    const mwSize dims[1] = {(const mwSize) total_length};
    mxArray *output_array = mxCreateCellArray(1, dims);
    for (int_t index = 0; index < length_seq1; index++)
        mxSetCell(output_array, index, mxGetCell(seq1, index));
    for (int_t index = 0; index < length_seq2; index++)
        mxSetCell(output_array, length_seq1 + index, mxGetCell(seq2, index));
    return output_array;
#elif defined(SWIGPYTHON)
    return PySequence_Concat((PyObject *) seq1, (PyObject *) seq2);
#endif
}

bool is_named_tuple(const LangObject *object) {
#if defined(SWIGMATLAB)
    return mxIsStruct(object);
#elif defined(SWIGPYTHON)
    return PyTuple_Check((PyObject *) object);
#endif
}

LangObject *new_states_controls_output_tuple(LangObject *states, LangObject *controls) {
    const char *fieldnames[2] = {"states", "controls"};
#if defined(SWIGMATLAB)
    const mwSize dims[1] = {(const mwSize) 1};
    mxArray *named_tuple = mxCreateStructArray(1, dims, 2, fieldnames);
    mxSetField(named_tuple, 0, fieldnames[0], states);
    mxSetField(named_tuple, 0, fieldnames[1], controls);
    return named_tuple;
#elif defined(SWIGPYTHON)
    // The list of field names in named tuples must be NULL-terminated in Python
    PyStructSequence_Field fields[3];
    fields[0].name = (char *) fieldnames[0];
    fields[0].doc = NULL;
    fields[1].name = (char *) fieldnames[1];
    fields[1].doc = NULL;
    fields[2].name = NULL;
    fields[2].doc = NULL;
    PyStructSequence_Desc tuple_descriptor;
    tuple_descriptor.name = (char *) "output";
    tuple_descriptor.doc = NULL;
    tuple_descriptor.fields = fields;
    tuple_descriptor.n_in_sequence = 2;
    PyStructSequence_InitType2(&tuple_type, &tuple_descriptor);
    tuple_type.tp_flags = tuple_type.tp_flags | Py_TPFLAGS_HEAPTYPE;
    PyObject *named_tuple = PyStructSequence_New(&tuple_type);
    PyStructSequence_SetItem(named_tuple, 0, (PyObject *) states);
    PyStructSequence_SetItem(named_tuple, 1, (PyObject *) controls);
    return named_tuple;
#endif
}

%}
