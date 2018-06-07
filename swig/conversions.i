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
// %pythoncode %{
// from numpy import copy, copyto

// class sequence_of_arrays(list):
//     def __init__(self):
//         super().__init__()
//     def __init__(self, iterable):
//         super().__init__(iterable)
//     def __getitem__(self, k):
//         return copy(super().__getitem__(k))
//     def __setitem__(self, k, v):
//         try:
//             indices_to_set = range(*k.indices(len(self)))
//         except AttributeError:
//             # k is probably an integer
//             try:
//                 indices_to_set = [int(k)]
//             except TypeError:
//                 # k is probably tuple
//                 indices_to_set = k
//         for index in indices_to_set:
//             copyto(super().__getitem__(index), v)

// %}

%{
// Global variable for Python module
PyTypeObject *tuple_type = NULL;
PyStructSequence_Desc *tuple_descriptor = NULL;
PyObject *sequence_of_arrays_module = NULL;
PyObject *copy_module = NULL;

%}
#endif

%{
#include "swig/conversions.h"

#include <algorithm>
#include <stdexcept>
#include <utility>
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
    double input_as_double = mxGetScalar(input);
    int input_as_int = mxGetScalar(input);
    if ((double) input_as_int != input_as_double)
        return false;
    return true;
#elif defined(SWIGPYTHON)
    if (!PyLong_Check((PyObject *) input))
        return false;
    return true;
#endif
}

bool is_real(const LangObject *input) {
#if defined(SWIGMATLAB)
    if (!mxIsScalar(input) || !mxIsNumeric(input))
        return false;
    return true;
#elif defined(SWIGPYTHON)
    if (!PyFloat_Check((PyObject *) input))
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

real_t real_from(const LangObject *scalar) {
#if defined(SWIGMATLAB)
    return (real_t) mxGetScalar(scalar);
#elif defined(SWIGPYTHON)
    return (real_t) PyFloat_AsDouble((PyObject *) scalar);
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

int numRows(const LangObject *input) {
    if (!is_matrix(input))
        throw std::invalid_argument("Input is not a valid matrix.");
#if defined(SWIGMATLAB)
    const mwSize *dims = mxGetDimensions(input);
    return dims[0];
#elif defined(SWIGPYTHON)
    npy_intp *dims = PyArray_DIMS((PyArrayObject *) input);
    return dims[0];
#endif
}

int numColumns(const LangObject *input) {
    if (!is_matrix(input))
        throw std::invalid_argument("Input is not a valid matrix.");
#if defined(SWIGMATLAB)
    const mwSize *dims = mxGetDimensions(input);
    return dims[1];
#elif defined(SWIGPYTHON)
    int nb_dims = PyArray_NDIM((PyArrayObject *) input);
    npy_intp *dims = PyArray_DIMS((PyArrayObject *) input);
    if (nb_dims == 1) {
        // column vector
        return 1;
    } else {
        return dims[1];
    }
#endif
}

bool is_matrix(const LangObject *input, const int_t nb_rows, const int_t nb_columns) {
    if (!is_matrix(input))
        return false;
    if (nb_rows != numRows(input) || nb_columns != numColumns(input))
        return false;
    return true;
}

double *asDoublePointer(LangObject *input) {
    if (!is_matrix(input))
        throw std::invalid_argument("Input is not of a valid matrix type.");
#if defined(SWIGMATLAB)
    return (double *) mxGetData(input);
#elif defined(SWIGPYTHON)
    PyObject *matrix = PyArray_FROM_OTF(input, NPY_FLOAT64, NPY_ARRAY_FARRAY_RO);
    if (matrix == NULL) {
        PyErr_Print();
        throw std::runtime_error("Something went wrong while converting matrix");
    }
    return (double *) PyArray_DATA((PyArrayObject *) matrix);
#endif
}

template<typename T>
LangObject *new_matrix(std::pair<int, int> dimensions, const T *data) {
    int_t nb_rows = dimensions.first;
    int_t nb_cols = dimensions.second;
#if defined(SWIGMATLAB)
    mxArray *matrix = mxCreateNumericMatrix(nb_rows, nb_cols, get_numeric_type<T>(), mxREAL);
    double *new_array = (double *) mxCalloc(nb_rows*nb_cols, sizeof(double));
    for (int_t i = 0; i < nb_rows*nb_cols; i++)
        new_array[i] = (double) data[i];
    mxSetData(matrix, new_array);
    return matrix;
#elif defined(SWIGPYTHON)
    PyObject *matrix = NULL;
    if (nb_cols == 1) {
        double *data_copy = (double *) calloc(nb_rows, sizeof(double));
        std::copy_n(data, nb_rows, data_copy);
        npy_intp npy_dims[1] = {nb_rows};
        matrix = PyArray_NewFromDataF(1, npy_dims, data_copy);
    } else {
        double *data_copy = (double *) calloc(nb_rows * nb_cols, sizeof(double));
        std::copy_n(data, nb_rows * nb_cols, data_copy);
        npy_intp npy_dims[2] = {nb_rows, nb_cols};
        matrix = PyArray_NewFromDataF(2, npy_dims, data_copy);
    }
    if (matrix == NULL) {
        PyErr_Print();
        throw std::runtime_error("Something went wrong while copying array");
    }
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
    PyObject *item = PyList_GetItem((PyObject *) sequence, index);
    Py_INCREF(item);
    return item;
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
LangObject *new_sequence_from(T *array, const int_t length) {
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
            char err_msg[MAX_STR_LEN];
            snprintf(err_msg, sizeof(err_msg), "Input %s elements must be scalars",
                LANG_SEQUENCE_NAME);
            throw std::invalid_argument(err_msg);
        }
        array[index] = int_from(item);
    }
}

void fill_real_array_from(const LangObject *sequence, real_t *array, const int_t length) {
    for (int_t index = 0; index < length; index++) {
        LangObject *item = from(sequence, index);
        if (!is_real(item)) {
            char err_msg[MAX_STR_LEN];
            snprintf(err_msg, sizeof(err_msg), "Input %s elements must be scalars",
                LANG_SEQUENCE_NAME);
            throw std::invalid_argument(err_msg);
        }
        array[index] = real_from(item);
    }
}

bool is_map(const LangObject *object) {
#if defined(SWIGMATLAB)
    if (!mxIsStruct(object))
        return false;
    if (mxGetNumberOfElements(object) != 1)
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

int num_elems(const LangObject *map) {
#if defined(SWIGMATLAB)
    return mxGetNumberOfFields(map);
#elif defined(SWIGPYTHON)
    return PyDict_Size((PyObject *) map);
#endif
}

LangObject *from(const LangObject *map, const char *key) {
    if (!has(map, key)) {
        char err_msg[MAX_STR_LEN];
        snprintf(err_msg, sizeof(err_msg), "Input %s has no key %s", LANG_MAP_NAME, key);
        throw std::invalid_argument(err_msg);
    }
#if defined(SWIGMATLAB)
    return mxGetField(map, 0, key);
#elif defined(SWIGPYTHON)
    PyObject *item = PyDict_GetItemString((PyObject *) map, key);
    if (item)
        Py_INCREF(item);
    return item;
#endif
}

const char *char_from(const LangObject *map, const char *key) {
    LangObject *value = from(map, key);
#if defined(SWIGMATLAB)
    return (const char *) mxArrayToString(value);
#elif defined(SWIGPYTHON)
    return (const char *) PyUnicode_AsUTF8AndSize(value, NULL);
#endif
}

int_t int_from(const LangObject *map, const char *key) {
    LangObject *value = from(map, key);
#if defined(SWIGMATLAB)
    return (int_t) mxGetScalar(value);
#elif defined(SWIGPYTHON)
    return (int_t) PyLong_AsLong(value);
#endif
}

real_t real_from(const LangObject *map, const char *key) {
    LangObject *value = from(map, key);
#if defined(SWIGMATLAB)
    return (real_t) mxGetScalar(value);
#elif defined(SWIGPYTHON)
    return (real_t) PyFloat_AsDouble(value);
#endif
}

bool is_string(LangObject *input) {
#if defined(SWIGMATLAB)
    return mxIsChar(input);
#elif defined(SWIGPYTHON)
    return PyUnicode_Check(input);
#endif
}

std::string string_from(LangObject *input) {
    const char *string;
#if defined(SWIGMATLAB)
    string = mxArrayToUTF8String(input);
    if (string == NULL)
        throw std::runtime_error("Error during string conversion.");
#elif defined(SWIGPYTHON)
    string = PyUnicode_AsUTF8(input);
#endif
    return std::string(string);
}

bool is_boolean(LangObject *input) {
#if defined(SWIGMATLAB)
    return mxIsLogicalScalar(input);
#elif defined(SWIGPYTHON)
    return PyBool_Check(input);
#endif
}

bool boolean_from(LangObject *input) {
#if defined(SWIGMATLAB)
    return mxIsLogicalScalarTrue(input);
#elif defined(SWIGPYTHON)
    return PyObject_IsTrue(input);
#endif
}

bool is_valid_option_type(LangObject *input) {
    return is_integer(input) || is_real(input) || is_matrix(input) || is_map(input)
                             || is_string(input) || is_boolean(input);
}

void to(LangObject *sequence, const int_t index, LangObject *item) {
#if defined(SWIGMATLAB)
    mxSetCell(sequence, index, item);
#elif defined(SWIGPYTHON)
    Py_INCREF(item);
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
    if (sequence_of_arrays_module == NULL)
        sequence_of_arrays_module = PyImport_Import(PyString_FromString("acados"));
    // Check if loading was succesful
    if (sequence_of_arrays_module == NULL)
        SWIG_Error(SWIG_RuntimeError, "Something went wrong when importing Python module");
    PyObject *pDict = PyModule_GetDict(sequence_of_arrays_module);
    PyObject *pClass = PyDict_GetItemString(pDict, "sequence_of_arrays");
    if (pClass)
        Py_INCREF(pClass);
    PyObject *sequence = NULL;
    if (PyCallable_Check(pClass)) {
        PyObject *args = PyTuple_New(1);
        PyObject *list = PyList_New(length);
        for (int_t index = 0; index < length; index++)
            PyList_SetItem(list, index, PyLong_FromLong((long) index));  // fill list with dummies
        Py_INCREF(list);
        PyTuple_SetItem(args, 0, list);
        sequence = PyObject_CallObject(pClass, args);
        Py_DECREF(pClass);
    }
#endif
    if (sequence == NULL) {
        char err_msg[MAX_STR_LEN];
        snprintf(err_msg, sizeof(err_msg), "Something went wrong during construction of %s "
            "with length %d", LANG_SEQUENCE_NAME, length);
        SWIG_Error(SWIG_RuntimeError, err_msg);
    }
    return sequence;
}

template<typename T>
LangObject *new_sequence_from(T **data, const int_t length,
    const int_t *nb_rows, const int_t *nb_columns) {

    LangObject *sequence = new_sequence_of_arrays(length);
    for (int_t index = 0; index < length; index++) {
        auto dims = std::make_pair(nb_rows[index], nb_columns[index]);
        LangObject *item = new_matrix<T>(dims, data[index]);
        to(sequence, index, item);
    }
    return sequence;
}

template<typename T>
LangObject *new_sequence_from(T **data, const int_t length,
    const int_t *nb_elems) {

    int_t *nb_columns = (int_t *) calloc(length, sizeof(int_t));
    for (int_t i = 0; i < length; i++)
        nb_columns[i] = 1;
    LangObject *result = new_sequence_from(data, length, nb_elems, nb_columns);
    free(nb_columns);
    return result;
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
        char err_msg[MAX_STR_LEN];
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

template<typename T>
void fill_array_from(const LangObject *input, T **array, const int_t length) {

    int_t nb[length];
    for (int_t i = 0; i < length; i++)
        nb[i] = 1;
    fill_array_from(input, array, length, nb, nb);
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

LangObject *new_output_tuple(int_t num_fields, const char **field_names, LangObject **content) {
#if defined(SWIGMATLAB)
    const mwSize dims[1] = {(const mwSize) 1};
    mxArray *named_tuple = mxCreateStructArray(1, dims, num_fields, field_names);
    for (int_t index = 0; index < num_fields; index++)
        mxSetField(named_tuple, 0, field_names[index], content[index]);
    return named_tuple;
#elif defined(SWIGPYTHON)
    PyObject *content_copy[num_fields];
    // Try loading Python module into global variable
    if (copy_module == NULL)
        copy_module = PyImport_Import(PyString_FromString("copy"));
    // Check if loading was succesful
    if (copy_module == NULL)
        SWIG_Error(SWIG_RuntimeError, "Something went wrong when importing Python module 'copy'");
    PyObject *pDict = PyModule_GetDict(copy_module);
    PyObject *pFunction = PyDict_GetItemString(pDict, "deepcopy");
    if (!PyCallable_Check(pFunction)) {
        SWIG_Error(SWIG_RuntimeError, "Function is not callable");
    }
    for (int_t index = 0; index < num_fields; index++) {
        PyObject *args = PyTuple_New(1);
        PyTuple_SetItem(args, 0, content[index]);
        content_copy[index] = PyObject_CallObject(pFunction, args);
    }

    PyStructSequence_Field *fields;
    fields = (PyStructSequence_Field *) calloc(num_fields+1, sizeof(PyStructSequence_Field));
    for (int_t index = 0; index < num_fields; index++) {
        fields[index].name = (char *) field_names[index];
        fields[index].doc = NULL;
    }
    // The list of field names in named tuples must be NULL-terminated in Python
    fields[num_fields].name = NULL;
    fields[num_fields].doc = NULL;

    tuple_descriptor = (PyStructSequence_Desc *) malloc(sizeof(PyStructSequence_Desc));
    tuple_descriptor->name = (char *) "output";
    tuple_descriptor->doc = NULL;
    tuple_descriptor->fields = fields;
    tuple_descriptor->n_in_sequence = num_fields;

    tuple_type = (PyTypeObject *) malloc(sizeof(PyTypeObject));
    PyStructSequence_InitType(tuple_type, tuple_descriptor);

    PyObject *named_tuple = PyStructSequence_New(tuple_type);
    for (int_t index = 0; index < num_fields; index++)
        PyStructSequence_SetItem(named_tuple, index, (PyObject *) content_copy[index]);
    return named_tuple;
#endif
}

LangObject *new_ocp_output_tuple(LangObject *states, LangObject *controls) {
    const char *field_names[2] = {"states", "controls"};
    LangObject *fields[2] = {states, controls};
    return new_output_tuple(2, field_names, fields);
}

LangObject *new_sim_output_tuple(LangObject *final_state, LangObject *forward_sensitivities) {
    const char *field_names[2] = {"final_state", "forward_sensitivities"};
    LangObject *fields[2] = {final_state, forward_sensitivities};
    return new_output_tuple(2, field_names, fields);
}

LangObject *new_ocp_nlp_function_output_tuple(
    LangObject * y, LangObject * jac_y, LangObject * hess_y) {
    const char *field_names[3] = {"y", "jac_y", "hess_y"};
    LangObject *fields[3] = {y, jac_y, hess_y};
    return new_output_tuple(3, field_names, fields);
}

void fill_array_from(const LangObject *input, int_t *array, const int_t length) {
    if (is_integer(input)) {
        int_t number = int_from(input);
        for (int_t i = 0; i < length; i++)
            array[i] = number;
    } else if (is_sequence(input, length)) {
        fill_int_array_from(input, array, length);
    } else {
        char err_msg[MAX_STR_LEN];
        snprintf(err_msg, sizeof(err_msg), \
            "Expected scalar or %s of length %d", LANG_SEQUENCE_NAME, length);
        throw std::invalid_argument(err_msg);
    }
}

void fill_array_from(const LangObject *input, real_t *array, const int_t length) {
    if (is_real(input)) {
        real_t number = real_from(input);
        for (int_t i = 0; i < length; i++)
            array[i] = number;
    } else if (is_sequence(input, length)) {
        fill_real_array_from(input, array, length);
    } else if (is_matrix(input, length, 1)) {
        copy_from(input, array, length);
    } else {
        char err_msg[MAX_STR_LEN];
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

#include "acados_cpp/options.hpp"

namespace acados {

template<typename T>
option_t *as_option_ptr(T val) {
    return new option<T>(val);
}

option_t *make_option_map(LangObject *val);

template<>
option_t *as_option_ptr(LangObject *val) {
    if (is_integer(val))
        return new option<int>(int_from(val));
    else if (is_real(val))
        return new option<double>(real_from(val));
    else if (is_boolean(val))
        return new option<bool>(boolean_from(val));
    else if (is_string(val))
        return new option<std::string>(string_from(val));
    else if (is_map(val))
        return make_option_map(val);
    else
        throw std::invalid_argument("Option does not have a valid type");
}

option_t *make_option_map(LangObject *val) {
    std::map<std::string, option_t *> option_map;
#if defined(SWIGMATLAB)
    int num_fields = mxGetNumberOfFields(val);
    for (int i = 0; i < num_fields; ++i) {
        std::string field_name {mxGetFieldNameByNumber(val, i)};
        option_map[field_name] = as_option_ptr(mxGetField(val, 0, field_name.c_str()));
    }
#elif defined(SWIGPYTHON)
    PyObject *key, *value;
    Py_ssize_t pos = 0;

    while (PyDict_Next(val, &pos, &key, &value)) {
        std::string field_name {PyUnicode_AsUTF8AndSize(key, NULL)};
        option_map[field_name] = as_option_ptr(value);
    }
#endif
    return new option<std::map<std::string, option_t *>>(option_map);
}

}  // namespace acados

%}
