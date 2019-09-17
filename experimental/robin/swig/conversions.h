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

#ifndef SWIG_CONVERSIONS_H_
#define SWIG_CONVERSIONS_H_

#if defined(SWIGMATLAB)

typedef mxArray LangObject;
#define LANG_SEQUENCE_NAME "cell array"
#define LANG_MAP_NAME "struct"
#define LANG_MATRIX_NAME "matrix"

#elif defined(SWIGPYTHON)

#define NO_IMPORT_ARRAY
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"

typedef PyObject LangObject;
#define LANG_SEQUENCE_NAME "list"
#define LANG_MAP_NAME "dictionary"
#define LANG_MATRIX_NAME "ndarray"

#define PyArray_NewFromDataF(nd, dims, data)                                  \
    PyArray_New(&PyArray_Type, nd, dims, NPY_FLOAT64, NULL, (void *) data, 0, \
                NPY_ARRAY_FARRAY | NPY_ARRAY_OWNDATA, NULL)
#endif

#include "acados/utils/types.h"

#endif  // SWIG_CONVERSIONS_H_
