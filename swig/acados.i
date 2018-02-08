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

%module acados

#if defined(SWIGMATLAB)
typedef mxArray LangObject;
%{
typedef mxArray LangObject;
%}
#define NONE NULL
#elif defined(SWIGPYTHON)
typedef PyObject LangObject;
%{
typedef PyObject LangObject;
%}
#define NONE Py_None
%ignore SwigPyIterator;
%ignore SwigPyIterator_swigregister;
#endif

%include "std_string.i"

%include "std_vector.i"
namespace std {
    %template(vector_i) vector<int>;
    %template(vector_string) vector<string>;
    %template(vector_O) vector<LangObject *>;
};

%include "std_pair.i"
namespace std {
    %template(pair_ii) pair<int, int>;
}

#if defined(SWIGPYTHON)
%{
#define SWIG_FILE_WITH_INIT
%}
#endif

%include "exception.i"

%exception {
    try {
        $action
    } catch (const std::exception& e) {
        SWIG_exception(SWIG_RuntimeError, e.what());
    }
}

%include "conversions.i"
%include "ocp_typemaps.i"

%{
// TODO(dimitris): support compilation with visual studio
#if (defined _WIN32 || defined _WIN64 || defined __MINGW32__ || defined __MINGW64__)
#include <windows.h>
char compiler[16] = "gcc";
#else
#include <dlfcn.h>
char compiler[16] = "cc";
#endif
// #include <xmmintrin.h>  // for floating point exceptions
// _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
%}

%feature("autodoc", "3");
%include "ocp_qp.i"
