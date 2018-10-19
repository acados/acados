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
#endif

%include "acados/utils/types.h"

%include "std_string.i"

%include "std_vector.i"
namespace std {
    %template(vector_i) vector<int>;
    %template(vector_d) vector<double>;
    %template(vector_string) vector<string>;
    %template(vector_O) vector<LangObject *>;
};

%include "std_pair.i"
namespace std {
    %template(pair_ii) pair<int, int>;
}

%include "std_map.i"
namespace std {
    %template(map_si) map< string, vector<int> >;
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

%feature("autodoc", "3");
%include "ocp_qp.i"
%include "ocp_nlp.i"

%include "integrator.i"
