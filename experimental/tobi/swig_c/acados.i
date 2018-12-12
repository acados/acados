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

%include "cpointer.i"
%pointer_class(int, intp)
%pointer_class(bool, boolp)

%{
#include "acados/utils/types.h"
#include "acados/utils/external_function_generic.h"
#include "acados_c/external_function_interface.h"
#include "acados_c/sim_interface.h"

//#include "acados_c/ocp_nlp_interface.h"
%}
%include "acados/utils/types.h"
%include "acados/utils/external_function_generic.h"
%include "acados_c/external_function_interface.h"
%include "acados_c/sim_interface.h"

//%include "acados_c/ocp_nlp_interface.h"


#if defined(SWIGPYTHON)
%{
#define SWIG_FILE_WITH_INIT
%}
#endif

//%include "conversions.i"

%feature("autodoc", "3");
