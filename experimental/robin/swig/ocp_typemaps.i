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

%include "std_map.i";

%typemap(in) std::map<std::string, acados::option_t *> {
    std::map<std::string, acados::option_t *> tmp;
#if defined(SWIGMATLAB)
    for (int i = 0; i < num_elems($input); ++i) {
        const char *key = mxGetFieldNameByNumber($input, i);
        mxArray *value = mxGetField($input, 0, key);
        tmp[std::string(key)] = acados::as_option_ptr(value);
    }
#elif defined(SWIGPYTHON)
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    while (PyDict_Next($input, &pos, &key, &value)) {
        if (!PyUnicode_Check(key))
            throw std::invalid_argument("Key must be a string");
        tmp[std::string(PyUnicode_AsUTF8AndSize(key, NULL))] = acados::as_option_ptr(value);
    }
#endif
    $1 = tmp;
}

%typemap(typecheck) std::map<std::string, acados::option_t *> {
    $1 = is_map($input) ? 1 : 0;
}

%typemap(out) std::map< std::string, std::vector<int> > {
    LangObject *result_map;

    std::vector<std::string> names;
    for (auto elem : $1) {
        names.push_back(elem.first);
    }
    const char **fieldnames = (const char **) calloc(names.size(), sizeof(char *));
    for (int i = 0; i < names.size(); ++i) {
        fieldnames[i] = names[i].c_str();
    }
#if defined(SWIGMATLAB)
    result_map = mxCreateStructMatrix(1, 1, names.size(), fieldnames);
    for (int i = 0; i < names.size(); ++i) {
        auto v = $1[names[i]];
        mxSetField(result_map, 0, fieldnames[i], new_matrix(std::make_pair(v.size(), 1), v.data()));
    }
#elif defined(SWIGPYTHON)
    result_map = PyDict_New();
    for (int i = 0; i < names.size(); ++i) {
        auto v = $1[names[i]];
        PyDict_SetItemString(result_map,
                             fieldnames[i],
                             new_matrix(std::make_pair(v.size(), 1), v.data()));
    }
#endif
    $result = result_map;
}

%typemap(in) std::vector<double> {
    if (!is_matrix($input)) {
        SWIG_exception(SWIG_ValueError, "Input is not of valid matrix type.");
        SWIG_fail;
    }
    int nbRows = numRows($input), nbCols = numColumns($input);
    int nbElems = nbRows * nbCols;
    std::vector<double> tmp(nbElems);
    std::copy_n(asDoublePointer($input), nbElems, tmp.begin());
    $1 = tmp;
}

%typemap(typecheck, precedence = SWIG_TYPECHECK_POINTER) std::vector<double> {
#if defined(SWIGMATLAB)
    $1 = (mxIsNumeric($input) ? 1 : 0);
#elif defined(SWIGPYTHON)
    $1 = (PyArray_Check($input) ? 1 : 0);
#endif
}

%typemap(out) std::vector< std::vector<double> > {
    std::vector<LangObject *> tmp;
    for (int i = 0; i < $1.size(); ++i)
        tmp.push_back(new_matrix(std::make_pair($1.at(i).size(), 1), $1.at(i).data()));
    $result = swig::from(tmp);
}


%typemap(out) ocp_qp_info  {
    const char *fields[5] = {"num_iter", "qp_solver_time", "condensing_time", "interface_time",
                             "total_time"};
#if defined(SWIGMATLAB)
    mxArray *mat_struct = mxCreateStructMatrix(1, 1, 5, fields);
    mxSetField(mat_struct, 0, fields[0], mxCreateDoubleScalar($1.num_iter));
    mxSetField(mat_struct, 0, fields[1], mxCreateDoubleScalar($1.solve_QP_time));
    mxSetField(mat_struct, 0, fields[2], mxCreateDoubleScalar($1.condensing_time));
    mxSetField(mat_struct, 0, fields[3], mxCreateDoubleScalar($1.interface_time));
    mxSetField(mat_struct, 0, fields[4], mxCreateDoubleScalar($1.total_time));
    $result = mat_struct;
#elif defined(SWIGPYTHON)
    PyObject *dict = PyDict_New();
    PyDict_SetItemString(dict, fields[0], PyLong_FromLong($1.num_iter));
    PyDict_SetItemString(dict, fields[1], PyFloat_FromDouble($1.solve_QP_time));
    PyDict_SetItemString(dict, fields[2], PyFloat_FromDouble($1.condensing_time));
    PyDict_SetItemString(dict, fields[3], PyFloat_FromDouble($1.interface_time));
    PyDict_SetItemString(dict, fields[4], PyFloat_FromDouble($1.total_time));
    $result = dict;
#endif
}

%typemap(in) int_t N {
    $1 = ($1_ltype) arg1->$1_name;
    SWIG_Error(SWIG_ValueError, "It's not allowed to change number of stages");
    SWIG_fail;
}

%typemap(out) int_t N {
    $result = to_scalar($1);
}

%typemap(in) const int_t * nx {
    SWIG_Error(SWIG_ValueError, "It's not allowed to change dimension of state vector");
    SWIG_fail;
}

%typemap(out) const int_t * nx {
    $result = new_sequence_from($1, arg1->N+1);
}

%typemap(in) const int_t * nu {
    SWIG_Error(SWIG_ValueError, "It's not allowed to change dimension of vector of controls");
    SWIG_fail;
}

%typemap(out) const int_t * nu {
    $result = new_sequence_from($1, arg1->N);
}

%typemap(in) const int_t * nb {
    SWIG_Error(SWIG_ValueError, "It's not allowed to change number of bounds");
    SWIG_fail;
}

%typemap(out) const int_t * nb {
    $result = new_sequence_from($1, arg1->N+1);
}

%typemap(in) const int_t * nc {
    SWIG_Error(SWIG_ValueError, "It's not allowed to change number of polytopic constraints");
    SWIG_fail;
}

%typemap(out) const int_t * nc {
    $result = new_sequence_from($1, arg1->N+1);
}

%typemap(in) const real_t ** B {
    $1 = ($1_ltype) arg1->$1_name;
    fill_array_from($input, $1, arg1->N, &(arg1->nx[1]), &(arg1->nu[0]));
}

%typemap(out) const real_t ** B {
    const int_t *nb_rows = &arg1->nx[1];
    const int_t *nb_cols = &arg1->nu[0];
    $result = new_sequence_from<$1_basetype>(($1_type) $1, arg1->N, nb_rows, nb_cols);
}

%typemap(in) const real_t ** b {
    $1 = ($1_ltype) arg1->$1_name;
    fill_array_from($input, $1, arg1->N, &(arg1->nx[1]));
}

%typemap(out) const real_t ** b {
    const int_t *nb_elems = &(arg1->nx[1]);
    $result = new_sequence_from<$1_basetype>(($1_type) $1, arg1->N, nb_elems);
}

%typemap(in) const real_t ** Q {
    $1 = ($1_ltype) arg1->$1_name;
    fill_array_from($input, $1, arg1->N+1, arg1->nx, arg1->nx);
}

%typemap(out) const real_t ** Q {
    const int_t *nb_rows = arg1->nx;
    const int_t *nb_cols = arg1->nx;
    $result = new_sequence_from<$1_basetype>(($1_type) $1, arg1->N+1, nb_rows, nb_cols);
}

%typemap(in) const real_t ** S {
    $1 = ($1_ltype) arg1->$1_name;
    fill_array_from($input, $1, arg1->N, arg1->nu, arg1->nx);
}

%typemap(out) const real_t ** S {
    const int_t *nb_rows = arg1->nu;
    const int_t *nb_cols = arg1->nx;
    $result = new_sequence_from<$1_basetype>(($1_type) $1, arg1->N, nb_rows, nb_cols);
}

%typemap(in) const real_t ** R {
    $1 = ($1_ltype) arg1->$1_name;
    fill_array_from($input, $1, arg1->N, arg1->nu, arg1->nu);
}

%typemap(out) const real_t ** R {
    const int_t *nb_rows = arg1->nu;
    const int_t *nb_cols = arg1->nu;
    $result = new_sequence_from<$1_basetype>(($1_type) $1, arg1->N, nb_rows, nb_cols);
}

%typemap(in) const real_t ** q {
    $1 = ($1_ltype) arg1->$1_name;
    fill_array_from($input, $1, arg1->N+1, arg1->nx);
}

%typemap(out) const real_t ** q {
    const int_t *nb_elems = arg1->nx;
    $result = new_sequence_from<$1_basetype>(($1_type) $1, arg1->N+1, nb_elems);
}

%typemap(in) const real_t ** r {
    $1 = ($1_ltype) arg1->$1_name;
    fill_array_from($input, $1, arg1->N, arg1->nu);
}

%typemap(out) const real_t ** r {
    const int_t *nb_elems = arg1->nu;
    $result = new_sequence_from<$1_basetype>(($1_type) $1, arg1->N, nb_elems);
}

%typemap(in) const int_t ** idxb {
    $1 = ($1_ltype) arg1->$1_name;
    fill_array_from($input, $1, arg1->N+1, arg1->nb);
}

%typemap(out) const int_t ** idxb {
    const int_t *nb_elems = arg1->nb;
    $result = new_sequence_from<$1_basetype>(($1_type) $1, arg1->N+1, nb_elems);
}

%typemap(in) const real_t ** lb {
    $1 = ($1_ltype) arg1->$1_name;
    fill_array_from($input, $1, arg1->N+1, arg1->nb);
}

%typemap(out) const real_t ** lb {
    const int_t *nb_elems = arg1->nb;
    $result = new_sequence_from<$1_basetype>(($1_type) $1, arg1->N+1, nb_elems);
}

%typemap(in) const real_t ** ub {
    $1 = ($1_ltype) arg1->$1_name;
    fill_array_from($input, $1, arg1->N+1, arg1->nb);
}

%typemap(out) const real_t ** ub {
    const int_t *nb_elems = arg1->nb;
    $result = new_sequence_from<$1_basetype>(($1_type) $1, arg1->N+1, nb_elems);
}

%typemap(in) const real_t ** Cx {
    $1 = ($1_ltype) arg1->$1_name;
    fill_array_from($input, $1, arg1->N+1, arg1->nc, arg1->nx);
}

%typemap(out) const real_t ** Cx {
    const int_t *nb_rows = arg1->nc;
    const int_t *nb_cols = arg1->nx;
    $result = new_sequence_from<$1_basetype>(($1_type) $1, arg1->N+1, nb_rows, nb_cols);
}

%typemap(in) const real_t ** Cu {
    $1 = ($1_ltype) arg1->$1_name;
    fill_array_from($input, $1, arg1->N, arg1->nc, arg1->nu);
}

%typemap(out) const real_t ** Cu {
    const int_t *nb_rows = arg1->nc;
    const int_t *nb_cols = arg1->nu;
    $result = new_sequence_from<$1_basetype>(($1_type) $1, arg1->N, nb_rows, nb_cols);
}

%typemap(in) const real_t ** lc {
    $1 = ($1_ltype) arg1->$1_name;
    fill_array_from($input, $1, arg1->N+1, arg1->nc);
}

%typemap(out) const real_t ** lc {
    const int_t *nb_elems = arg1->nc;
    $result = new_sequence_from<$1_basetype>(($1_type) $1, arg1->N+1, nb_elems);
}

%typemap(in) const real_t ** uc {
    $1 = ($1_ltype) arg1->$1_name;
    fill_array_from($input, $1, arg1->N+1, arg1->nc);
}

%typemap(out) const real_t ** uc {
    const int_t *nb_elems = arg1->nc;
    $result = new_sequence_from<$1_basetype>(($1_type) $1, arg1->N+1, nb_elems);
}

%typemap(in) const real_t ** lg {
    $1 = ($1_ltype) arg1->$1_name;
    fill_array_from($input, $1, arg1->N+1, arg1->ng);
}

%typemap(out) const real_t ** lg {
    const int_t *nb_elems = arg1->ng;
    $result = new_sequence_from<$1_basetype>(($1_type) $1, arg1->N+1, nb_elems);
}

%typemap(in) const real_t ** ug {
    $1 = ($1_ltype) arg1->$1_name;
    fill_array_from($input, $1, arg1->N+1, arg1->ng);
}

%typemap(out) const real_t ** ug {
    const int_t *nb_elems = arg1->ng;
    $result = new_sequence_from<$1_basetype>(($1_type) $1, arg1->N+1, nb_elems);
}
