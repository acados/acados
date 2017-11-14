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

%typemap(in) const real_t ** A {
    $1 = ($1_ltype) arg1->$1_name;
    fill_array_from($input, $1, arg1->N, &(arg1->nx[1]), &(arg1->nx[0]));
}

%typemap(out) const real_t ** A {
    const int_t *nb_rows = &arg1->nx[1];
    const int_t *nb_cols = &arg1->nx[0];
    $result = new_sequence_from<$1_basetype>(($1_type) $1, arg1->N, nb_rows, nb_cols);
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
