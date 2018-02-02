
%{

#include <iostream>
#include <numeric>
#include <sstream>
#include <vector>

#include "acados_c/ocp_qp.h"
#include "acados_cpp/ocp_qp.hpp"
#include "acados_cpp/ocp_qp_solution.hpp"

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/print.h"

// void acados_OcpQp_A_sequence_set(acados::OcpQp *qp, LangObject *input) {
//     for (int i = 0; i < qp->N; i++) {
//         int nbElems = numRows(input) * numColumns(input);
//         std::vector<double> tmp(nbElems);
//         std::copy_n(asDoublePointer(input), nbElems, tmp.begin());
//         qp->update(qp->A, tmp);
//     }
// }

// LangObject *acados_OcpQp_A_sequence_get(acados::OcpQp *qp) {
//     std::vector<LangObject *> list_of_matrices;
//     for (int i = 0; i < qp->N; i++) {
//         int dims[2] = {qp->numRowsA(i), qp->numColsA(i)};
//         list_of_matrices.push_back(new_matrix(dims, qp->getA(i).data()));
//     }
//     return swig::from(list_of_matrices);
// }

bool is_valid_ocp_dimensions_map(const LangObject *input) {
    if (!is_map(input))
        return false;
    int_t N = int_from(input, "N");
    LangObject *nx = from(input, "nx");
    if (!is_integer(nx) && !is_sequence(nx, N+1)) {
        return false;
    }
    LangObject *nu = from(input, "nu");
    if (!is_integer(nu) && !is_sequence(nu, N)) {
        return false;
    }
    return true;
}

bool qp_dimensions_equal(const ocp_qp_dims *qp1, const ocp_qp_dims *qp2) {
    if (qp1->N != qp2->N)
        return false;
    int_t N = qp1->N;
    for (int_t i = 0; i < N; i++) {
        if (qp1->nx[i] != qp2->nx[i])
            return false;
        else if (qp1->nu[i] != qp2->nu[i])
            return false;
        else if (qp1->nb[i] != qp2->nb[i])
            return false;
        else if (qp1->ng[i] != qp2->ng[i])
            return false;
    }
    if (qp1->nx[N] != qp2->nx[N])
        return false;
    else if (qp1->nb[N] != qp2->nb[N])
        return false;
    else if (qp1->ng[N] != qp2->ng[N])
        return false;
    return true;
}

ocp_qp_dims *map_to_ocp_qp_dims(const LangObject *map) {
    if (!is_valid_ocp_dimensions_map(map)) {
        std::string err_msg =
            std::string("Input must be a valid OCP %s that specifies at least N, nx, nu")
            + std::string(LANG_MAP_NAME);
        throw std::invalid_argument(err_msg);
    }

    int_t N = int_from(map, "N");
    ocp_qp_dims *qp_dims = create_ocp_qp_dims(N);
    fill_array_from(map, "nx", qp_dims->nx, N+1);
    fill_array_from(map, "nu", qp_dims->nu, N+1);
    fill_array_from(map, "nb", qp_dims->nb, N+1);
    fill_array_from(map, "nc", qp_dims->ng, N+1);
    qp_dims->nu[N] = 0;
    // Default behavior is that initial state is fixed
    if (!has(map, "nb")) {
        qp_dims->nb[0] = qp_dims->nx[0];
    }
    return qp_dims;
}

LangObject *ocp_qp_output(const ocp_qp_in *in, const ocp_qp_out *out) {
    real_t **states_copy, **controls_copy;
    ocp_qp_dims *dims = in->dim;
    states_copy = (real_t **) malloc((dims->N+1) * sizeof(real_t *));
    controls_copy = (real_t **) malloc((dims->N+1) * sizeof(real_t *));
    for (int_t i = 0; i <= dims->N; i++) {
        states_copy[i] = (real_t *) calloc(dims->nx[i], sizeof(real_t));
        for (int_t j = 0; j < dims->nx[i]; j++)
            states_copy[i][j] = DVECEL_LIBSTR(out->ux, dims->nu[i] + j);
        controls_copy[i] = (real_t *) calloc(dims->nu[i], sizeof(real_t));
        for (int_t j = 0; j < dims->nu[i]; j++)
            controls_copy[i][j] = DVECEL_LIBSTR(out->ux, j);
    }

    LangObject *x_star = new_sequence_from(states_copy, dims->N+1, dims->nx);
    LangObject *u_star = new_sequence_from(controls_copy, dims->N+1, dims->nu);
    return new_ocp_output_tuple(x_star, u_star);
}

%}

%include "acados_cpp/options.hpp"

%rename("$ignore", %$isconstructor) ocp_qp_solution;
%include "acados_cpp/ocp_qp_solution.hpp"

%ignore make_dimensions_ptr;
%ignore extract;
%rename("$ignore", %$isconstructor) ocp_qp;
%include "acados_cpp/ocp_qp.hpp"

%rename("%s", %$isconstructor) ocp_qp;
%rename("%s") extract;

%{
using std::vector;
using std::string;
%}

%extend acados::ocp_qp {

    ocp_qp(uint N = 10, uint nx = 2, uint nu = 1, uint nbx = 0, uint nbu = 0, uint ng = 0, bool fix_x0 = true) {
        if (fix_x0 == false)
            return new acados::ocp_qp(N, nx, nu, nbx, nbu, ng);
        vector<uint> nbx_v(N+1, 0);
        nbx_v.at(0) = nx;
        acados::ocp_qp *qp = new acados::ocp_qp(vector<uint>(N+1, nx), vector<uint>(N+1, nu), nbx_v,
                                  vector<uint>(N+1, nbu), vector<uint>(N+1, ng));
        std::vector<uint> idx(nx);
        std::iota(std::begin(idx), std::end(idx), 0);
        qp->state_bounds_indices(0, idx);
        return qp;
    }

    LangObject *extract(string field) {
        vector<vector<double>> tmp = $self->extract(field);
        vector<LangObject *> result;
        for (int i = 0; i < tmp.size(); ++i)
            result.push_back(new_matrix($self->dimensions(field, i), tmp.at(i).data()));
        return swig::from(result);
    }

    vector<string> fields() {
        return vector<string>({"Q", "S", "R", "q", "r", "A", "B", "b", "lbx", "ubx", "lbu", "ubu", "C", "D", "lg", "ug"});
    }

    char *__str__() {
        static char tmp[10000];
        std::ostringstream stream;
        stream << *($self);
        std::string a = stream.str();
        std::copy(std::begin(a), std::end(a), tmp);
        return tmp;
    }
}

// %include "hpipm/include/hpipm_d_ocp_qp_sol.h"

// %include "acados/ocp_qp/ocp_qp_common.h"

%ignore ocp_qp_solve;
%ignore ocp_qp_solver;
%include "acados_c/ocp_qp.h"
