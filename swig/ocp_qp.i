
%{

#include <iostream>
#include <vector>

#include "acados_cpp/ocp_qp.hpp"

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados_c/ocp_qp.h"

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

%ignore extract;
%include "acados_cpp/ocp_qp.hpp"

%rename("%s") extract;

%{
using std::vector;
using std::string;
%}

%extend acados::ocp_qp {

    LangObject *extract(string field) {
        vector<vector<double>> tmp = $self->extract(field);
        vector<LangObject *> result;
        int dims_array[2] = {2, 3};
        double data[6] = {0, 0, 0, 0, 0, 0};
        result.push_back(new_matrix(dims_array, data));
        result.push_back(new_matrix(dims_array, data));
        result.push_back(new_matrix(dims_array, data));
        result.push_back(new_matrix(dims_array, data));
        LangObject *obj = swig::from(result);
        Py_INCREF(obj);
        return obj;
    }

    vector<string> fields() {
        return vector<string>({"Q", "S", "R", "q", "r", "A", "B", "b", "lbx", "ubx", "lbu", "ubu", "C", "D", "lg", "ug"});
    }

    char *__str__() {
        static char tmp[1];
        std::cout << *($self);
        return tmp;
    }
}

%include "hpipm/include/hpipm_d_ocp_qp_sol.h"

%include "acados/ocp_qp/ocp_qp_common.h"

// %include "acados_c/ocp_qp.h"

// %extend ocp_qp_solver {

//     ocp_qp_solver(ocp_qp_solver_t solver_name, const acados::OcpQp& qp, LangObject *options = NONE) {

//         ocp_qp_solver_plan plan;
//         plan.qp_solver = solver_name;
        
//         void *args = ocp_qp_create_args(&plan, qp.dimensions);
//         ocp_qp_solver *solver = ocp_qp_create(&plan, qp.dimensions, args);
//         return solver;
//     }

//     d_ocp_qp_sol *evaluate(const acados::OcpQp& input) {
//         d_ocp_qp_sol *result = create_ocp_qp_out(input.dimensions);
//         int_t return_code = ocp_qp_solve($self, input.qp, result);
//         if (return_code != 0)
//             throw std::runtime_error("qp solver failed!");
//         return result;
//     }

// }
