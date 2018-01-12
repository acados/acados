
%{

#include <iostream>
#include <vector>

#include "acados_cpp/ocp_qp.h"

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados_c/ocp_qp.h"

void OcpQp_A_set(OcpQp *qp, LangObject *input) {
    for (int i = 0; i < qp->N; i++)
        qp->setA(i, asDoublePointer(input), numRows(input), numColumns(input));
}

LangObject *OcpQp_A_get(OcpQp *qp) {
    std::vector<LangObject *> list_of_matrices;
    for (int i = 0; i < qp->N; i++) {
        int dims[2] = {qp->numRowsA(i), qp->numColsA(i)};
        list_of_matrices.push_back(new_matrix(dims, qp->getA(i).data()));
    }
    return swig::from(list_of_matrices);
}

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

%rename("%s") OcpQp;
%include "acados_cpp/ocp_qp.h"

%extend OcpQp {

    LangObject *A;

    void setQ(LangObject *input) {
        for (int i = 0; i <= $self->N; i++)
            $self->setQ(i, asDoublePointer(input));
    }

    void setS(LangObject *input) {
        for (int i = 0; i <= $self->N; i++)
            $self->setS(i, asDoublePointer(input));
    }

    void setR(LangObject *input) {
        for (int i = 0; i <= $self->N; i++)
            $self->setR(i, asDoublePointer(input));
    }

    void setq(LangObject *input) {
        for (int i = 0; i <= $self->N; i++)
            $self->setq(i, asDoublePointer(input));
    }

    void setr(LangObject *input) {
        for (int i = 0; i <= $self->N; i++)
            $self->setr(i, asDoublePointer(input));
    }

    void setA(LangObject *input) {
        for (int i = 0; i < $self->N; i++)
            $self->setA(i, asDoublePointer(input));
    }

    void setB(LangObject *input) {
        for (int i = 0; i < $self->N; i++)
            $self->setB(i, asDoublePointer(input));
    }

    void setb(LangObject *input) {
        for (int i = 0; i < $self->N; i++)
            $self->setb(i, asDoublePointer(input));
    }

    void setlb(LangObject *input) {
        for (int i = 0; i <= $self->N; i++)
            $self->setlb(i, asDoublePointer(input));
    }

    void setub(LangObject *input) {
        for (int i = 0; i <= $self->N; i++)
            $self->setub(i, asDoublePointer(input));
    }

    void setC(LangObject *input) {
        for (int i = 0; i <= $self->N; i++)
            $self->setC(i, asDoublePointer(input));
    }

    void setD(LangObject *input) {
        for (int i = 0; i <= $self->N; i++)
            $self->setD(i, asDoublePointer(input));
    }

    void setlg(LangObject *input) {
        for (int i = 0; i <= $self->N; i++)
            $self->setlg(i, asDoublePointer(input));
    }

    void setug(LangObject *input) {
        for (int i = 0; i <= $self->N; i++)
            $self->setug(i, asDoublePointer(input));
    }

    char *__str__() {
        static char tmp[1];
        std::cout << *($self);
        return tmp;
    }
}

%rename("%s") ocp_qp_solver;
%include "acados_c/ocp_qp.h"

%extend ocp_qp_solver {

    ocp_qp_solver(ocp_qp_solver_t solver_name, LangObject *dimensions, LangObject *options = NONE) {

        ocp_qp_dims *dims = map_to_ocp_qp_dims(dimensions);

        ocp_qp_solver_plan plan;
        plan.qp_solver = solver_name;
        
        void *args = ocp_qp_create_args(&plan, dims);
        ocp_qp_solver *solver = ocp_qp_create(&plan, dims, args);
        return solver;
    }

    LangObject *evaluate(ocp_qp_in *input) {
        ocp_qp_out *result = create_ocp_qp_out(input->dim);
        int_t return_code = ocp_qp_solve($self, input, result);
        if (return_code != 0)
            throw std::runtime_error("qp solver failed!");
        return ocp_qp_output(input, result);
    }

}