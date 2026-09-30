#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.


cdef extern from "acados/ocp_nlp/ocp_nlp_common.h":
    ctypedef struct ocp_nlp_config:
        pass

    ctypedef struct ocp_nlp_dims:
        pass

    ctypedef struct ocp_nlp_in:
        pass

    ctypedef struct ocp_nlp_out:
        pass


cdef extern from "acados_c/ocp_nlp_interface.h":
    ctypedef enum ocp_nlp_solver_t:
        pass

    ctypedef enum ocp_nlp_cost_t:
        pass

    ctypedef enum ocp_nlp_dynamics_t:
        pass

    ctypedef enum ocp_nlp_constraints_t:
        pass

    ctypedef enum ocp_nlp_reg_t:
        pass

    ctypedef struct ocp_nlp_plan:
        pass

    ctypedef struct ocp_nlp_solver:
        pass

    int ocp_nlp_cost_model_set(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in_,
        int start_stage, const char *field, void *value)
    int ocp_nlp_constraints_model_set(ocp_nlp_config *config, ocp_nlp_dims *dims,
        ocp_nlp_in *in_, ocp_nlp_out *out_, int stage, const char *field, void *value)

    # out
    void ocp_nlp_out_set(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_out *out, ocp_nlp_in *in_,
        int stage, const char *field, void *value)
    void ocp_nlp_out_get(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_out *out,
        int stage, const char *field, void *value)
    void ocp_nlp_get_at_stage(ocp_nlp_solver *solver, int stage, const char *field, void *value)
    void ocp_nlp_get_from_iterate(ocp_nlp_solver *solver, int iter, int stage, const char *field, void *value)
    int ocp_nlp_dims_get_from_attr(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_out *out,
        int stage, const char *field)
    void ocp_nlp_constraint_dims_get_from_attr(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_out *out,
        int stage, const char *field, int *dims_out)
    void ocp_nlp_cost_dims_get_from_attr(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_out *out,
        int stage, const char *field, int *dims_out)
    void ocp_nlp_qp_dims_get_from_attr(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_out *out,
        int stage, const char *field, int *dims_out)

    # in
    void ocp_nlp_in_set(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *nlp_in,
        int stage, const char *field, void *value)
    void ocp_nlp_in_get(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *nlp_in,
        int stage, const char *field, void *value)

    # opts
    void ocp_nlp_solver_opts_set(ocp_nlp_config *config, void *opts_, const char *field, void* value)

    # solver
    void ocp_nlp_eval_residuals(ocp_nlp_solver *solver, ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out)
    void ocp_nlp_eval_param_sens(ocp_nlp_solver *solver, char *field, int stage, int index, ocp_nlp_out *sens_nlp_out)
    void ocp_nlp_eval_cost(ocp_nlp_solver *solver, ocp_nlp_in *nlp_in_, ocp_nlp_out *nlp_out)
    void ocp_nlp_eval_params_jac(ocp_nlp_solver *solver, ocp_nlp_in *nlp_in_, ocp_nlp_out *nlp_out)
    void ocp_nlp_eval_lagrange_grad_p(ocp_nlp_solver *solver, ocp_nlp_in *nlp_in_, const char *field, void* value)
    # get/set
    void ocp_nlp_get(ocp_nlp_solver *solver, const char *field, void *return_value_)
    void ocp_nlp_set(ocp_nlp_solver *solver, int stage, const char *field, void *value)
