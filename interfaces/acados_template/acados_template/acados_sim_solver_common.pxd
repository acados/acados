#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.


cdef extern from "acados/sim/sim_common.h":
    ctypedef struct sim_config:
        pass

    ctypedef struct sim_opts:
        pass

    ctypedef struct sim_in:
        pass

    ctypedef struct sim_out:
        pass


cdef extern from "acados_c/sim_interface.h":

    ctypedef struct sim_plan:
        pass

    ctypedef struct sim_solver:
        pass

    # out
    void sim_out_get(sim_config *config, void *dims, sim_out *out, const char *field, void *value)
    int sim_dims_get_from_attr(sim_config *config, void *dims, const char *field, void *dims_data)

    # mem
    void sim_memory_get(sim_config *config, void *dims, void *mem, const char *field, void *value)

    # opts
    void sim_opts_set(sim_config *config, void *opts_, const char *field, void *value)

    # get/set
    void sim_in_set(sim_config *config, void *dims, sim_in *sim_in, const char *field, void *value)
    void sim_solver_set(sim_solver *solver, const char *field, void *value)