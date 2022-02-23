# -*- coding: future_fstrings -*-
#
# Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
# Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
# Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
# Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
#
# This file is part of acados.
#
# The 2-Clause BSD License
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.;
#
# cython: language_level=3
# cython: profile=False
# distutils: language=c

cimport cython
from libc cimport string

cimport acados_solver_common
cimport acados_solver

cimport numpy as cnp

import os
import json
from datetime import datetime
import numpy as np


cdef class AcadosOcpSolverFast:
    """
    Class to interact with the acados ocp solver C object.

        :param model_name:
        :param N:
    """

    cdef acados_solver.nlp_solver_capsule *capsule
    cdef void *nlp_opts
    cdef acados_solver_common.ocp_nlp_dims *nlp_dims
    cdef acados_solver_common.ocp_nlp_config *nlp_config
    cdef acados_solver_common.ocp_nlp_out *nlp_out
    cdef acados_solver_common.ocp_nlp_in *nlp_in
    cdef acados_solver_common.ocp_nlp_solver *nlp_solver

    cdef str model_name
    cdef int N
    cdef bint solver_created

    def __cinit__(self, str model_name, int N):
        self.model_name = model_name
        self.N = N

        self.solver_created = False

        # create capsule
        self.capsule = acados_solver.acados_create_capsule()

        # create solver
        assert acados_solver.acados_create(self.capsule) == 0
        self.solver_created = True

        # get pointers solver
        self.nlp_opts = acados_solver.acados_get_nlp_opts(self.capsule)
        self.nlp_dims = acados_solver.acados_get_nlp_dims(self.capsule)
        self.nlp_config = acados_solver.acados_get_nlp_config(self.capsule)
        self.nlp_out = acados_solver.acados_get_nlp_out(self.capsule)
        self.nlp_in = acados_solver.acados_get_nlp_in(self.capsule)
        self.nlp_solver = acados_solver.acados_get_nlp_solver(self.capsule)


    def solve(self):
        """
        Solve the ocp with current input.
        """
        return acados_solver.acados_solve(self.capsule)


    def set_new_time_steps(self, new_time_steps):
        """
        Set new time steps before solving. Only reload library without code generation but with new time steps.

            :param new_time_steps: vector of new time steps for the solver

            .. note:: This allows for different use-cases: either set a new size of time-steps or a new distribution of
                      the shooting nodes without changing the number, e.g., to reach a different final time. Both cases
                      do not require a new code export and compilation.
        """
        raise NotImplementedError()


    def get(self, int stage, str field_):
        """
        Get the last solution of the solver:

            :param stage: integer corresponding to shooting node
            :param field: string in ['x', 'u', 'z', 'pi', 'lam', 't', 'sl', 'su',]

            .. note:: regarding lam, t: \n
                    the inequalities are internally organized in the following order: \n
                    [ lbu lbx lg lh lphi ubu ubx ug uh uphi; \n
                      lsbu lsbx lsg lsh lsphi usbu usbx usg ush usphi]

            .. note:: pi: multipliers for dynamics equality constraints \n
                      lam: multipliers for inequalities \n
                      t: slack variables corresponding to evaluation of all inequalities (at the solution) \n
                      sl: slack variables of soft lower inequality constraints \n
                      su: slack variables of soft upper inequality constraints \n
        """

        out_fields = ['x', 'u', 'z', 'pi', 'lam', 't', 'sl', 'su']
        field = field_.encode('utf-8')

        if field_ not in out_fields:
            raise Exception('AcadosOcpSolver.get(): {} is an invalid argument.\
                    \n Possible values are {}. Exiting.'.format(field_, out_fields))

        if stage < 0 or stage > self.N:
            raise Exception('AcadosOcpSolver.get(): stage index must be in [0, N], got: {}.'.format(self.N))

        if stage == self.N and field_ == 'pi':
            raise Exception('AcadosOcpSolver.get(): field {} does not exist at final stage {}.'\
                .format(field_, stage))

        cdef int dims = acados_solver_common.ocp_nlp_dims_get_from_attr(self.nlp_config,
            self.nlp_dims, self.nlp_out, stage, field)

        cdef cnp.ndarray[cnp.float64_t, ndim=1] out = np.zeros((dims,))
        acados_solver_common.ocp_nlp_out_get(self.nlp_config, \
            self.nlp_dims, self.nlp_out, stage, field, <void *> out.data)

        return out


    def print_statistics(self):
        """
        prints statistics of previous solver run as a table:
            - iter: iteration number
            - res_stat: stationarity residual
            - res_eq: residual wrt equality constraints (dynamics)
            - res_ineq: residual wrt inequality constraints (constraints)
            - res_comp: residual wrt complementarity conditions
            - qp_stat: status of QP solver
            - qp_iter: number of QP iterations
            - qp_res_stat: stationarity residual of the last QP solution
            - qp_res_eq: residual wrt equality constraints (dynamics) of the last QP solution
            - qp_res_ineq: residual wrt inequality constraints (constraints)  of the last QP solution
            - qp_res_comp: residual wrt complementarity conditions of the last QP solution
        """
        acados_solver.acados_print_stats(self.capsule)


    def store_iterate(self, filename='', overwrite=False):
        """
        Stores the current iterate of the ocp solver in a json file.

            :param filename: if not set, use model_name + timestamp + '.json'
            :param overwrite: if false and filename exists add timestamp to filename
        """
        if filename == '':
            filename += self.model_name + '_' + 'iterate' + '.json'

        if not overwrite:
            # append timestamp
            if os.path.isfile(filename):
                filename = filename[:-5]
                filename += datetime.utcnow().strftime('%Y-%m-%d-%H:%M:%S.%f') + '.json'

        # get iterate:
        solution = dict()

        for i in range(self.N+1):
            solution['x_'+str(i)] = self.get(i,'x')
            solution['u_'+str(i)] = self.get(i,'u')
            solution['z_'+str(i)] = self.get(i,'z')
            solution['lam_'+str(i)] = self.get(i,'lam')
            solution['t_'+str(i)] = self.get(i, 't')
            solution['sl_'+str(i)] = self.get(i, 'sl')
            solution['su_'+str(i)] = self.get(i, 'su')
        for i in range(self.N):
            solution['pi_'+str(i)] = self.get(i,'pi')

        # save
        with open(filename, 'w') as f:
            json.dump(solution, f, default=lambda x: x.tolist(), indent=4, sort_keys=True)
        print("stored current iterate in ", os.path.join(os.getcwd(), filename))


    def load_iterate(self, filename):
        """
        Loads the iterate stored in json file with filename into the ocp solver.
        """
        if not os.path.isfile(filename):
            raise Exception('load_iterate: failed, file does not exist: ' + os.path.join(os.getcwd(), filename))

        with open(filename, 'r') as f:
            solution = json.load(f)

        for key in solution.keys():
            (field, stage) = key.split('_')
            self.set(int(stage), field, np.array(solution[key]))


    def get_stats(self, field_):
        """
        Get the information of the last solver call.

            :param field: string in ['statistics', 'time_tot', 'time_lin', 'time_sim', 'time_sim_ad', 'time_sim_la', 'time_qp', 'time_qp_solver_call', 'time_reg', 'sqp_iter']
        """
        fields = ['time_tot',  # total cpu time previous call
                  'time_lin',  # cpu time for linearization
                  'time_sim',  # cpu time for integrator
                  'time_sim_ad',  # cpu time for integrator contribution of external function calls
                  'time_sim_la',  # cpu time for integrator contribution of linear algebra
                  'time_qp',   # cpu time qp solution
                  'time_qp_solver_call',  # cpu time inside qp solver (without converting the QP)
                  'time_qp_xcond',
                  'time_glob',  # cpu time globalization
                  'time_reg',  # cpu time regularization
                  'sqp_iter',  # number of SQP iterations
                  'qp_iter',  # vector of QP iterations for last SQP call
                  'statistics',  # table with info about last iteration
                  'stat_m',
                  'stat_n',]

        field = field_
        field = field.encode('utf-8')

        if (field_ not in fields):
            raise Exception('AcadosOcpSolver.get_stats(): {} is not a valid argument.\
                    \n Possible values are {}. Exiting.'.format(fields, fields))

        if field_ in ['sqp_iter', 'stat_m', 'stat_n']:
            return self.__get_stat_int(field)

        elif field_.startswith('time'):
            return self.__get_stat_double(field)

        elif field_ == 'statistics':
            sqp_iter = self.get_stats("sqp_iter")
            stat_m = self.get_stats("stat_m")
            stat_n = self.get_stats("stat_n")
            min_size = min([stat_m, sqp_iter+1])
            return self.__get_stat_matrix(field, stat_n+1, min_size)

        elif field_ == 'qp_iter':
            NotImplementedError("TODO! Cython not aware if SQP or SQP_RTI")
            full_stats = self.get_stats('statistics')
            if self.acados_ocp.solver_options.nlp_solver_type == 'SQP':
                out = full_stats[6, :]
            elif self.acados_ocp.solver_options.nlp_solver_type == 'SQP_RTI':
                out = full_stats[2, :]
        else:
            NotImplementedError("TODO!")


    def __get_stat_int(self, field):
        cdef int out
        acados_solver_common.ocp_nlp_get(self.nlp_config, self.nlp_solver, field, <void *> &out)
        return out

    def __get_stat_double(self, field):
        cdef cnp.ndarray[cnp.float64_t, ndim=1] out = np.zeros((1,))
        acados_solver_common.ocp_nlp_get(self.nlp_config, self.nlp_solver, field, <void *> out.data)
        return out

    def __get_stat_matrix(self, field, n, m):
        cdef cnp.ndarray[cnp.float64_t, ndim=2] out_mat = np.ascontiguousarray(np.zeros((n, m)), dtype=np.float64)
        acados_solver_common.ocp_nlp_get(self.nlp_config, self.nlp_solver, field, <void *> out_mat.data)
        return out_mat


    def get_cost(self):
        """
        Returns the cost value of the current solution.
        """
        # compute cost internally
        acados_solver_common.ocp_nlp_eval_cost(self.nlp_solver, self.nlp_in, self.nlp_out)

        # create output
        cdef double out

        # call getter
        acados_solver_common.ocp_nlp_get(self.nlp_config, self.nlp_solver, "cost_value", <void *> &out)

        return out


    def get_residuals(self):
        """
        Returns an array of the form [res_stat, res_eq, res_ineq, res_comp].
        """
        # TODO: check if RTI, only eval then
        # if self.acados_ocp.solver_options.nlp_solver_type == 'SQP_RTI':
        # compute residuals if RTI
        acados_solver_common.ocp_nlp_eval_residuals(self.nlp_solver, self.nlp_in, self.nlp_out)

        # create output array
        cdef cnp.ndarray[cnp.float64_t, ndim=1] out = np.ascontiguousarray(np.zeros((4,), dtype=np.float64))
        cdef double double_value

        field = "res_stat".encode('utf-8')
        acados_solver_common.ocp_nlp_get(self.nlp_config, self.nlp_solver, field, <void *> &double_value)
        out[0] = double_value

        field = "res_eq".encode('utf-8')
        acados_solver_common.ocp_nlp_get(self.nlp_config, self.nlp_solver, field, <void *> &double_value)
        out[1] = double_value

        field = "res_ineq".encode('utf-8')
        acados_solver_common.ocp_nlp_get(self.nlp_config, self.nlp_solver, field, <void *> &double_value)
        out[2] = double_value

        field = "res_comp".encode('utf-8')
        acados_solver_common.ocp_nlp_get(self.nlp_config, self.nlp_solver, field, <void *> &double_value)
        out[3] = double_value

        return out


    # Note: this function should not be used anymore, better use cost_set, constraints_set
    def set(self, int stage, str field_, value_):

        """
        Set numerical data inside the solver.

            :param stage: integer corresponding to shooting node
            :param field: string in ['x', 'u', 'pi', 'lam', 't', 'p']

            .. note:: regarding lam, t: \n
                    the inequalities are internally organized in the following order: \n
                    [ lbu lbx lg lh lphi ubu ubx ug uh uphi; \n
                      lsbu lsbx lsg lsh lsphi usbu usbx usg ush usphi]

            .. note:: pi: multipliers for dynamics equality constraints \n
                      lam: multipliers for inequalities \n
                      t: slack variables corresponding to evaluation of all inequalities (at the solution) \n
                      sl: slack variables of soft lower inequality constraints \n
                      su: slack variables of soft upper inequality constraints \n
        """
        cost_fields = ['y_ref', 'yref']
        constraints_fields = ['lbx', 'ubx', 'lbu', 'ubu']
        out_fields = ['x', 'u', 'pi', 'lam', 't', 'z', 'sl', 'su']

        field = field_.encode('utf-8')

        cdef cnp.ndarray[cnp.float64_t, ndim=1] value = np.ascontiguousarray(value_, dtype=np.float64)

        # treat parameters separately
        if field_ == 'p':
            assert acados_solver.acados_update_params(self.capsule, stage, <double *> value.data, value.shape[0]) == 0
        else:
            if field_ not in constraints_fields + cost_fields + out_fields:
                raise Exception("AcadosOcpSolver.set(): {} is not a valid argument.\
                    \nPossible values are {}. Exiting.".format(field, \
                    constraints_fields + cost_fields + out_fields + ['p']))

            dims = acados_solver_common.ocp_nlp_dims_get_from_attr(self.nlp_config,
                self.nlp_dims, self.nlp_out, stage, field)

            if value_.shape[0] != dims:
                msg = 'AcadosOcpSolver.set(): mismatching dimension for field "{}" '.format(field_)
                msg += 'with dimension {} (you have {})'.format(dims, value_.shape[0])
                raise Exception(msg)

            if field_ in constraints_fields:
                acados_solver_common.ocp_nlp_constraints_model_set(self.nlp_config,
                    self.nlp_dims, self.nlp_in, stage, field, <void *> value.data)
            elif field_ in cost_fields:
                acados_solver_common.ocp_nlp_cost_model_set(self.nlp_config,
                    self.nlp_dims, self.nlp_in, stage, field, <void *> value.data)
            elif field_ in out_fields:
                acados_solver_common.ocp_nlp_out_set(self.nlp_config,
                    self.nlp_dims, self.nlp_out, stage, field, <void *> value.data)


    def cost_set(self, int stage, str field_, value_):
        """
        Set numerical data in the cost module of the solver.

            :param stage: integer corresponding to shooting node
            :param field: string, e.g. 'yref', 'W', 'ext_cost_num_hess'
            :param value: of appropriate size
        """
        field = field_.encode('utf-8')

        cdef int dims[2]
        acados_solver_common.ocp_nlp_cost_dims_get_from_attr(self.nlp_config, \
            self.nlp_dims, self.nlp_out, stage, field, &dims[0])

        cdef double[::1,:] value

        value_shape = value_.shape
        if len(value_shape) == 1:
            value_shape = (value_shape[0], 0)
            value = np.asfortranarray(value_[None,:])

        elif len(value_shape) == 2:
            # Get elements in column major order
            value = np.asfortranarray(value_)

        if value_shape[0] != dims[0] or value_shape[1] != dims[1]:
            raise Exception('AcadosOcpSolver.cost_set(): mismatching dimension', \
                ' for field "{}" with dimension {} (you have {})'.format( \
                field_, tuple(dims), value_shape))

        acados_solver_common.ocp_nlp_cost_model_set(self.nlp_config, \
            self.nlp_dims, self.nlp_in, stage, field, <void *> &value[0][0])


    def constraints_set(self, int stage, str field_, value_):
        """
        Set numerical data in the constraint module of the solver.

            :param stage: integer corresponding to shooting node
            :param field: string in ['lbx', 'ubx', 'lbu', 'ubu', 'lg', 'ug', 'lh', 'uh', 'uphi', 'C', 'D']
            :param value: of appropriate size
        """
        field = field_.encode('utf-8')

        cdef int dims[2]
        acados_solver_common.ocp_nlp_constraint_dims_get_from_attr(self.nlp_config, \
            self.nlp_dims, self.nlp_out, stage, field, &dims[0])

        cdef double[::1,:] value

        value_shape = value_.shape
        if len(value_shape) == 1:
            value_shape = (value_shape[0], 0)
            value = np.asfortranarray(value_[None,:])

        elif len(value_shape) == 2:
            # Get elements in column major order
            value = np.asfortranarray(value_)

        if value_shape[0] != dims[0] or value_shape[1] != dims[1]:
            raise Exception('AcadosOcpSolver.constraints_set(): mismatching dimension' \
                ' for field "{}" with dimension {} (you have {})'.format(field_, tuple(dims), value_shape))

        acados_solver_common.ocp_nlp_constraints_model_set(self.nlp_config, \
            self.nlp_dims, self.nlp_in, stage, field, <void *> &value[0][0])

        return


    def dynamics_get(self, int stage, str field_):
        """
        Get numerical data from the dynamics module of the solver:

            :param stage: integer corresponding to shooting node
            :param field: string, e.g. 'A'
        """
        field = field_.encode('utf-8')

        # get dims
        cdef int[2] dims
        acados_solver_common.ocp_nlp_dynamics_dims_get_from_attr(self.nlp_config, self.nlp_dims, self.nlp_out, stage, field, &dims[0])

        # create output data
        cdef cnp.ndarray[cnp.float64_t, ndim=2] out = np.zeros((dims[0], dims[1]), order='F')

        # call getter
        acados_solver_common.ocp_nlp_get_at_stage(self.nlp_config, self.nlp_dims, self.nlp_solver, stage, field, <void *> out.data)

        return out


    def options_set(self, str field_, value_):
        """
        Set options of the solver.

            :param field: string, e.g. 'print_level', 'rti_phase', 'initialize_t_slacks', 'step_length', 'alpha_min', 'alpha_reduction'
            :param value: of type int, float
        """
        int_fields = ['print_level', 'rti_phase', 'initialize_t_slacks']
        double_fields = ['step_length', 'tol_eq', 'tol_stat', 'tol_ineq', 'tol_comp', 'alpha_min', 'alpha_reduction']
        string_fields = ['globalization']

        # encode
        field = field_.encode('utf-8')

        cdef int int_value
        cdef double double_value
        cdef unsigned char[::1] string_value

        # check field availability and type
        if field_ in int_fields:
            if not isinstance(value_, int):
                raise Exception('solver option {} must be of type int. You have {}.'.format(field_, type(value_)))

            if field_ == 'rti_phase':
                if value_ < 0 or value_ > 2:
                    raise Exception('AcadosOcpSolver.solve(): argument \'rti_phase\' can '
                        'take only values 0, 1, 2 for SQP-RTI-type solvers')
                if self.acados_ocp.solver_options.nlp_solver_type != 'SQP_RTI' and value_ > 0:
                    raise Exception('AcadosOcpSolver.solve(): argument \'rti_phase\' can '
                        'take only value 0 for SQP-type solvers')

            int_value = value_
            acados_solver_common.ocp_nlp_solver_opts_set(self.nlp_config, self.nlp_opts, field, <void *> &int_value)

        elif field_ in double_fields:
            if not isinstance(value_, float):
                raise Exception('solver option {} must be of type float. You have {}.'.format(field_, type(value_)))

            double_value = value_
            acados_solver_common.ocp_nlp_solver_opts_set(self.nlp_config, self.nlp_opts, field, <void *> &double_value)

        elif field_ in string_fields:
            if not isinstance(value_, bytes):
                raise Exception('solver option {} must be of type str. You have {}.'.format(field_, type(value_)))

            string_value = value_.encode('utf-8')
            acados_solver_common.ocp_nlp_solver_opts_set(self.nlp_config, self.nlp_opts, field, <void *> &string_value[0])

        else:
            raise Exception('AcadosOcpSolver.options_set() does not support field {}.'\
                '\n Possible values are {}.'.format(field_, ', '.join(int_fields + double_fields + string_fields)))


    def __del__(self):
        if self.solver_created:
            acados_solver.acados_free(self.capsule)
            acados_solver.acados_free_capsule(self.capsule)
