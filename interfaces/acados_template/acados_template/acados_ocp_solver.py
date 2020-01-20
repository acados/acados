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

from ctypes import *
import numpy as np

class acados_ocp_solver:
    def __init__(self, acados_ocp, shared_lib):
        self.shared_lib = CDLL(shared_lib)
        self.shared_lib.acados_create()

        self.shared_lib.acados_get_nlp_opts.restype = c_void_p
        self.nlp_opts = self.shared_lib.acados_get_nlp_opts()

        self.shared_lib.acados_get_nlp_dims.restype = c_void_p
        self.nlp_dims = self.shared_lib.acados_get_nlp_dims()

        self.shared_lib.acados_get_nlp_config.restype = c_void_p
        self.nlp_config = self.shared_lib.acados_get_nlp_config()

        self.shared_lib.acados_get_nlp_out.restype = c_void_p
        self.nlp_out = self.shared_lib.acados_get_nlp_out()

        self.shared_lib.acados_get_nlp_in.restype = c_void_p
        self.nlp_in = self.shared_lib.acados_get_nlp_in()

        self.shared_lib.acados_get_nlp_solver.restype = c_void_p
        self.nlp_solver = self.shared_lib.acados_get_nlp_solver()

        self.acados_ocp = acados_ocp

    def solve(self):
        status = self.shared_lib.acados_solve()
        return status


    def get(self, stage_, field_):

        out_fields = ['x', 'u', 'z']
        field = field_
        field = field.encode('utf-8')

        if (field_ not in out_fields):
            raise Exception("acados_solver: {} is not a valid key for method `set(value)`.\
                    \n Possible values are {}. Exiting.".format(out_fields))

        self.shared_lib.ocp_nlp_dims_get_from_attr.argtypes = \
            [c_void_p, c_void_p, c_void_p, c_int, c_char_p]
        self.shared_lib.ocp_nlp_dims_get_from_attr.restype = c_int

        dims = self.shared_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, \
            self.nlp_dims, self.nlp_out, stage_, field)

        out = np.ascontiguousarray(np.zeros((dims,)), dtype=np.float64)
        out_data = cast(out.ctypes.data, POINTER(c_double))

        self.shared_lib.ocp_nlp_out_get.argtypes = \
            [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]
        self.shared_lib.ocp_nlp_out_get(self.nlp_config, \
            self.nlp_dims, self.nlp_out, stage_, field, out_data)

        # out = cast((out), POINTER(c_double))

        return out

    def get_stats(self, field_):
        fields = ['time_tot']
        field = field_
        field = field.encode('utf-8')
        if (field_ not in fields):
            raise Exception("acados_solver: {} is not a valid key for method `set(value)`.\
                    \n Possible values are {}. Exiting.".format(fields, fields))
        out = np.ascontiguousarray(np.zeros((1,)), dtype=np.float64)
        out_data = cast(out.ctypes.data, POINTER(c_double))

        self.shared_lib.ocp_nlp_get.argtypes = [c_void_p, c_void_p, c_char_p, c_void_p]
        self.shared_lib.ocp_nlp_get(self.nlp_config, self.nlp_solver, field, out_data)

        return out

    # Note: this function should not be used anymore, better use cost_set, constraints_set
    def set(self, stage_, field_, value_):

        cost_fields = ['y_ref', 'yref']
        constraints_fields = ['lbx', 'ubx', 'lbu', 'ubu']
        out_fields = ['x', 'u']

        # cast value_ to avoid conversion issues
        value_ = value_.astype(float)

        field = field_
        field = field.encode('utf-8')

        stage = c_int(stage_)

        # treat parameters separately
        if field_ is 'p':
            # not setting parameters
            self.shared_lib.acados_update_params.argtypes = [c_int, POINTER(c_double)]
            self.shared_lib.acados_update_params.restype = c_int
            value_data = cast(value_.ctypes.data, POINTER(c_double))
            self.shared_lib.acados_update_params(stage, value_data, value_.shape[0])
        else:
            if (field_ not in constraints_fields) and \
                    (field_ not in cost_fields) and (field_ not in out_fields):
                raise Exception("acados_solver: {} is not a valid key for method `set(value)`.\
                    \nPossible values are {} and {}. Exiting.".format(field, \
                    cost_fields, constraints_fields, out_fields))

            self.shared_lib.ocp_nlp_dims_get_from_attr.argtypes = \
                [c_void_p, c_void_p, c_void_p, c_int, c_char_p]
            self.shared_lib.ocp_nlp_dims_get_from_attr.restype = c_int

            dims = self.shared_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, \
                self.nlp_dims, self.nlp_out, stage_, field)

            if value_.shape[0] != dims:
                msg = 'acados_solver.set(): mismatching dimension for field "{}"'.format(field_)
                msg += 'with dimension {} (you have {})'.format(dims, value_.shape[0])
                raise Exception(msg)

            value_data = cast(value_.ctypes.data, POINTER(c_double))
            value_data_p = cast((value_data), c_void_p)

            if field_ in constraints_fields:
                self.shared_lib.ocp_nlp_constraints_model_set.argtypes = \
                    [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]
                self.shared_lib.ocp_nlp_constraints_model_set(self.nlp_config, \
                    self.nlp_dims, self.nlp_in, stage, field, value_data_p)
            elif field_ in cost_fields:
                self.shared_lib.ocp_nlp_cost_model_set.argtypes = \
                    [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]
                self.shared_lib.ocp_nlp_cost_model_set(self.nlp_config, \
                    self.nlp_dims, self.nlp_in, stage, field, value_data_p)
            elif field_ in out_fields:
                self.shared_lib.ocp_nlp_out_set.argtypes = \
                    [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]
                self.shared_lib.ocp_nlp_out_set(self.nlp_config, \
                    self.nlp_dims, self.nlp_out, stage, field, value_data_p)

        return


    def cost_set(self, stage_, field_, value_):
        # cast value_ to avoid conversion issues
        value_ = value_.astype(float)

        field = field_
        field = field.encode('utf-8')

        stage = c_int(stage_)
        self.shared_lib.ocp_nlp_cost_dims_get_from_attr.argtypes = \
            [c_void_p, c_void_p, c_void_p, c_int, c_char_p, POINTER(c_int)]
        self.shared_lib.ocp_nlp_cost_dims_get_from_attr.restype = c_int

        dims = np.ascontiguousarray(np.zeros((2,)), dtype=np.intc)
        dims_data = cast(dims.ctypes.data, POINTER(c_int))

        self.shared_lib.ocp_nlp_cost_dims_get_from_attr(self.nlp_config, \
            self.nlp_dims, self.nlp_out, stage_, field, dims_data)

        value_shape = value_.shape
        if len(value_shape) == 1:
            value_shape = (value_shape[0], 0)

        if value_shape != tuple(dims):
            raise Exception('acados_solver.set(): mismatching dimension', \
                ' for field "{}" with dimension {} (you have {})'.format( \
                field_, tuple(dims), value_shape))

        value_data = cast(value_.ctypes.data, POINTER(c_double))
        value_data_p = cast((value_data), c_void_p)

        self.shared_lib.ocp_nlp_cost_model_set.argtypes = \
            [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]
        self.shared_lib.ocp_nlp_cost_model_set(self.nlp_config, \
            self.nlp_dims, self.nlp_in, stage, field, value_data_p)

        return


    def constraints_set(self, stage_, field_, value_):
        # cast value_ to avoid conversion issues
        value_ = value_.astype(float)

        field = field_
        field = field.encode('utf-8')

        stage = c_int(stage_)
        self.shared_lib.ocp_nlp_constraint_dims_get_from_attr.argtypes = \
            [c_void_p, c_void_p, c_void_p, c_int, c_char_p, POINTER(c_int)]
        self.shared_lib.ocp_nlp_constraint_dims_get_from_attr.restype = c_int

        dims = np.ascontiguousarray(np.zeros((2,)), dtype=np.intc)
        dims_data = cast(dims.ctypes.data, POINTER(c_int))

        self.shared_lib.ocp_nlp_constraint_dims_get_from_attr(self.nlp_config, \
            self.nlp_dims, self.nlp_out, stage_, field, dims_data)

        value_shape = value_.shape
        if len(value_shape) == 1:
            value_shape = (value_shape[0], 0)

        if value_shape != tuple(dims):
            raise Exception('acados_solver.set(): mismatching dimension' \
                ' for field "{}" with dimension {} (you have {})'.format(field_, tuple(dims), value_shape))

        value_data = cast(value_.ctypes.data, POINTER(c_double))
        value_data_p = cast((value_data), c_void_p)

        self.shared_lib.ocp_nlp_constraints_model_set.argtypes = \
            [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]
        self.shared_lib.ocp_nlp_constraints_model_set(self.nlp_config, \
            self.nlp_dims, self.nlp_in, stage, field, value_data_p)

        return


    def __del__(self):
        self.shared_lib.acados_free()
