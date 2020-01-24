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

class AcadosSimSolver:
    def __init__(self, acados_sim, shared_lib):

        self.sim_struct = acados_sim

        model_name = self.sim_struct.model.name

        self.shared_lib = CDLL(shared_lib)
        getattr(self.shared_lib, f"{model_name}_acados_sim_create")()

        getattr(self.shared_lib, f"{model_name}_acados_get_sim_opts").restype = c_void_p
        self.sim_opts = getattr(self.shared_lib, f"{model_name}_acados_get_sim_opts")()

        getattr(self.shared_lib, f"{model_name}_acados_get_sim_dims").restype = c_void_p
        self.sim_dims = getattr(self.shared_lib, f"{model_name}_acados_get_sim_dims")()

        getattr(self.shared_lib, f"{model_name}_acados_get_sim_config").restype = c_void_p
        self.sim_config = getattr(self.shared_lib, f"{model_name}_acados_get_sim_config")()

        getattr(self.shared_lib, f"{model_name}_acados_get_sim_out").restype = c_void_p
        self.sim_out = getattr(self.shared_lib, f"{model_name}_acados_get_sim_out")()

        getattr(self.shared_lib, f"{model_name}_acados_get_sim_in").restype = c_void_p
        self.sim_in = getattr(self.shared_lib, f"{model_name}_acados_get_sim_in")()

        nu = self.sim_struct.dims.nu
        nx = self.sim_struct.dims.nx
        self.gettable = {
            'x': nx,
            'xn': nx,
            'u': nu,
            'S_forw': nx*(nx+nu),
            'Sx': nx*nx
        }

        self.settable = ['S_adj', 'S_forw', 'T', 'x', 'u', 'xdot', 'z', 'Su', 'Sx', 'p']
        self.model_name = model_name

    def solve(self):
        status = getattr(self.shared_lib, f"{self.model_name}_acados_sim_solve")()
        return status

    def get(self, field_):

        field = field_
        field = field.encode('utf-8')

        if field_ in self.gettable.keys():

            # allocate array
            dims = self.gettable[field_]
            out = np.ascontiguousarray(np.zeros((dims,)), dtype=np.float64)
            out_data = cast(out.ctypes.data, POINTER(c_double))

            self.shared_lib.sim_out_get.argtypes = [c_void_p, c_void_p, c_void_p, c_char_p, c_void_p]
            self.shared_lib.sim_out_get(self.sim_config, self.sim_dims, self.sim_out, field, out_data)

        else:
            raise Exception(f'acados_solver.set(): Unknown field {field}, available fiels are {",".join(self.gettable.keys())}')

        # out = cast((out), POINTER(c_double))

        return out

    def set(self, field_, value_):

        # cast value_ to avoid conversion issues
        if type(value_) == float:
            value_ = np.array([value_])

        value_ = value_.astype(float)
        value_data = cast(value_.ctypes.data, POINTER(c_double))
        value_data_p = cast((value_data), c_void_p)

        field = field_
        field = field.encode('utf-8')

        # treat parameters separately
        if field_ is 'p':
            model_name = self.sim_struct.model.name
            getattr(self.shared_lib, f"{model_name}_acados_sim_update_params").argtypes = [POINTER(c_double)]
            value_data = cast(value_.ctypes.data, POINTER(c_double))
            getattr(self.shared_lib, f"{model_name}_acados_sim_update_params")(value_data, value_.shape[0])

        elif field_ in self.settable:
            self.shared_lib.sim_in_set.argtypes = [c_void_p, c_void_p, c_void_p, c_char_p, c_void_p]
            self.shared_lib.sim_in_set(self.sim_config, self.sim_dims, self.sim_in, field, value_data_p)
        else:
            raise Exception(f'acados_solver.set(): Unknown field {field}, available fiels are {",".join(self.settable)}')

        return

    def __del__(self):
        getattr(self.shared_lib, f"{self.model_name}_acados_sim_free")()
