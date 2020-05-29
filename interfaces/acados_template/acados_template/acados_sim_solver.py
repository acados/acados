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

import sys, os, json

import numpy as np

from ctypes import *
from copy import deepcopy

from .generate_c_code_explicit_ode import generate_c_code_explicit_ode
from .generate_c_code_implicit_ode import generate_c_code_implicit_ode
from .generate_c_code_gnsf import generate_c_code_gnsf
from .acados_sim import AcadosSim
from .acados_ocp import AcadosOcp
from .acados_model import acados_model_strip_casadi_symbolics
from .utils import is_column, render_template, format_class_dict, np_array_to_list,\
     make_model_consistent, set_up_imported_gnsf_model


def make_sim_dims_consistent(acados_sim):
    dims = acados_sim.dims
    model = acados_sim.model
    # nx
    if is_column(model.x):
        dims.nx = model.x.shape[0]
    else:
        raise Exception("model.x should be column vector!")

    # nu
    if is_column(model.u):
        dims.nu = model.u.shape[0]
    elif model.u == None or model.u == []:
        dims.nu = 0
    else:
        raise Exception("model.u should be column vector or None!")

    # nz
    if is_column(model.z):
        dims.nz = model.z.shape[0]
    elif model.z == None or model.z == []:
        dims.nz = 0
    else:
        raise Exception("model.z should be column vector or None!")

    # np
    if is_column(model.p):
        dims.np = model.p.shape[0]
    elif model.p == None or model.p == []:
        dims.np = 0
    else:
        raise Exception("model.p should be column vector or None!")


def get_sim_layout():
    current_module = sys.modules[__name__]
    acados_path = os.path.dirname(current_module.__file__)
    with open(acados_path + '/acados_sim_layout.json', 'r') as f:
        sim_layout = json.load(f)
    return sim_layout



def sim_formulation_json_dump(acados_sim, json_file='acados_sim.json'):
    # Load acados_sim structure description
    sim_layout = get_sim_layout()

    # Copy input sim object dictionary
    sim_dict = dict(deepcopy(acados_sim).__dict__)

    for key, v in sim_layout.items():
        # skip non dict attributes
        if not isinstance(v, dict): continue
        # Copy sim object attributes dictionaries
        sim_dict[key]=dict(getattr(acados_sim, key).__dict__)

    sim_dict['model'] = acados_model_strip_casadi_symbolics(sim_dict['model'])
    sim_json = format_class_dict(sim_dict)

    with open(json_file, 'w') as f:
        json.dump(sim_json, f, default=np_array_to_list, indent=4, sort_keys=True)


def sim_render_templates(json_file, model_name):
    # setting up loader and environment
    json_path = '{cwd}/{json_file}'.format(
        cwd=os.getcwd(),
        json_file=json_file)

    if not os.path.exists(json_path):
        raise Exception("{} not found!".format(json_path))

    template_dir = 'c_generated_code/'

    ## Render templates
    in_file = 'acados_sim_solver.in.c'
    out_file = 'acados_sim_solver_{}.c'.format(model_name)
    render_template(in_file, out_file, template_dir, json_path)

    in_file = 'acados_sim_solver.in.h'
    out_file = 'acados_sim_solver_{}.h'.format(model_name)
    render_template(in_file, out_file, template_dir, json_path)

    in_file = 'Makefile.in'
    out_file = 'Makefile'
    render_template(in_file, out_file, template_dir, json_path)

    in_file = 'main_sim.in.c'
    out_file = 'main_sim_{}.c'.format(model_name)
    render_template(in_file, out_file, template_dir, json_path)

    ## folder model
    template_dir = 'c_generated_code/{}_model/'.format(model_name)

    in_file = 'model.in.h'
    out_file = '{}_model.h'.format(model_name)
    render_template(in_file, out_file, template_dir, json_path)


def sim_generate_casadi_functions(acados_sim):
    model = acados_sim.model
    model = make_model_consistent(model)

    integrator_type = acados_sim.solver_options.integrator_type
    # generate external functions
    if integrator_type == 'ERK':
        # explicit model -- generate C code
        generate_c_code_explicit_ode(model)
    elif integrator_type == 'IRK':
        # implicit model -- generate C code
        opts = dict(generate_hess=1)
        generate_c_code_implicit_ode(model, opts)
    elif integrator_type == 'GNSF':
        generate_c_code_gnsf(model)

class AcadosSimSolver:
    def __init__(self, acados_sim_, json_file='acados_sim.json'):

        if isinstance(acados_sim_, AcadosOcp):
            # set up acados_sim_
            acados_sim = AcadosSim()
            acados_sim.model = acados_sim_.model
            acados_sim.dims.nx = acados_sim_.dims.nx
            acados_sim.dims.nu = acados_sim_.dims.nu
            acados_sim.dims.nz = acados_sim_.dims.nz
            acados_sim.dims.np = acados_sim_.dims.np
            acados_sim.solver_options.integrator_type = acados_sim_.solver_options.integrator_type

        elif isinstance(acados_sim_, AcadosSim):
            acados_sim = acados_sim_

        acados_sim.__problem_class = 'SIM'

        model_name = acados_sim.model.name
        make_sim_dims_consistent(acados_sim)

        # reuse existing json and casadi functions, when creating integrator from ocp
        if isinstance(acados_sim_, AcadosSim):
            if acados_sim.solver_options.integrator_type == 'GNSF':
                set_up_imported_gnsf_model(acados_sim)

            sim_generate_casadi_functions(acados_sim)
            sim_formulation_json_dump(acados_sim, json_file)

        # render templates
        sim_render_templates(json_file, model_name)

        ## Compile solver
        os.chdir('c_generated_code')
        os.system('make sim_shared_lib')
        os.chdir('..')

        # Ctypes
        shared_lib = 'c_generated_code/libacados_sim_solver_' + model_name + '.so'

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
            'Sx': nx*nx,
            'Su': nx*nu,
            'S_adj': nx+nu,
            'S_hess': (nx+nu)*(nx+nu),
        }

        self.settable = ['S_adj', 'T', 'x', 'u', 'xdot', 'z', 'p'] # S_forw
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

            if field_ == 'S_forw':
                nu = self.sim_struct.dims.nu
                nx = self.sim_struct.dims.nx
                out = out.reshape(nx, nx+nu, order='F')
            elif field_ == 'Sx':
                nx = self.sim_struct.dims.nx
                out = out.reshape(nx, nx, order='F')
            elif field_ == 'Sx':
                nx = self.sim_struct.dims.nx
                nu = self.sim_struct.dims.nu
                out = out.reshape(nx, nu, order='F')
            elif field_ == 'S_hess':
                nx = self.sim_struct.dims.nx
                nu = self.sim_struct.dims.nu
                out = out.reshape(nx+nu, nx+nu, order='F')
        else:
            raise Exception(f'AcadosSimSolver.get(): Unknown field {field_},' \
                f' available fields are {", ".join(self.gettable.keys())}')

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
            raise Exception(f'AcadosSimSolver.set(): Unknown field {field_},' \
                f' available fields are {", ".join(self.settable)}')

        return


    def __del__(self):
        getattr(self.shared_lib, f"{self.model_name}_acados_sim_free")()
