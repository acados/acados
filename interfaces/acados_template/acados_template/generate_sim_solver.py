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

from .generate_c_code_explicit_ode import *
from .generate_c_code_implicit_ode import *
from .acados_sim import *
from .acados_sim_solver import acados_sim_solver
from ctypes import *
from copy import deepcopy
from .utils import ACADOS_PATH, is_column, render_template, dict2json, np_array_to_list


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

    sim_dict['model'] = acados_dae_strip_non_num(sim_dict['model'])
    sim_json = dict2json(sim_dict)

    with open(json_file, 'w') as f:
        json.dump(sim_json, f, default=np_array_to_list, indent=4, sort_keys=True)


def sim_render_templates(json_file, model_name):
    # setting up loader and environment
    template_glob = os.path.join(
        ACADOS_PATH,
        'interfaces/acados_template/acados_template/c_templates_tera/*')
        # TODO(andrea): this breaks when running main_test.py...
        # 'interfaces/acados_template/acados_template/c_templates_tera/*!(swp.*)')
    acados_template_path = os.path.join(
        ACADOS_PATH,
        'interfaces/acados_template/acados_template/c_templates_tera')

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

    ## folder model
    template_dir = 'c_generated_code/{}_model/'.format(model_name)

    in_file = 'model.in.h'
    out_file = '{}_model.h'.format(model_name)
    render_template(in_file, out_file, template_dir, json_path)


def sim_generate_casadi_functions(acados_sim):
    model = acados_sim.model
    integrator_type = acados_sim.solver_options.integrator_type
    # generate external functions
    if integrator_type == 'ERK':
        # explicit model -- generate C code
        generate_c_code_explicit_ode(model)
    elif integrator_type == 'IRK':
        # implicit model -- generate C code
        opts = dict(generate_hess=1)
        generate_c_code_implicit_ode(model, opts)


def generate_sim_solver(acados_sim, json_file='acados_sim.json'):

    model_name = acados_sim.model.name

    make_sim_dims_consistent(acados_sim)

    sim_formulation_json_dump(acados_sim, json_file)

    # render templates
    sim_render_templates(json_file, model_name)
    # generate casadi functions
    sim_generate_casadi_functions(acados_sim)

    ## Compile solver
    os.chdir('c_generated_code')
    os.system('make sim_shared_lib')
    os.chdir('..')

    # get
    sim_solver = acados_sim_solver(acados_sim, 'c_generated_code/libacados_sim_solver_' + model_name + '.so')
    return sim_solver


def generate_sim_solver_from_ocp(acados_ocp, json_file='acados_ocp_nlp.json'):

    model_name = acados_ocp.model.name

    # set up acados_sim_
    acados_sim_ = acados_sim()
    acados_sim_.model = acados_ocp.model
    acados_sim_.dims.nx = acados_ocp.dims.nx
    acados_sim_.dims.nu = acados_ocp.dims.nu
    acados_sim_.dims.nz = acados_ocp.dims.nz
    acados_sim_.dims.np = acados_ocp.dims.np

    acados_sim_.solver_options.integrator_type = acados_ocp.solver_options.integrator_type

    # render templates
    sim_render_templates(json_file, model_name)
    # generate casadi functions
    sim_generate_casadi_functions(acados_sim_)

    ## Compile solver
    os.chdir('c_generated_code')
    os.system('make sim_shared_lib')
    os.chdir('..')

    # get
    sim_solver = acados_sim_solver(acados_sim_, 'c_generated_code/libacados_sim_solver_' + model_name + '.so')
    return sim_solver
