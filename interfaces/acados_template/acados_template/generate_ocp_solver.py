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
from .generate_c_code_constraint import *
from .generate_c_code_nls_cost import *
from .acados_ocp_nlp import *
from .acados_ocp_solver import acados_ocp_solver
from ctypes import *
from copy import deepcopy
from .utils import ACADOS_PATH, is_column, render_template, dict2json, np_array_to_list


def make_ocp_dims_consistent(acados_ocp):
    dims = acados_ocp.dims
    cost = acados_ocp.cost
    constraints = acados_ocp.constraints
    model = acados_ocp.model

    # nbx_0
    if (constraints.lbx_0 == [] and constraints.ubx_0 == []):
        dims.nbx_0 = 0
    elif not (constraints.lbx_0 == [] and constraints.ubx_0 == []):
        this_shape = constraints.lbx_0.shape
        other_shape = constraints.ubx_0.shape
        if not this_shape == other_shape:
            raise Exception("lbx_0, ubx_0 have different shapes!")
        if not is_column(constraints.lbx_0):
            raise Exception("lbx_0, ubx_0 must be column vectors!")

        dims.nbx_0 = constraints.lbx_0.size
    else:
        raise Exception("lbx_0, ubx_0 have different shapes!")

    # nx
    if is_column(model.x):
        dims.nx = model.x.shape[0]
    else:
        raise Exception("model.x should be column vector!")

    # # nu
    # if not model.u == None:
    #     if is_column(model.u):
    #         dims.nu = model.u.shape[0]
    #     else:
    #         raise Exception("model.u should be column vector!")
    # else:
    #     dims.nu = 0

    # # nz
    # if not model.z == None:
    #     print(model.z)
    #     if is_column(model.z):
    #         dims.nz = model.z.shape[0]
    #     else:
    #         raise Exception("model.z should be column vector!")
    # else:
    #     dims.nz = 0


def get_ocp_nlp_layout():
    current_module = sys.modules[__name__]
    acados_path = os.path.dirname(current_module.__file__)
    with open(acados_path + '/acados_layout.json', 'r') as f:
        ocp_nlp_layout = json.load(f)
    return ocp_nlp_layout


def ocp_formulation_json_dump(acados_ocp, json_file='acados_ocp_nlp.json'):
    # Load acados_ocp_nlp structure description
    ocp_layout = get_ocp_nlp_layout()

    # Copy input ocp object dictionary
    ocp_nlp_dict = dict(deepcopy(acados_ocp).__dict__)

    for acados_struct, v in ocp_layout.items():
        # skip non dict attributes
        if not isinstance(v, dict): continue
        #  setattr(ocp_nlp, acados_struct, dict(getattr(acados_ocp, acados_struct).__dict__))
        # Copy ocp object attributes dictionaries
        ocp_nlp_dict[acados_struct]=dict(getattr(acados_ocp, acados_struct).__dict__)

    ocp_nlp_dict['model'] = acados_ocp_model_strip_casadi_symbolics(ocp_nlp_dict['model'])

    ocp_nlp_json = dict2json(ocp_nlp_dict)

    with open(json_file, 'w') as f:
        json.dump(ocp_nlp_json, f, default=np_array_to_list, indent=4, sort_keys=True)



def ocp_formulation_json_load(json_file='acados_ocp_nlp.json'):
    # Load acados_ocp_nlp structure description
    ocp_layout = get_ocp_nlp_layout()

    with open(json_file, 'r') as f:
        ocp_nlp_json = json.load(f)

    ocp_nlp_dict = json2dict(ocp_nlp_json, ocp_nlp_json['dims'])

    # Instantiate acados_ocp_nlp object
    acados_ocp = acados_ocp_nlp()

    # load class dict
    acados_ocp.__dict__ = ocp_nlp_dict

    # laod class attributes dict, dims, constraints, etc
    for acados_struct, v  in ocp_layout.items():
        # skip non dict attributes
        if not isinstance(v, dict): continue
        setattr(acados_ocp, acados_struct, ocp_nlp_as_object(ocp_nlp_dict[acados_struct]))

    return acados_ocp



def generate_ocp_solver(acados_ocp, json_file='acados_ocp_nlp.json',
                        con_h=None, con_hN=None, con_p=None, con_pN=None):

    model = acados_ocp.model
    name = model.name

    # make dims consistent
    make_ocp_dims_consistent(acados_ocp)

    # set integrator time automatically
    acados_ocp.solver_options.Tsim = acados_ocp.solver_options.tf / acados_ocp.dims.N

    # generate external functions
    if acados_ocp.solver_options.integrator_type == 'ERK':
        # explicit model -- generate C code
        generate_c_code_explicit_ode(model)
    else:
        # implicit model -- generate C code
        opts = dict(generate_hess=1)
        generate_c_code_implicit_ode(model, opts)

    if acados_ocp.dims.nphi > 0 or acados_ocp.dims.nh > 0:
        generate_c_code_constraint(acados_ocp.model, name, False)

    if acados_ocp.dims.nphi_e > 0 or acados_ocp.dims.nh_e > 0:
        generate_c_code_constraint(acados_ocp.model, name, True)

    if acados_ocp.cost.cost_type == 'NONLINEAR_LS':
        acados_ocp.cost.Vx = np.zeros((acados_ocp.dims.ny, acados_ocp.dims.nx))
        acados_ocp.cost.Vu = np.zeros((acados_ocp.dims.ny, acados_ocp.dims.nu))
        generate_c_code_nls_cost(acados_ocp.model, name, False)

    if acados_ocp.cost.cost_type_e == 'NONLINEAR_LS':
        acados_ocp.cost.Vx_e = np.zeros((acados_ocp.dims.ny_e, acados_ocp.dims.nx))
        generate_c_code_nls_cost(acados_ocp.model, name, True)

    # dump to json
    ocp_formulation_json_dump(acados_ocp, json_file)

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
    in_file = 'main.in.c'
    out_file = 'main_{}.c'.format(model.name)
    render_template(in_file, out_file, template_dir, json_path)

    in_file = 'acados_solver.in.c'
    out_file = 'acados_solver_{}.c'.format(model.name)
    render_template(in_file, out_file, template_dir, json_path)

    in_file = 'acados_solver.in.h'
    out_file = 'acados_solver_{}.h'.format(model.name)
    render_template(in_file, out_file, template_dir, json_path)

    in_file = 'Makefile.in'
    out_file = 'Makefile'
    render_template(in_file, out_file, template_dir, json_path)

    in_file = 'acados_solver_sfun.in.c'
    out_file = 'acados_solver_sfunction_{}.c'.format(model.name)
    render_template(in_file, out_file, template_dir, json_path)

    in_file = 'make_sfun.in.m'
    out_file = 'make_sfun.m'
    render_template(in_file, out_file, template_dir, json_path)

    in_file = 'acados_sim_solver.in.c'
    out_file = 'acados_sim_solver_{}.c'.format(model.name)
    render_template(in_file, out_file, template_dir, json_path)

    in_file = 'acados_sim_solver.in.h'
    out_file = 'acados_sim_solver_{}.h'.format(model.name)
    render_template(in_file, out_file, template_dir, json_path)

    ## folder model
    template_dir = 'c_generated_code/{}_model/'.format(model.name)

    in_file = 'model.in.h'
    out_file = '{}_model.h'.format(model.name)
    render_template(in_file, out_file, template_dir, json_path)

    # constraints on convex over nonlinear fuction
    if acados_ocp.constraints.constr_type == 'BGP' and acados_ocp.dims.nphi > 0:
        # constraints on outer fuction
        template_dir = 'c_generated_code/{}_constraints/'.format(name)
        in_file = 'phi_constraint.in.h'
        out_file =  '{}_phi_constraint.h'.format(name)
        render_template(in_file, out_file, template_dir, json_path)

    # terminal constraints on convex over nonlinear fuction
    if acados_ocp.constraints.constr_type_e == 'BGP' and acados_ocp.dims.nphi_e > 0:
        # terminal constraints on outer fuction
        template_dir = 'c_generated_code/{}_constraints/'.format(name)
        in_file = 'phi_e_constraint.in.h'
        out_file =  '{}_phi_e_constraint.h'.format(name)
        render_template(in_file, out_file, template_dir, json_path)

    # nonlinear constraints
    if acados_ocp.constraints.constr_type == 'BGH' and acados_ocp.dims.nh > 0:
        template_dir = 'c_generated_code/{}_constraints/'.format(acados_ocp.model.name)
        in_file = 'h_constraint.in.h'
        out_file = '{}_h_constraint.h'.format(acados_ocp.model.name)
        render_template(in_file, out_file, template_dir, json_path)

    # terminal nonlinear constraints
    if acados_ocp.constraints.constr_type_e == 'BGH' and acados_ocp.dims.nh_e > 0:
        template_dir = 'c_generated_code/{}_constraints/'.format(acados_ocp.model.name)
        in_file = 'h_e_constraint.in.h'
        out_file = '{}_h_e_constraint.h'.format(acados_ocp.model.name)
        render_template(in_file, out_file, template_dir, json_path)

    # nonlinear cost function
    if acados_ocp.cost.cost_type == 'NONLINEAR_LS':
        template_dir = 'c_generated_code/{}_r_cost/'.format(name)
        in_file = 'r_cost.in.h'
        out_file = '{}_r_cost.h'.format(name)
        render_template(in_file, out_file, template_dir, json_path)

    # terminal nonlinear cost function
    if acados_ocp.cost.cost_type_e == 'NONLINEAR_LS':
        template_dir = 'c_generated_code/{}_r_e_cost/'.format(name)
        in_file = 'r_e_cost.in.h'
        out_file = '{}_r_e_cost.h'.format(name)
        render_template(in_file, out_file, template_dir, json_path)

    ## Compile solver
    os.chdir('c_generated_code')
    os.system('make clean')
    os.system('make ocp_shared_lib')
    os.chdir('..')

    # get
    ocp_solver = acados_ocp_solver(acados_ocp, 'c_generated_code/libacados_ocp_solver_' + model.name + '.so')
    return ocp_solver
