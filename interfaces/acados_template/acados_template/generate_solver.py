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
from .generate_c_code_constraint_e import *
from .generate_c_code_nls_cost import *
from .generate_c_code_nls_cost_e import *
from .acados_ocp_nlp import *
from .acados_ocp_solver import acados_ocp_solver
from .acados_sim_solver import acados_sim_solver
from ctypes import *
from copy import deepcopy
from .utils import ACADOS_PATH, get_tera, is_column

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

def store_ocp_solver(acados_ocp, json_file='acados_ocp_nlp.json'):
    # Load acados_ocp_nlp structure description
    ocp_layout = get_ocp_nlp_layout()

    # Instatiate acados_ocp_nlp object
    ocp_nlp = acados_ocp_nlp()

    # Copy input ocp object dictionary
    ocp_nlp_dict = dict(acados_ocp.__dict__)

    for acados_struct, v  in ocp_layout.items():
        # skip non dict attributes
        if not isinstance(v, dict): continue
        #  setattr(ocp_nlp, acados_struct, dict(getattr(acados_ocp, acados_struct).__dict__))
        # Copy ocp object attributes dictionaries
        ocp_nlp_dict[acados_struct]=dict(getattr(acados_ocp, acados_struct).__dict__)

    #  ocp_nlp = acados_ocp
    #  ocp_nlp.cost = acados_ocp.cost.__dict__
    #  ocp_nlp.constraints = acados_ocp.constraints.__dict__
    #  ocp_nlp.solver_config = acados_ocp.solver_config.__dict__
    #  ocp_nlp.dims = acados_ocp.dims.__dict__
    #  ocp_nlp = ocp_nlp.__dict__

    ocp_nlp_json = dict2json(ocp_nlp_dict)

    with open(json_file, 'w') as f:
        json.dump(ocp_nlp_json, f, default=np_array_to_list)

def load_ocp_solver(json_file='acados_ocp_nlp.json'):
    # Load acados_ocp_nlp structure description
    ocp_layout = get_ocp_nlp_layout()

    with open(json_file, 'r') as f:
        ocp_nlp_json = json.load(f)

    ocp_nlp_dict = json2dict(ocp_nlp_json, ocp_nlp_json['dims'])

    # Instantiate acados_ocp_nlp object
    acados_ocp = acados_ocp_nlp()

    # load class dict
    acados_ocp.__dict__ = ocp_nlp_dict

    # laod class attributes dict, dims, constraints, ecc
    for acados_struct, v  in ocp_layout.items():
        # skip non dict attributes
        if not isinstance(v, dict): continue
        setattr(acados_ocp, acados_struct, ocp_nlp_as_object(ocp_nlp_dict[acados_struct]))

    #  acados_ocp = ocp_nlp_as_object(ocp_nlp_dict)
    #  acados_ocp.cost = ocp_nlp_as_object(acados_ocp.cost)
    #  acados_ocp.constraints = ocp_nlp_as_object(acados_ocp.constraints)
    #  acados_ocp.solver_config = ocp_nlp_as_object(acados_ocp.solver_config)
    #  acados_ocp.dims = ocp_nlp_as_object(acados_ocp.dims)

    return acados_ocp

def generate_solvers(acados_ocp, model, con_h=None, con_hN=None, con_p=None, con_pN=None):

    model = acados_ocp.model
    name = model.name

    # make dims consistent
    make_ocp_dims_consistent(acados_ocp)

    acados_ocp = ocp_nlp_as_object(ocp_nlp_dict)
    acados_ocp.cost = ocp_nlp_as_object(acados_ocp.cost)
    acados_ocp.constraints = ocp_nlp_as_object(acados_ocp.constraints)
    acados_ocp.solver_config = ocp_nlp_as_object(acados_ocp.solver_config)
    acados_ocp.dims = ocp_nlp_as_object(acados_ocp.dims)

    # setting up loader and environment
    acados_path = os.path.dirname(os.path.abspath(__file__))

    template_glob = acados_path + '/c_templates_tera/*'
    acados_template_path = acados_path + '/c_templates_tera'

    if acados_ocp.solver_config.integrator_type == 'ERK':
        # explicit model -- generate C code
        generate_c_code_explicit_ode(model)
    else:
        # implicit model -- generate C code
        opts = dict(generate_hess=1)
        generate_c_code_implicit_ode(model, opts)

    if con_p is not None and con_h is None:
        raise Exception('h constraint is missing!')

    if con_h is not None:
        # nonlinear part of nonlinear constraints
        generate_c_code_constraint_e(acados_ocp.con_phi_e, name)
    elif acados_ocp.constraints.constr_type_e  == 'BGH' and acados_ocp.dims.nh_e > 0:
        generate_c_code_constraint_e(acados_ocp.con_h_e, name)

    if acados_ocp.cost.cost_type == 'NONLINEAR_LS':
        acados_ocp.cost.Vx = np.zeros((acados_ocp.dims.ny, acados_ocp.dims.nx))
        acados_ocp.cost.Vu = np.zeros((acados_ocp.dims.ny, acados_ocp.dims.nu))
        generate_c_code_nls_cost(acados_ocp.cost_r, name)

    if acados_ocp.cost.cost_type_e == 'NONLINEAR_LS':
        acados_ocp.cost.Vx_e = np.zeros((acados_ocp.dims.ny_e, acados_ocp.dims.nx))
        generate_c_code_nls_cost_e(acados_ocp.cost_r_e, name)

    ocp_nlp = deepcopy(acados_ocp)
    ocp_nlp.cost = acados_ocp.cost.__dict__
    ocp_nlp.constraints = acados_ocp.constraints.__dict__
    ocp_nlp.solver_options = acados_ocp.solver_options.__dict__
    ocp_nlp.dims = acados_ocp.dims.__dict__
    ocp_nlp.con_h = acados_ocp.con_h.__dict__
    ocp_nlp.con_h_e = acados_ocp.con_h_e.__dict__
    ocp_nlp.con_phi = acados_ocp.con_phi.__dict__
    ocp_nlp.con_phi_e = acados_ocp.con_phi_e.__dict__
    ocp_nlp.cost_r = acados_ocp.cost_r.__dict__
    ocp_nlp.cost_r_e = acados_ocp.cost_r_e.__dict__
    ocp_nlp.model = acados_ocp.model.__dict__
    ocp_nlp = ocp_nlp.__dict__

    # need to strip non-numerical stuff from expressions for now
    ocp_nlp['con_h'] = acados_constraint_strip_non_num(ocp_nlp['con_h'])
    ocp_nlp['con_h_e'] = acados_constraint_strip_non_num(ocp_nlp['con_h_e'])
    ocp_nlp['con_phi'] = acados_constraint_strip_non_num(ocp_nlp['con_phi'])
    ocp_nlp['con_phi_e'] = acados_constraint_strip_non_num(ocp_nlp['con_phi_e'])

    ocp_nlp['model'] = acados_dae_strip_non_num(ocp_nlp['model'])

    ocp_nlp['cost_r'] = acados_cost_strip_non_num(ocp_nlp['cost_r'])
    ocp_nlp['cost_r_e'] = acados_cost_strip_non_num(ocp_nlp['cost_r_e'])

    ocp_nlp = dict2json(ocp_nlp)

    with open(json_file, 'w') as f:
        json.dump(ocp_nlp, f, default=np_array_to_list, indent=4, sort_keys=True)

    with open(json_file, 'r') as f:
        ocp_nlp_json = json.load(f)

    ocp_nlp_dict = json2dict(ocp_nlp_json, ocp_nlp_json['dims'])

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

    def render_template(in_file, out_file, template_dir):
        cwd = os.getcwd()
        if not os.path.exists(template_dir):
            os.mkdir(template_dir)
        os.chdir(template_dir)

        # call tera as system cmd
        os_cmd = "{tera_path} '{template_glob}' '{in_file}' '{json_path}' '{out_file}'".format(
            tera_path=tera_path,
            template_glob=template_glob,
            json_path=json_path,
            in_file=in_file,
            out_file=out_file
        )
        status = os.system(os_cmd)
        if (status != 0):
            raise Exception('Rendering of {} failed! Exiting.\n'.format(in_file))

        os.chdir(cwd)

    ## folder: c_generated_code
    template_dir = 'c_generated_code/'

    in_file = 'main.in.c'
    out_file = 'main_{}.c'.format(model.name)
    render_template(in_file, out_file, template_dir)

    in_file = 'acados_solver.in.c'
    out_file = 'acados_solver_{}.c'.format(model.name)
    render_template(in_file, out_file, template_dir)

    in_file = 'acados_solver.in.h'
    out_file = 'acados_solver_{}.h'.format(model.name)
    render_template(in_file, out_file, template_dir)

    in_file = 'acados_sim_solver.in.c'
    out_file = 'acados_sim_solver_{}.c'.format(model.name)
    render_template(in_file, out_file, template_dir)

    in_file = 'acados_sim_solver.in.h'
    out_file = 'acados_sim_solver_{}.h'.format(model.name)
    render_template(in_file, out_file, template_dir)

    in_file = 'Makefile.in'
    out_file = 'Makefile'
    render_template(in_file, out_file, template_dir)

    in_file = 'acados_solver_sfun.in.c'
    out_file = 'acados_solver_sfunction_{}.c'.format(model.name)
    render_template(in_file, out_file, template_dir)

    in_file = 'make_sfun.in.m'
    out_file = 'make_sfun.m'
    render_template(in_file, out_file, template_dir)

    ## folder model
    template_dir = 'c_generated_code/{}_model/'.format(model.name)

    in_file = 'model.in.h'
    out_file = '{}_model.h'.format(model.name)
    render_template(in_file, out_file, template_dir)

    # constraints on convex over nonlinear fuction
    if acados_ocp.constraints.constr_type == 'BGP' and acados_ocp.dims.nphi > 0:
        # constraints on outer fuction
        template_dir = 'c_generated_code/{}_constraints/'.format(name)
        in_file = 'phi_constraint.in.h'
        out_file =  '{}_phi_constraint.h'.format(name)
        render_template(in_file, out_file, template_dir)

        # constraints on inner fuction
        template_dir = 'c_generated_code/{}_constraints/'.format(name)
        in_file = 'r_constraint.in.h'
        out_file = '{}_r_constraint.h'.format(name)
        render_template(in_file, out_file, template_dir)

    # terminal constraints on convex over nonlinear fuction
    if acados_ocp.constraints.constr_type_e == 'BGP' and acados_ocp.dims.nphi_e > 0:
        # terminal constraints on outer fuction
        template_dir = 'c_generated_code/{}_constraints/'.format(name)
        in_file = 'phi_e_constraint.in.h'
        out_file =  '{}_phi_e_constraint.h'.format(name)
        render_template(in_file, out_file, template_dir)

        # terminal constraints on inner function
        template_dir = 'c_generated_code/{}_constraints/'.format(name)
        in_file = 'r_e_constraint.in.h'
        out_file = '{}_r_e_constraint.h'.format(name)
        render_template(in_file, out_file, template_dir)

    # nonlinear constraints
    if acados_ocp.constraints.constr_type == 'BGH' and acados_ocp.dims.nh > 0:
        template_dir = 'c_generated_code/{}_constraints/'.format(acados_ocp.model.name)
        in_file = 'h_constraint.in.h'
        out_file = '{}_h_constraint.h'.format(acados_ocp.model.name)
        render_template(in_file, out_file, template_dir)

    # terminal nonlinear constraints
    if acados_ocp.constraints.constr_type_e == 'BGH' and acados_ocp.dims.nh_e > 0:
        template_dir = 'c_generated_code/{}_constraints/'.format(acados_ocp.model.name)
        in_file = 'h_e_constraint.in.h'
        out_file = '{}_h_e_constraint.h'.format(acados_ocp.model.name)
        render_template(in_file, out_file, template_dir)

    # nonlinear cost function
    if acados_ocp.cost.cost_type == 'NONLINEAR_LS':
        template_dir = 'c_generated_code/{}_r_cost/'.format(name)
        in_file = 'r_cost.in.h'
        out_file = '{}_r_cost.h'.format(name)
        render_template(in_file, out_file, template_dir)

    # terminal nonlinear cost function
    if acados_ocp.cost.cost_type_e == 'NONLINEAR_LS':
        template_dir = 'c_generated_code/{}_r_e_cost/'.format(name)
        in_file = 'r_e_cost.in.h'
        out_file = '{}_r_e_cost.h'.format(name)
        render_template(in_file, out_file, template_dir)

    ## Compile solver
    os.chdir('c_generated_code')
    os.system('make')
    os.system('make shared_lib')
    os.chdir('..')

    ocp_solver = acados_ocp_solver(acados_ocp, 'c_generated_code/libacados_ocp_solver_' + model.name + '.so')
    sim_solver = acados_sim_solver(acados_ocp, 'c_generated_code/libacados_sim_solver_' + model.name + '.so')
    return ocp_solver, sim_solver
