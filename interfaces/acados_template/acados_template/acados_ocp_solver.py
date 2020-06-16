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
from casadi import CasadiMeta, Function, SX

from copy import deepcopy

from .generate_c_code_explicit_ode import generate_c_code_explicit_ode
from .generate_c_code_implicit_ode import generate_c_code_implicit_ode
from .generate_c_code_gnsf import generate_c_code_gnsf
from .generate_c_code_constraint import generate_c_code_constraint
from .generate_c_code_nls_cost import generate_c_code_nls_cost
from .generate_c_code_external_cost import generate_c_code_external_cost
from .acados_ocp import AcadosOcp
from .acados_model import acados_model_strip_casadi_symbolics
from .utils import is_column, is_empty, casadi_length, render_template, acados_class2dict,\
     format_class_dict, ocp_check_against_layout, np_array_to_list, make_model_consistent,\
     set_up_imported_gnsf_model


def make_ocp_dims_consistent(acados_ocp):
    dims = acados_ocp.dims
    cost = acados_ocp.cost
    constraints = acados_ocp.constraints
    model = acados_ocp.model
    opts = acados_ocp.solver_options

    # nx
    if is_column(model.x):
        dims.nx = casadi_length(model.x)
    else:
        raise Exception('model.x should be column vector!')

    # nu
    if is_empty(model.u):
        dims.nu = 0
    else:
        dims.nu = casadi_length(model.u)

    # nz
    if is_empty(model.z):
        dims.nz = 0
    else:
        dims.nz = casadi_length(model.z)

    # np
    if is_empty(model.p):
        dims.np = 0
    else:
        dims.np = casadi_length(model.p)
    if acados_ocp.parameter_values.shape[0] != dims.np:
        raise Exception('inconsistent dimension np, regarding model.p and parameter_values.')

    ## cost
    # path
    if cost.cost_type == 'LINEAR_LS':
        ny = cost.W.shape[0]
        if cost.Vx.shape[0] != ny or cost.Vu.shape[0] != ny:
            raise Exception('inconsistent dimension ny, regarding W, Vx, Vu.' + \
                            f'\nGot W[{cost.W.shape}], Vx[{cost.Vx.shape}], Vu[{cost.Vu.shape}]\n')
        if dims.nz != 0 and cost.Vz.shape[0] != ny:
            raise Exception('inconsistent dimension ny, regarding W, Vx, Vu, Vz.' + \
                            f'\nGot W[{cost.W.shape}], Vx[{cost.Vx.shape}], Vu[{cost.Vu.shape}], Vz[{cost.Vz.shape}]\n')
        if cost.Vx.shape[1] != dims.nx and ny != 0:
            raise Exception('inconsistent dimension: Vx should have nx columns.')
        if cost.Vu.shape[1] != dims.nu and ny != 0:
            raise Exception('inconsistent dimension: Vu should have nu columns.')
        if cost.yref.shape[0] != ny:
            raise Exception('inconsistent dimension: regarding W, yref.' + \
                            f'\nGot W[{cost.W.shape}], yref[{cost.yref.shape}]\n')
        dims.ny = ny

    elif cost.cost_type == 'NONLINEAR_LS':
        ny = cost.W.shape[0]
        if is_empty(model.cost_y_expr) and ny != 0:
            raise Exception('inconsistent dimension ny: regarding W, cost_y_expr.')
        elif casadi_length(model.cost_y_expr) != ny:
            raise Exception('inconsistent dimension ny: regarding W, cost_y_expr.')
        if cost.yref.shape[0] != ny:
            raise Exception('inconsistent dimension: regarding W, yref.' + \
                            f'\nGot W[{cost.W.shape}], yref[{cost.yref.shape}]\n')
        dims.ny = ny

    # terminal
    if cost.cost_type_e == 'LINEAR_LS':
        ny_e = cost.W_e.shape[0]
        if cost.Vx_e.shape[0] != ny_e:
            raise Exception('inconsistent dimension ny_e: regarding W_e, cost_y_expr_e.'  + \
                f'\nGot W_e[{cost.W_e.shape}], Vx_e[{cost.Vx_e.shape}]')
        if cost.Vx_e.shape[1] != dims.nx and ny_e != 0:
            raise Exception('inconsistent dimension: Vx_e should have nx columns.')
        if cost.yref_e.shape[0] != ny_e:
            raise Exception('inconsistent dimension: regarding W_e, yref_e.')
        dims.ny_e = ny_e

    elif cost.cost_type_e == 'NONLINEAR_LS':
        ny_e = cost.W_e.shape[0]
        if is_empty(model.cost_y_expr_e) and ny_e != 0:
            raise Exception('inconsistent dimension ny_e: regarding W_e, cost_y_expr_e.')
        elif casadi_length(model.cost_y_expr_e) != ny_e:
            raise Exception('inconsistent dimension ny_e: regarding W_e, cost_y_expr_e.')
        if cost.yref_e.shape[0] != ny_e:
            raise Exception('inconsistent dimension: regarding W_e, yref_e.')
        dims.ny_e = ny_e


    ## constraints
    # initial
    if (constraints.lbx_0 == [] and constraints.ubx_0 == []):
        dims.nbx_0 = 0
    else:
        this_shape = constraints.lbx_0.shape
        other_shape = constraints.ubx_0.shape
        if not this_shape == other_shape:
            raise Exception('lbx_0, ubx_0 have different shapes!')
        if not is_column(constraints.lbx_0):
            raise Exception('lbx_0, ubx_0 must be column vectors!')
        dims.nbx_0 = constraints.lbx_0.size

    if all(constraints.lbx_0 == constraints.ubx_0):
        dims.nbxe_0 = dims.nbx_0

    # path
    nbx = constraints.idxbx.shape[0]
    if constraints.ubx.shape[0] != nbx or constraints.lbx.shape[0] != nbx:
        raise Exception('inconsistent dimension nbx, regarding idxbx, ubx, lbx.')
    else:
        dims.nbx = nbx

    nbu = constraints.idxbu.shape[0]
    if constraints.ubu.shape[0] != nbu or constraints.lbu.shape[0] != nbu:
        raise Exception('inconsistent dimension nbu, regarding idxbu, ubu, lbu.')
    else:
        dims.nbu = nbu

    ng = constraints.lg.shape[0]
    if constraints.ug.shape[0] != ng or constraints.C.shape[0] != ng \
       or constraints.D.shape[0] != ng:
        raise Exception('inconsistent dimension ng, regarding lg, ug, C, D.')
    else:
        dims.ng = ng

    if not is_empty(model.con_h_expr):
        nh = casadi_length(model.con_h_expr)
    else:
        nh = 0

    if constraints.uh.shape[0] != nh or constraints.lh.shape[0] != nh:
        raise Exception('inconsistent dimension nh, regarding lh, uh, con_h_expr.')
    else:
        dims.nh = nh

    if is_empty(model.con_phi_expr):
        dims.nphi = 0
        dims.nr = 0
    else:
        dims.nphi = casadi_length(model.con_phi_expr)
        if is_empty(model.con_r_expr):
            raise Exception('convex over nonlinear constraints: con_r_expr but con_phi_expr is nonempty')
        else:
            dims.nr = casadi_length(model.con_r_expr)

    # terminal
    nbx_e = constraints.idxbx_e.shape[0]
    if constraints.ubx_e.shape[0] != nbx_e or constraints.lbx_e.shape[0] != nbx_e:
        raise Exception('inconsistent dimension nbx_e, regarding idxbx_e, ubx_e, lbx_e.')
    else:
        dims.nbx_e = nbx_e

    ng_e = constraints.lg_e.shape[0]
    if constraints.ug_e.shape[0] != ng_e or constraints.C_e.shape[0] != ng_e:
        raise Exception('inconsistent dimension ng_e, regarding_e lg_e, ug_e, C_e.')
    else:
        dims.ng_e = ng_e

    if not is_empty(model.con_h_expr_e):
        nh_e = casadi_length(model.con_h_expr_e)
    else:
        nh_e = 0

    if constraints.uh_e.shape[0] != nh_e or constraints.lh_e.shape[0] != nh_e:
        raise Exception('inconsistent dimension nh_e, regarding lh_e, uh_e, con_h_expr_e.')
    else:
        dims.nh_e = nh_e

    if is_empty(model.con_phi_expr_e):
        dims.nphi_e = 0
        dims.nr_e = 0
    else:
        dims.nphi_e = casadi_length(model.con_phi_expr_e)
        if is_empty(model.con_r_expr_e):
            raise Exception('convex over nonlinear constraints: con_r_expr_e but con_phi_expr_e is nonempty')
        else:
            dims.nr_e = casadi_length(model.con_r_expr_e)

    # Slack dimensions
    nsbx = constraints.idxsbx.shape[0]
    if is_empty(constraints.lsbx):
        constraints.lsbx = np.zeros((nsbx,))
    elif constraints.lsbx.shape[0] != nsbx:
        raise Exception('inconsistent dimension nsbx, regarding idxsbx, lsbx.')
    if is_empty(constraints.usbx):
        constraints.usbx = np.zeros((nsbx,))
    elif constraints.usbx.shape[0] != nsbx:
        raise Exception('inconsistent dimension nsbx, regarding idxsbx, usbx.')
    dims.nsbx = nsbx

    nsbu = constraints.idxsbu.shape[0]
    if is_empty(constraints.lsbu):
        constraints.lsbu = np.zeros((nsbu,))
    elif constraints.lsbu.shape[0] != nsbu:
        raise Exception('inconsistent dimension nsbu, regarding idxsbu, lsbu.')
    if is_empty(constraints.usbu):
        constraints.usbu = np.zeros((nsbu,))
    elif constraints.usbu.shape[0] != nsbu:
        raise Exception('inconsistent dimension nsbu, regarding idxsbu, usbu.')
    dims.nsbu = nsbu

    nsh = constraints.idxsh.shape[0]
    if is_empty(constraints.lsh):
        constraints.lsh = np.zeros((nsh,))
    elif constraints.lsh.shape[0] != nsh:
        raise Exception('inconsistent dimension nsh, regarding idxsh, lsh.')
    if is_empty(constraints.ush):
        constraints.ush = np.zeros((nsh,))
    elif constraints.ush.shape[0] != nsh:
        raise Exception('inconsistent dimension nsh, regarding idxsh, ush.')
    dims.nsh = nsh

    nsphi = constraints.idxsphi.shape[0]
    if is_empty(constraints.lsphi):
        constraints.lsphi = np.zeros((nsphi,))
    elif constraints.lsphi.shape[0] != nsphi:
        raise Exception('inconsistent dimension nsphi, regarding idxsphi, lsphi.')
    if is_empty(constraints.usphi):
        constraints.usphi = np.zeros((nsphi,))
    elif constraints.usphi.shape[0] != nsphi:
        raise Exception('inconsistent dimension nsphi, regarding idxsphi, usphi.')
    dims.nsphi = nsphi

    nsg = constraints.idxsg.shape[0]
    if is_empty(constraints.lsg):
        constraints.lsg = np.zeros((nsg,))
    elif constraints.lsg.shape[0] != nsg:
        raise Exception('inconsistent dimension nsg, regarding idxsg, lsg.')
    if is_empty(constraints.usg):
        constraints.usg = np.zeros((nsg,))
    elif constraints.usg.shape[0] != nsg:
        raise Exception('inconsistent dimension nsg, regarding idxsg, usg.')
    dims.nsg = nsg

    ns = nsbx + nsbu + nsh + nsg + nsphi
    wrong_field = ""
    if cost.Zl.shape[0] != ns:
        wrong_field = "Zl"
        dim = cost.Zl.shape[0]
    elif cost.Zu.shape[0] != ns:
        wrong_field = "Zu"
        dim = cost.Zu.shape[0]
    elif cost.zl.shape[0] != ns:
        wrong_field = "zl"
        dim = cost.zl.shape[0]
    elif cost.zu.shape[0] != ns:
        wrong_field = "zu"
        dim = cost.zu.shape[0]

    if wrong_field != "":
        raise Exception(f'Inconsistent size for field {wrong_field}, with dimension {dim}, \n\t'\
            + f'Detected ns = {ns} = nsbx + nsbu + nsg + nsh + nsphi.\n\t'\
            + f'With nsbx = {nsbx}, nsbu = {nsbu}, nsg = {nsg}, nsh = {nsh}, nsphi = {nsphi}')

    dims.ns = ns

    nsbx_e = constraints.idxsbx_e.shape[0]
    if is_empty(constraints.lsbx_e):
        constraints.lsbx_e = np.zeros((nsbx_e,))
    elif constraints.lsbx_e.shape[0] != nsbx_e:
        raise Exception('inconsistent dimension nsbx_e, regarding idxsbx_e, lsbx_e.')
    if is_empty(constraints.usbx_e):
        constraints.usbx_e = np.zeros((nsbx_e,))
    elif constraints.usbx_e.shape[0] != nsbx_e:
        raise Exception('inconsistent dimension nsbx_e, regarding idxsbx_e, usbx_e.')
    dims.nsbx_e = nsbx_e

    nsh_e = constraints.idxsh_e.shape[0]
    if is_empty(constraints.lsh_e):
        constraints.lsh_e = np.zeros((nsh_e,))
    elif constraints.lsh_e.shape[0] != nsh_e:
        raise Exception('inconsistent dimension nsh_e, regarding idxsh_e, lsh_e.')
    if is_empty(constraints.ush_e):
        constraints.ush_e = np.zeros((nsh_e,))
    elif constraints.ush_e.shape[0] != nsh_e:
        raise Exception('inconsistent dimension nsh_e, regarding idxsh_e, ush_e.')
    dims.nsh_e = nsh_e

    nsg_e = constraints.idxsg_e.shape[0]
    if is_empty(constraints.lsg_e):
        constraints.lsg_e = np.zeros((nsg_e,))
    elif constraints.lsg_e.shape[0] != nsg_e:
        raise Exception('inconsistent dimension nsg_e, regarding idxsg_e, lsg_e.')
    if is_empty(constraints.usg_e):
        constraints.usg_e = np.zeros((nsg_e,))
    elif constraints.usg_e.shape[0] != nsg_e:
        raise Exception('inconsistent dimension nsg_e, regarding idxsg_e, usg_e.')
    dims.nsg_e = nsg_e

    nsphi_e = constraints.idxsphi_e.shape[0]
    if is_empty(constraints.lsphi_e):
        constraints.lsphi_e = np.zeros((nsphi_e,))
    elif constraints.lsphi_e.shape[0] != nsphi_e:
        raise Exception('inconsistent dimension nsphi_e, regarding idxsphi_e, lsphi_e.')
    if is_empty(constraints.usphi_e):
        constraints.usphi_e = np.zeros((nsphi_e,))
    elif constraints.usphi_e.shape[0] != nsphi_e:
        raise Exception('inconsistent dimension nsphi_e, regarding idxsphi_e, usphi_e.')
    dims.nsphi_e = nsphi_e

    # terminal
    ns_e = nsbx_e + nsh_e + nsg_e + nsphi_e
    wrong_field = ""
    if cost.Zl_e.shape[0] != ns_e:
        wrong_field = "Zl_e"
        dim = cost.Zl_e.shape[0]
    elif cost.Zu_e.shape[0] != ns_e:
        wrong_field = "Zu_e"
        dim = cost.Zu_e.shape[0]
    elif cost.zl_e.shape[0] != ns_e:
        wrong_field = "zl_e"
        dim = cost.zl_e.shape[0]
    elif cost.zu_e.shape[0] != ns_e:
        wrong_field = "zu_e"
        dim = cost.zu_e.shape[0]

    if wrong_field != "":
        raise Exception(f'Inconsistent size for field {wrong_field}, with dimension {dim}, \n\t'\
            + f'Detected ns_e = {ns_e} = nsbx_e + nsg_e + nsh_e + nsphi_e.\n\t'\
            + f'With nsbx_e = {nsbx_e}, nsg_e = {nsg_e}, nsh_e = {nsh_e}, nsphi_e = {nsphi_e}')

    dims.ns_e = ns_e

    # discretization
    if is_empty(opts.time_steps) and is_empty(opts.shooting_nodes):
        # uniform discretization
        opts.time_steps = opts.tf / dims.N * np.ones((dims.N,))

    elif not is_empty(opts.shooting_nodes):
        if np.shape(opts.shooting_nodes)[0] != dims.N+1:
            raise Exception('inconsistent dimension N, regarding shooting_nodes.')

        time_steps = np.zeros((dims.N,))
        for i in range(dims.N):
            time_steps[i] = opts.shooting_nodes[i+1] - opts.shooting_nodes[i]
        opts.time_steps = time_steps

    elif (not is_empty(opts.time_steps)) and (not is_empty(opts.shooting_nodes)):
        Exception('Please provide either time_steps or shooting_nodes for nonuniform discretization')

    tf = np.sum(opts.time_steps)
    if (tf - opts.tf) / tf > 1e-15:
        raise Exception(f'Inconsistent discretization: {opts.tf}'\
            f' = tf != sum(opts.time_steps) = {tf}.')



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
    # TODO: maybe make one funciton with formatting

    for acados_struct, v in ocp_layout.items():
        # skip non dict attributes
        if not isinstance(v, dict): continue
        #  setattr(ocp_nlp, acados_struct, dict(getattr(acados_ocp, acados_struct).__dict__))
        # Copy ocp object attributes dictionaries
        ocp_nlp_dict[acados_struct]=dict(getattr(acados_ocp, acados_struct).__dict__)

    ocp_nlp_dict = format_class_dict(ocp_nlp_dict)

    # strip symbolics
    ocp_nlp_dict['model'] = acados_model_strip_casadi_symbolics(ocp_nlp_dict['model'])

    # strip shooting_nodes
    ocp_nlp_dict['solver_options'].pop('shooting_nodes', None)

    dims_dict = acados_class2dict(acados_ocp.dims)

    ocp_check_against_layout(ocp_nlp_dict, dims_dict)

    with open(json_file, 'w') as f:
        json.dump(ocp_nlp_dict, f, default=np_array_to_list, indent=4, sort_keys=True)



def ocp_formulation_json_load(json_file='acados_ocp_nlp.json'):
    # Load acados_ocp_nlp structure description
    ocp_layout = get_ocp_nlp_layout()

    with open(json_file, 'r') as f:
        ocp_nlp_json = json.load(f)

    ocp_nlp_dict = json2dict(ocp_nlp_json, ocp_nlp_json['dims'])

    # Instantiate AcadosOcp object
    acados_ocp = AcadosOcp()

    # load class dict
    acados_ocp.__dict__ = ocp_nlp_dict

    # laod class attributes dict, dims, constraints, etc
    for acados_struct, v  in ocp_layout.items():
        # skip non dict attributes
        if not isinstance(v, dict): continue
        acados_attribute = getattr(acados_ocp, acados_struct)
        acados_attribute.__dict__ = ocp_nlp_dict[acados_struct]
        setattr(acados_ocp, acados_struct, acados_attribute)

    return acados_ocp


def ocp_generate_external_functions(acados_ocp, model):

    model = make_model_consistent(model)
    if acados_ocp.solver_options.integrator_type == 'ERK':
        # explicit model -- generate C code
        generate_c_code_explicit_ode(model)
    elif acados_ocp.solver_options.integrator_type == 'IRK':
        # implicit model -- generate C code
        opts = dict(generate_hess=1)
        generate_c_code_implicit_ode(model, opts)
    elif acados_ocp.solver_options.integrator_type == 'GNSF':
        generate_c_code_gnsf(model)
    else:
        raise Exception("ocp_generate_external_functions: unknown integrator type.")


    if acados_ocp.solver_options.hessian_approx == 'EXACT':
        opts = dict(generate_hess=1)
    else:
        opts = dict(generate_hess=0)

    if acados_ocp.dims.nphi > 0 or acados_ocp.dims.nh > 0:
        generate_c_code_constraint(model, model.name, False, opts)

    if acados_ocp.dims.nphi_e > 0 or acados_ocp.dims.nh_e > 0:
        generate_c_code_constraint(model, model.name, True, opts)

    # dummy matrices
    if not acados_ocp.cost.cost_type == 'LINEAR_LS':
        acados_ocp.cost.Vx = np.zeros((acados_ocp.dims.ny, acados_ocp.dims.nx))
        acados_ocp.cost.Vu = np.zeros((acados_ocp.dims.ny, acados_ocp.dims.nu))
    if not acados_ocp.cost.cost_type_e == 'LINEAR_LS':
        acados_ocp.cost.Vx_e = np.zeros((acados_ocp.dims.ny_e, acados_ocp.dims.nx))


    if acados_ocp.cost.cost_type == 'NONLINEAR_LS':
        generate_c_code_nls_cost(model, model.name, False)
    elif acados_ocp.cost.cost_type == 'EXTERNAL':
        generate_c_code_external_cost(model, False)

    if acados_ocp.cost.cost_type_e == 'NONLINEAR_LS':
        generate_c_code_nls_cost(model, model.name, True)
    elif acados_ocp.cost.cost_type_e == 'EXTERNAL':
        generate_c_code_external_cost(model, True)


def ocp_render_templates(acados_ocp, json_file):

    name = acados_ocp.model.name

    # setting up loader and environment
    json_path = '{cwd}/{json_file}'.format(
        cwd=os.getcwd(),
        json_file=json_file)

    if not os.path.exists(json_path):
        raise Exception('{} not found!'.format(json_path))

    template_dir = 'c_generated_code/'

    ## Render templates
    in_file = 'main.in.c'
    out_file = 'main_{}.c'.format(name)
    render_template(in_file, out_file, template_dir, json_path)

    in_file = 'acados_solver.in.c'
    out_file = 'acados_solver_{}.c'.format(name)
    render_template(in_file, out_file, template_dir, json_path)

    in_file = 'acados_solver.in.h'
    out_file = 'acados_solver_{}.h'.format(name)
    render_template(in_file, out_file, template_dir, json_path)

    in_file = 'Makefile.in'
    out_file = 'Makefile'
    render_template(in_file, out_file, template_dir, json_path)

    in_file = 'acados_solver_sfun.in.c'
    out_file = 'acados_solver_sfunction_{}.c'.format(name)
    render_template(in_file, out_file, template_dir, json_path)

    in_file = 'make_sfun.in.m'
    out_file = 'make_sfun.m'
    render_template(in_file, out_file, template_dir, json_path)

    in_file = 'acados_sim_solver.in.c'
    out_file = 'acados_sim_solver_{}.c'.format(name)
    render_template(in_file, out_file, template_dir, json_path)

    in_file = 'acados_sim_solver.in.h'
    out_file = 'acados_sim_solver_{}.h'.format(name)
    render_template(in_file, out_file, template_dir, json_path)

    ## folder model
    template_dir = 'c_generated_code/{}_model/'.format(name)

    in_file = 'model.in.h'
    out_file = '{}_model.h'.format(name)
    render_template(in_file, out_file, template_dir, json_path)

    # constraints on convex over nonlinear function
    if acados_ocp.constraints.constr_type == 'BGP' and acados_ocp.dims.nphi > 0:
        # constraints on outer function
        template_dir = 'c_generated_code/{}_constraints/'.format(name)
        in_file = 'phi_constraint.in.h'
        out_file =  '{}_phi_constraint.h'.format(name)
        render_template(in_file, out_file, template_dir, json_path)

    # terminal constraints on convex over nonlinear function
    if acados_ocp.constraints.constr_type_e == 'BGP' and acados_ocp.dims.nphi_e > 0:
        # terminal constraints on outer function
        template_dir = 'c_generated_code/{}_constraints/'.format(name)
        in_file = 'phi_e_constraint.in.h'
        out_file =  '{}_phi_e_constraint.h'.format(name)
        render_template(in_file, out_file, template_dir, json_path)

    # nonlinear constraints
    if acados_ocp.constraints.constr_type == 'BGH' and acados_ocp.dims.nh > 0:
        template_dir = 'c_generated_code/{}_constraints/'.format(name)
        in_file = 'h_constraint.in.h'
        out_file = '{}_h_constraint.h'.format(name)
        render_template(in_file, out_file, template_dir, json_path)

    # terminal nonlinear constraints
    if acados_ocp.constraints.constr_type_e == 'BGH' and acados_ocp.dims.nh_e > 0:
        template_dir = 'c_generated_code/{}_constraints/'.format(name)
        in_file = 'h_e_constraint.in.h'
        out_file = '{}_h_e_constraint.h'.format(name)
        render_template(in_file, out_file, template_dir, json_path)

    # nonlinear cost function
    if acados_ocp.cost.cost_type == 'NONLINEAR_LS':
        template_dir = 'c_generated_code/{}_cost/'.format(name)
        in_file = 'cost_y_fun.in.h'
        out_file = '{}_cost_y_fun.h'.format(name)
        render_template(in_file, out_file, template_dir, json_path)

    # terminal nonlinear cost function
    if acados_ocp.cost.cost_type_e == 'NONLINEAR_LS':
        template_dir = 'c_generated_code/{}_cost/'.format(name)
        in_file = 'cost_y_e_fun.in.h'
        out_file = '{}_cost_y_e_fun.h'.format(name)
        render_template(in_file, out_file, template_dir, json_path)

    # external cost
    if acados_ocp.cost.cost_type == 'EXTERNAL':
        template_dir = 'c_generated_code/{}_cost/'.format(name)
        in_file = 'external_cost.in.h'
        out_file = '{}_external_cost.h'.format(name)
        render_template(in_file, out_file, template_dir, json_path)

    # external cost - terminal
    if acados_ocp.cost.cost_type_e == 'EXTERNAL':
        template_dir = 'c_generated_code/{}_cost/'.format(name)
        in_file = 'external_cost_e.in.h'
        out_file = '{}_external_cost_e.h'.format(name)
        render_template(in_file, out_file, template_dir, json_path)




class AcadosOcpSolver:
    """
    class to interact with the acados ocp solver C object
    """
    def __init__(self, acados_ocp, json_file='acados_ocp_nlp.json'):

        self.solver_created = False
        self.N = acados_ocp.dims.N
        model = acados_ocp.model

        # make dims consistent
        make_ocp_dims_consistent(acados_ocp)

        if acados_ocp.solver_options.integrator_type == 'GNSF':
            set_up_imported_gnsf_model(acados_ocp)

        # set integrator time automatically
        acados_ocp.solver_options.Tsim = acados_ocp.solver_options.time_steps[0]

        # generate external functions
        ocp_generate_external_functions(acados_ocp, model)

        # dump to json
        ocp_formulation_json_dump(acados_ocp, json_file)

        # render templates
        ocp_render_templates(acados_ocp, json_file)

        ## Compile solver
        os.chdir('c_generated_code')
        os.system('make clean_ocp_shared_lib')
        os.system('make ocp_shared_lib')
        os.chdir('..')

        self.shared_lib_name = 'c_generated_code/libacados_ocp_solver_' + model.name + '.so'

        # get
        self.shared_lib = CDLL(self.shared_lib_name)
        self.shared_lib.acados_create()
        self.solver_created = True

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
        """
        solve the ocp with current input
        """
        status = self.shared_lib.acados_solve()
        return status


    def get(self, stage_, field_):
        """
        get the last solution of the solver:
            :param stage: integer corresponding to shooting node
            :param field_: string in ['x', 'u', 'z', 'pi', 'lam', 't', 'sl', 'su',]

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

        out_fields = ['x', 'u', 'z', 'pi', 'lam', 't']
        mem_fields = ['sl', 'su']
        field = field_
        field = field.encode('utf-8')

        if (field_ not in out_fields + mem_fields):
            raise Exception('AcadosOcpSolver.get(): {} is an invalid argument.\
                    \n Possible values are {}. Exiting.'.format(field_, out_fields + mem_fields))

        if not isinstance(stage_, int):
            raise Exception('AcadosOcpSolver.get(): stage index must be Integer.')

        if stage_ < 0 or stage_ > self.N:
            raise Exception('AcadosOcpSolver.get(): stage index must be in [0, N], got: {}.'.format(self.N))

        if stage_ == self.N and field_ == 'pi':
            raise Exception('AcadosOcpSolver.get(): field {} does not exist at final stage {}.'\
                .format(field_, stage_))

        self.shared_lib.ocp_nlp_dims_get_from_attr.argtypes = \
            [c_void_p, c_void_p, c_void_p, c_int, c_char_p]
        self.shared_lib.ocp_nlp_dims_get_from_attr.restype = c_int

        dims = self.shared_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, \
            self.nlp_dims, self.nlp_out, stage_, field)

        out = np.ascontiguousarray(np.zeros((dims,)), dtype=np.float64)
        out_data = cast(out.ctypes.data, POINTER(c_double))

        if (field_ in out_fields):
            self.shared_lib.ocp_nlp_out_get.argtypes = \
                [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]
            self.shared_lib.ocp_nlp_out_get(self.nlp_config, \
                self.nlp_dims, self.nlp_out, stage_, field, out_data)
        elif field_ in mem_fields:
            self.shared_lib.ocp_nlp_get_at_stage.argtypes = \
                [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]
            self.shared_lib.ocp_nlp_get_at_stage(self.nlp_config, \
                self.nlp_dims, self.nlp_solver, stage_, field, out_data)

        return out


    def print_statistics(self):
        stat = self.get_stats("statistics")

        if self.acados_ocp.solver_options.nlp_solver_type == 'SQP':
            print('\niter\tres_stat\tres_eq\t\tres_ineq\tres_comp\tqp_stat\tqp_iter')
            if stat.shape[0]>7:
                print('\tqp_res_stat\tqp_res_eq\tqp_res_ineq\tqp_res_comp')
            for jj in range(stat.shape[1]):
                print('{:d}\t{:e}\t{:e}\t{:e}\t{:e}\t{:d}\t{:d}'.format( \
                     int(stat[0][jj]), stat[1][jj], stat[2][jj], \
                     stat[3][jj], stat[4][jj], int(stat[5][jj]), int(stat[6][jj])))
                if stat.shape[0]>7:
                    print('\t{:e}\t{:e}\t{:e}\t{:e}'.format( \
                        stat[7][jj], stat[8][jj], stat[9][jj], stat[10][jj]))
            print('\n')
        elif self.acados_ocp.solver_options.nlp_solver_type == 'SQP_RTI':
            print('\niter\tqp_stat\tqp_iter')
            if stat.shape[0]>3:
                print('\tqp_res_stat\tqp_res_eq\tqp_res_ineq\tqp_res_comp')
            for jj in range(stat.shape[1]):
                print('{:d}\t{:d}\t{:d}'.format( int(stat[0][jj]), int(stat[1][jj]), int(stat[2][jj])))
                if stat.shape[0]>3:
                    print('\t{:e}\t{:e}\t{:e}\t{:e}'.format( \
                         stat[3][jj], stat[4][jj], stat[5][jj], stat[6][jj]))
            print('\n')

        return


    def get_stats(self, field_):
        """
        get the information of the last solver call:
            :param field_: string in ['statistics', 'time_tot', 'time_lin', 'time_sim', 'time_sim_ad', 'time_sim_la', 'time_qp', 'time_qp_solver_call', 'time_reg', 'sqp_iter']
        """

        fields = ['time_tot',  # total cpu time previous call
                  'time_lin',  # cpu time for linearization
                  'time_sim',  # cpu time for integrator
                  'time_sim_ad',  # cpu time for integrator contribution of external function calls
                  'time_sim_la',  # cpu time for integrator contribution of linear algebra
                  'time_qp',   # cpu time qp solution
                  'time_qp_solver_call',  # cpu time inside qp solver (without converting the QP)
                  'time_qp_xcond',
                  'time_reg',  # cpu time regularization
                  'sqp_iter',  # number of SQP iterations
                  'statistics',  # table with info about last iteration
                  'stat_m',
                  'stat_n',
                ]

        field = field_
        field = field.encode('utf-8')
        if (field_ not in fields):
            raise Exception('AcadosOcpSolver.get_stats(): {} is not a valid argument.\
                    \n Possible values are {}. Exiting.'.format(fields, fields))

        if field_ in ['sqp_iter', 'stat_m', 'stat_n']:
            out = np.ascontiguousarray(np.zeros((1,)), dtype=np.int64)
            out_data = cast(out.ctypes.data, POINTER(c_int64))

        elif field_ == 'statistics':
            sqp_iter = self.get_stats("sqp_iter")
            stat_m = self.get_stats("stat_m")
            stat_n = self.get_stats("stat_n")

            min_size = min([stat_m, sqp_iter+1])

            out = np.ascontiguousarray(
                        np.zeros( (stat_n[0]+1, min_size[0]) ), dtype=np.float64)
            out_data = cast(out.ctypes.data, POINTER(c_double))

        else:
            out = np.ascontiguousarray(np.zeros((1,)), dtype=np.float64)
            out_data = cast(out.ctypes.data, POINTER(c_double))

        self.shared_lib.ocp_nlp_get.argtypes = [c_void_p, c_void_p, c_char_p, c_void_p]
        self.shared_lib.ocp_nlp_get(self.nlp_config, self.nlp_solver, field, out_data)

        return out

    # Note: this function should not be used anymore, better use cost_set, constraints_set
    def set(self, stage_, field_, value_):
        """
        set numerical data inside the solver:
            :param stage_: integer corresponding to shooting node
            :param field_: string in ['x', 'u', 'pi', 'lam', 't']

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
        out_fields = ['x', 'u', 'pi', 'lam', 't']

        # cast value_ to avoid conversion issues
        value_ = value_.astype(float)

        field = field_
        field = field.encode('utf-8')

        stage = c_int(stage_)

        # treat parameters separately
        if field_ is 'p':
            self.shared_lib.acados_update_params.argtypes = [c_int, POINTER(c_double)]
            self.shared_lib.acados_update_params.restype = c_int
            value_data = cast(value_.ctypes.data, POINTER(c_double))
            self.shared_lib.acados_update_params(stage, value_data, value_.shape[0])
        else:
            if field_ not in constraints_fields + cost_fields + out_fields:
                raise Exception("AcadosOcpSolver.set(): {} is not a valid argument.\
                    \nPossible values are {}. Exiting.".format(field, \
                    constraints_fields + cost_fields + out_fields + ['p']))

            self.shared_lib.ocp_nlp_dims_get_from_attr.argtypes = \
                [c_void_p, c_void_p, c_void_p, c_int, c_char_p]
            self.shared_lib.ocp_nlp_dims_get_from_attr.restype = c_int

            dims = self.shared_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, \
                self.nlp_dims, self.nlp_out, stage_, field)

            if value_.shape[0] != dims:
                msg = 'AcadosOcpSolver.set(): mismatching dimension for field "{}" '.format(field_)
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
        """
        set numerical data in the cost module of the solver:
            :param stage_: integer corresponding to shooting node
            :param field_: string, e.g. 'yref', 'W', 'ext_cost_num_hess'
            :param value_: of appropriate size
        """
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
            raise Exception('AcadosOcpSolver.cost_set(): mismatching dimension', \
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
        """
        set numerical data in the constraint module of the solver:
        Parameters:
            :param stage_: integer corresponding to shooting node
            :param field_: string, e.g. 'lbx'
            :param value_: of appropriate size
        """
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
            raise Exception('AcadosOcpSolver.constraints_set(): mismatching dimension' \
                ' for field "{}" with dimension {} (you have {})'.format(field_, tuple(dims), value_shape))

        value_data = cast(value_.ctypes.data, POINTER(c_double))
        value_data_p = cast((value_data), c_void_p)

        self.shared_lib.ocp_nlp_constraints_model_set.argtypes = \
            [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]
        self.shared_lib.ocp_nlp_constraints_model_set(self.nlp_config, \
            self.nlp_dims, self.nlp_in, stage, field, value_data_p)

        return


    def options_set(self, field_, value_):
        """
        set options of the solver:
        Parameters:
            :param field_: string, e.g. 'print_level', 'rti_phase', 'initialize_t_slacks', 'step_length'
            :param value_: of type int, float
        """
        int_fields = ['print_level', 'rti_phase', 'initialize_t_slacks']
        double_fields = ['step_length']
        string_fields = ['globalization']

        if field_ in int_fields:
            if not isinstance(value_, int):
                raise Exception('solver option {} must be of type int. You have {}.'.format(field_, type(value_)))
            else:
                value_ctypes = c_int(value_)

        elif field_ in double_fields:
            if not isinstance(value_, float):
                raise Exception('solver option {} must be of type float. You have {}.'.format(field_, type(value_)))
            else:
                value_ctypes = c_double(value_)

        elif field_ in string_fields:
            if not isinstance(value_, str):
                raise Exception('solver option {} must be of type str. You have {}.'.format(field_, type(value_)))
            else:
                value_ctypes = value_.encode('utf-8')

        if field_ == 'rti_phase':
            if value_ < 0 or value_ > 2:
                raise Exception('AcadosOcpSolver.solve(): argument \'rti_phase\' can '
                    'take only values 0, 1, 2 for SQP-RTI-type solvers')
            if self.acados_ocp.solver_options.nlp_solver_type != 'SQP_RTI' and value_ > 0:
                raise Exception('AcadosOcpSolver.solve(): argument \'rti_phase\' can '
                    'take only value 0 for SQP-type solvers')

        field = field_
        field = field.encode('utf-8')

        if field_ in string_fields:
            self.shared_lib.ocp_nlp_solver_opts_set.argtypes = \
                [c_void_p, c_void_p, c_char_p, c_char_p]
            self.shared_lib.ocp_nlp_solver_opts_set(self.nlp_config, \
                self.nlp_opts, field, value_ctypes)
        else:
            self.shared_lib.ocp_nlp_solver_opts_set.argtypes = \
                [c_void_p, c_void_p, c_char_p, c_void_p]
            self.shared_lib.ocp_nlp_solver_opts_set(self.nlp_config, \
                self.nlp_opts, field, byref(value_ctypes))

        return


    def __del__(self):
        if self.solver_created:
            self.shared_lib.acados_free()
            del self.shared_lib

        # NOTE: DLL cannot be easily unloaded!!!
        # see https://stackoverflow.com/questions/359498/how-can-i-unload-a-dll-using-ctypes-in-python
        # while isLoaded(self.shared_lib_name):
        #     dlclose(handle)
