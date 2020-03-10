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
from casadi import CasadiMeta, Function

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
     format_class_dict, ocp_check_json_against_layout, np_array_to_list, make_model_consistent


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
            raise Exception('lbx_0, ubx_0 have different shapes!')
        if not is_column(constraints.lbx_0):
            raise Exception('lbx_0, ubx_0 must be column vectors!')

        dims.nbx_0 = constraints.lbx_0.size
    else:
        raise Exception('lbx_0, ubx_0 have different shapes!')

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


def set_up_imported_gnsf_model(acados_ocp):

    gnsf = acados_ocp.gnsf_model

    # check CasADi version
    dump_casadi_version = gnsf['casadi_version']
    casadi_version = CasadiMeta.version()

    if not casadi_version == dump_casadi_version:
        raise Exception("GNSF model was dumped with another CasADi version.\n"
                + "Please use the same version for compatibility, serialize version:"
                + dump_casadi_version + " current Python CasADi verison: " + casadi_version)

    # load model
    phi_fun = Function.deserialize(gnsf['phi_fun'])
    phi_fun_jac_y = Function.deserialize(gnsf['phi_fun_jac_y'])
    phi_jac_y_uhat = Function.deserialize(gnsf['phi_jac_y_uhat'])
    f_lo_fun_jac_x1k1uz = Function.deserialize(gnsf['f_lo_fun_jac_x1k1uz'])
    get_matrices_fun = Function.deserialize(gnsf['get_matrices_fun'])

    # obtain gnsf dimensions
    size_gnsf_A = get_matrices_fun.size_out(0)
    acados_ocp.dims.gnsf_nx1 = size_gnsf_A[1]
    acados_ocp.dims.gnsf_nz1 = size_gnsf_A[0] - size_gnsf_A[1]
    acados_ocp.dims.gnsf_nuhat = max(phi_fun.size_in(1))
    acados_ocp.dims.gnsf_ny = max(phi_fun.size_in(0))
    acados_ocp.dims.gnsf_nout = max(phi_fun.size_out(0))

    # save gnsf functions in model
    acados_ocp.model.phi_fun = phi_fun
    acados_ocp.model.phi_fun_jac_y = phi_fun_jac_y
    acados_ocp.model.phi_jac_y_uhat = phi_jac_y_uhat
    acados_ocp.model.f_lo_fun_jac_x1k1uz = f_lo_fun_jac_x1k1uz
    acados_ocp.model.get_matrices_fun = get_matrices_fun

    del acados_ocp.gnsf_model


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

    # strip symbolics
    ocp_nlp_dict['model'] = acados_model_strip_casadi_symbolics(ocp_nlp_dict['model'])

    ocp_nlp_dict = format_class_dict(ocp_nlp_dict)

    dims_dict = acados_class2dict(acados_ocp.dims)

    ocp_check_json_against_layout(ocp_nlp_dict, dims_dict)

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


    if acados_ocp.dims.nphi > 0 or acados_ocp.dims.nh > 0:
        generate_c_code_constraint(model, model.name, False)

    if acados_ocp.dims.nphi_e > 0 or acados_ocp.dims.nh_e > 0:
        generate_c_code_constraint(model, model.name, True)

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

        model = acados_ocp.model

        # make dims consistent
        make_ocp_dims_consistent(acados_ocp)

        if acados_ocp.solver_options.integrator_type == 'GNSF':
            set_up_imported_gnsf_model(acados_ocp)

        # set integrator time automatically
        acados_ocp.solver_options.Tsim = acados_ocp.solver_options.tf / acados_ocp.dims.N

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


    def solve(self, rti_phase=0):
        """
        solve the ocp with current input
        :param rti_phase: 0 = preparation + feedback, 1 = preparation only,
         2 = feedback only (if SQP_RTI is used, otherwise only 0 (default) is allowed)
        """
        if isinstance(rti_phase, int) == False or rti_phase < 0 or rti_phase > 2: 
            raise Exception('AcadosOcpSolver.solve(): argument \'rti_phase\' can ' 
                'take only values 0, 1, 2 for SQP-RTI-type solvers')
        if self.acados_ocp.solver_options.nlp_solver_type != 'SQP_RTI' and rti_phase > 0:
            raise Exception('AcadosOcpSolver.solve(): argument \'rti_phase\' can ' 
                'take only value 0 for SQP-type solvers')
        self.shared_lib.acados_solve.argtypes = [c_int]

        status = self.shared_lib.acados_solve(rti_phase)
        return status


    def get(self, stage_, field_):
        """
        get the last solution of the solver:
            :param stage: integer corresponding to shooting node
            :param field_: string in ['x', 'u', 'z', 'pi']
        """

        out_fields = ['x', 'u', 'z', 'pi']
        field = field_
        field = field.encode('utf-8')

        if (field_ not in out_fields):
            raise Exception('AcadosOcpSolver.get(): {} is an invalid argument.\
                    \n Possible values are {}. Exiting.'.format(out_fields))

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
            :param field_: string in ['time_tot', 'time_lin', 'time_qp', 'time_reg', 'sqp_iter']
        """

        fields = ['time_tot',  # total cpu time previous call
                  'time_lin',  # cpu time for linearization
                  'time_qp',   # cpu time qp solution
                  'time_qp_solver_call',  # cpu time inside qp solver (without converting the QP)
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

        cost_fields = ['y_ref', 'yref']
        constraints_fields = ['lbx', 'ubx', 'lbu', 'ubu']
        out_fields = ['x', 'u', 'pi']

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
            if (field_ not in constraints_fields) and \
                    (field_ not in cost_fields) and (field_ not in out_fields):
                raise Exception("AcadosOcpSolver.set(): {} is not a valid argument.\
                    \nPossible values are {} and {}. Exiting.".format(field, \
                    cost_fields, constraints_fields, out_fields))

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
            :param field_: string, e.g. 'yref', 'W'
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
            raise Exception('acados_solver.set(): mismatching dimension' \
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
            :param field_: string, e.g. 'print_level'
            :param value_: of type int
        """
        # cast value_ to avoid conversion issues
        if isinstance(value_, int) is not True:
            raise Exception('solver options must be of type int. You have {}'.format())

        field = field_
        field = field.encode('utf-8')

        value_ctypes = c_int(value_)

        self.shared_lib.ocp_nlp_solver_opts_set.argtypes = \
            [c_void_p, c_void_p, c_char_p, c_void_p]
        self.shared_lib.ocp_nlp_solver_opts_set(self.nlp_config, \
            self.nlp_opts, field, byref(value_ctypes))

        return


    def __del__(self):
        self.shared_lib.acados_free()
        del self.shared_lib

        # NOTE: DLL cannot be easily unloaded!!!
        # see https://stackoverflow.com/questions/359498/how-can-i-unload-a-dll-using-ctypes-in-python
        # while isLoaded(self.shared_lib_name):
        #     dlclose(handle)
