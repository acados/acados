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
from ctypes import *
from copy import deepcopy
from .utils import ACADOS_PATH, get_tera

def generate_solver(acados_ocp, json_file='acados_ocp_nlp.json'):

    # get tera renderer
    tera_path = get_tera()

    model = acados_ocp.model
    name = model.name
    if acados_ocp.solver_options.integrator_type == 'ERK':
        # explicit model -- generate C code
        generate_c_code_explicit_ode(model)
    else:
        # implicit model -- generate C code
        opts = dict(generate_hess=1)
        generate_c_code_implicit_ode(model, opts)

    if acados_ocp.constraints.constr_type == 'BGP' and acados_ocp.dims.nphi > 0:
        # nonlinear part of nonlinear constraints
        generate_c_code_constraint(acados_ocp.con_phi, name)
    elif acados_ocp.constraints.constr_type  == 'BGH' and acados_ocp.dims.nh > 0:
        generate_c_code_constraint(acados_ocp.con_h, name)

    if acados_ocp.constraints.constr_type_e  == 'BGP' and acados_ocp.dims.nphi_e > 0:
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

    ## Load solver
    solver = acados_solver(
        acados_ocp,
        'c_generated_code/libacados_solver_{model_name}.so'.format(model_name=model.name)
    )
    return solver

class acados_solver:
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
            self.nlp_dims, self.nlp_in, stage, field, value_data_p);

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
            self.nlp_dims, self.nlp_in, stage, field, value_data_p);


    def solve(self):
        status = self.shared_lib.acados_solve()
        return status

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
        self.shared_lib.ocp_nlp_get(self.nlp_config, self.nlp_solver, field, out_data);

        return out

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

    def __del__(self):
        self.shared_lib.acados_free()


