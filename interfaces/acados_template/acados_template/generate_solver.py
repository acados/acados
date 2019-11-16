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

from jinja2 import Environment, FileSystemLoader
from .generate_c_code_explicit_ode import *
from .generate_c_code_implicit_ode import *
from .generate_c_code_constraint import *
from .generate_c_code_constraint_e import *
from .generate_c_code_nls_cost import *
from .generate_c_code_nls_cost_e import *
from .acados_ocp_nlp import *
from ctypes import *
from copy import deepcopy

def generate_solver(acados_ocp, json_file='acados_ocp_nlp.json'):

    if os.environ.get('ACADOS_SOURCE_DIR') is None: 
        acados_template_path = os.path.dirname(os.path.abspath(__file__))
        acados_path = acados_template_path + '/../../../'
        tera_path = acados_path + '/bin/' 
    else:
        acados_path = os.environ.get('ACADOS_SOURCE_DIR')
        tera_path = acados_path + '/bin/'
    t_renderer_path = tera_path + 't_renderer'
    if (os.path.exists(t_renderer_path) is False):
        raise Exception('t_renderer binaries not found. In order to be able to ' + \
                'successfully render C code templates, you need to download the ' + \
                't_renderer binaries for your platform from ' + \
                'https://github.com/acados/tera_renderer/releases/ and ' + \
                'place them in <acados_root>/bin (please strip the ' + \
                'version and platform from the binaries e.g. ' + \
                't_renderer-v0.0.20 -> t_renderer).')

    model = acados_ocp.model
    if acados_ocp.solver_options.integrator_type == 'ERK':
        # explicit model -- generate C code
        generate_c_code_explicit_ode(model)
    else:
        # implicit model -- generate C code
        opts = dict(generate_hess=1)
        generate_c_code_implicit_ode(model, opts)
    
    if acados_ocp.constraints.constr_type == 'BGP' and acados_ocp.dims.nphi > 0:
        # nonlinear part of nonlinear constraints 
        generate_c_code_constraint(acados_ocp.con_phi)
    elif acados_ocp.constraints.constr_type  == 'BGH' and acados_ocp.dims.nh > 0: 
        generate_c_code_constraint(acados_ocp.con_h)

    if acados_ocp.constraints.constr_type_e  == 'BGP' and acados_ocp.dims.nphi_e > 0:
        # nonlinear part of nonlinear constraints 
        generate_c_code_constraint_e(acados_ocp.con_phi_e)
    elif acados_ocp.constraints.constr_type_e  == 'BGH' and acados_ocp.dims.nh_e > 0: 
        generate_c_code_constraint_e(acados_ocp.con_h_e)
    
    if acados_ocp.cost.cost_type == 'NONLINEAR_LS':
        generate_c_code_nls_cost(acados_ocp.cost_r)

    if acados_ocp.cost.cost_type_e == 'NONLINEAR_LS':
        generate_c_code_nls_cost_e(acados_ocp.cost_r_e)

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
        json.dump(ocp_nlp, f, default=np_array_to_list)

    with open(json_file, 'r') as f:
        ocp_nlp_json = json.load(f)

    ocp_nlp_dict = json2dict(ocp_nlp_json, ocp_nlp_json['dims'])

    # setting up loader and environment
    template_glob = acados_path + '/interfaces/acados_template/acados_template/c_templates_tera/*'
    acados_template_path = acados_path + '/interfaces/acados_template/acados_template/c_templates_tera'

    # create c_generated_code folder
    if not os.path.exists('c_generated_code'):
        os.mkdir('c_generated_code')

    os.chdir('c_generated_code')
    # render source template
    template_file = 'main.in.c'
    out_file = 'main_' + model.name + '.c'
    # output file
    os_cmd = tera_path + 't_renderer ' + "\"" + template_glob + "\"" + ' ' + "\"" \
            + template_file + "\"" + ' ' + "\"" + '../' + json_file + \
            "\"" + ' ' + "\"" + out_file + "\""

    status = os.system(os_cmd)
    if (status != 0):
        raise Exception('Rendering of {} failed! Exiting.\n'.format(template_file))

    os.chdir('..')
        
    os.chdir('c_generated_code')
    # render source template
    template_file = 'acados_solver.in.c'
    out_file = 'acados_solver_' + model.name + '.c'
    # output file
    os_cmd = tera_path + 't_renderer ' + "\"" + template_glob + "\"" + ' ' + "\"" \
            + template_file + "\"" + ' ' + "\"" + '../' + json_file + \
            "\"" + ' ' + "\"" + out_file + "\""

    status = os.system(os_cmd)
    if (status != 0):
        raise Exception('Rendering of {} failed! Exiting.\n'.format(template_file))

    os.chdir('..')

    os.chdir('c_generated_code')
    # render source template
    template_file = 'acados_solver.in.h'
    out_file = 'acados_solver_' + model.name + '.h'
    # output file
    os_cmd = tera_path + 't_renderer ' + "\"" + template_glob + "\"" + ' ' + "\"" \
            + template_file + "\"" + ' ' + "\"" + '../' + json_file + \
            "\"" + ' ' + "\"" + out_file + "\""

    status = os.system(os_cmd)
    if (status != 0):
        raise Exception('Rendering of {} failed! Exiting.\n'.format(template_file))

    os.chdir('..')

    os.chdir('c_generated_code')
    # render source template
    template_file = 'acados_sim_solver.in.c'
    out_file = 'acados_sim_solver_' + model.name + '.c'
    # output file
    os_cmd = tera_path + 't_renderer ' + "\"" + template_glob + "\"" + ' ' + "\"" \
            + template_file + "\"" + ' ' + "\"" + '../' + json_file + \
            "\"" + ' ' + "\"" + out_file + "\""

    status = os.system(os_cmd)
    if (status != 0):
        raise Exception('Rendering of {} failed! Exiting.\n'.format(template_file))

    os.chdir('..')

    os.chdir('c_generated_code')
    # render source template
    template_file = 'acados_sim_solver.in.h'
    out_file = 'acados_sim_solver_' + model.name + '.h'
    # output file
    os_cmd = tera_path + 't_renderer ' + "\"" + template_glob + "\"" + ' ' + "\"" \
            + template_file + "\"" + ' ' + "\"" + '../' + json_file + \
            "\"" + ' ' + "\"" + out_file + "\""

    status = os.system(os_cmd)
    if (status != 0):
        raise Exception('Rendering of {} failed! Exiting.\n'.format(template_file))

    os.chdir('..')

    os.chdir('c_generated_code/' + model.name + '_model/')
    # render source template
    template_file = 'model.in.h'
    out_file = model.name + '_model.h'
    # output file
    os_cmd = tera_path + 't_renderer ' + "\"" + template_glob + "\"" + ' ' + "\"" \
            + template_file + "\"" + ' ' + "\"" + '../../' + json_file + \
            "\"" + ' ' + "\"" + out_file + "\""

    status = os.system(os_cmd)
    if (status != 0):
        raise Exception('Rendering of {} failed! Exiting.\n'.format(template_file))

    os.chdir('../..')

    if acados_ocp.constraints.constr_type == 'BGP' and acados_ocp.dims.nphi > 0:
        # create folder
        if not os.path.exists('c_generated_code/' + acados_ocp.con_phi.name + '_phi_constraint/'):
            os.mkdir('c_generated_code/' + acados_ocp.con_phi.name + '_phi_constraint/')
        os.chdir('c_generated_code/' + acados_ocp.con_phi.name + '_phi_constraint/')
        # render source template
        template_file = 'phi_constraint.in.h'
        out_file = acados_ocp.con_phi.name + '_phi_constraint.h'
        # output file
        os_cmd = tera_path + 't_renderer ' + "\"" + template_glob + "\"" + ' ' + "\"" \
                + template_file + "\"" + ' ' + "\"" + '../../' + json_file + \
                "\"" + ' ' + "\"" + out_file + "\""

        status = os.system(os_cmd)
        if (status != 0):
            raise Exception('Rendering of {} failed! Exiting.\n'.format(template_file))

        os.chdir('../..')
        # create folder
        if not os.path.exists('c_generated_code/' + acados_ocp.con_phi.name + '_r_constraint/'):
            os.mkdir('c_generated_code/' + acados_ocp.con_phi.name + '_r_constraint/')
        os.chdir('c_generated_code/' + acados_ocp.con_phi.name + '_r_constraint/')
        # render source template
        template_file = 'r_constraint.in.h'
        out_file = acados_ocp.con_phi.name + '_r_constraint.h'
        # output file
        os_cmd = tera_path + 't_renderer ' + "\"" + template_glob + "\"" + ' ' + "\"" \
                + template_file + "\"" + ' ' + "\"" + '../../' + json_file + \
                "\"" + ' ' + "\"" + out_file + "\""

        status = os.system(os_cmd)
        if (status != 0):
            raise Exception('Rendering of {} failed! Exiting.\n'.format(template_file))

        os.chdir('../..')

    if acados_ocp.constraints.constr_type_e == 'BGP' and acados_ocp.dims.nphi_e > 0:
        # create folder
        if not os.path.exists('c_generated_code/' + acados_ocp.con_phi_e.name + '_phi_e_constraint/'):
            os.mkdir('c_generated_code/' + acados_ocp.con_phi_e.name + '_phi_e_constraint/')
        os.chdir('c_generated_code/' + acados_ocp.con_phi_e.name + '_phi_e_constraint/')
        # render source template
        template_file = 'phi_e_constraint.in.h'
        out_file = acados_ocp.con_phi_e.name + '_phi_e_constraint.h'
        # output file
        os_cmd = tera_path + 't_renderer ' + "\"" + template_glob + "\"" + ' ' + "\"" \
                + template_file + "\"" + ' ' + "\"" + '../../' + json_file + \
                "\"" + ' ' + "\"" + out_file + "\""

        status = os.system(os_cmd)
        if (status != 0):
            raise Exception('Rendering of {} failed! Exiting.\n'.format(template_file))

        os.chdir('../..')
        # create folder
        if not os.path.exists('c_generated_code/' + acados_ocp.con_phi_e.name + '_r_e_constraint/'):
            os.mkdir('c_generated_code/' + acados_ocp.con_phi_e.name + '_r_e_constraint/')
        os.chdir('c_generated_code/' + acados_ocp.con_phi_e.name + '_r_e_constraint/')
        # render source template
        template_file = 'r_e_constraint.in.h'
        out_file = acados_ocp.con_phi_e.name + '_r_e_constraint.h'
        # output file
        os_cmd = tera_path + 't_renderer ' + "\"" + template_glob + "\"" + ' ' + "\"" \
                + template_file + "\"" + ' ' + "\"" + '../../' + json_file + \
                "\"" + ' ' + "\"" + out_file + "\""

        status = os.system(os_cmd)
        if (status != 0):
            raise Exception('Rendering of {} failed! Exiting.\n'.format(template_file))

        os.chdir('../..')

    if acados_ocp.constraints.constr_type == 'BGH' and acados_ocp.dims.nh > 0:
        # create folder
        if not os.path.exists('c_generated_code/' + acados_ocp.con_h.name + '_h_constraint/'):
            os.mkdir('c_generated_code/' + acados_ocp.con_h.name + '_h_constraint/')
        os.chdir('c_generated_code/' + acados_ocp.con_h.name + '_h_constraint/')
        # render source template
        template_file = 'h_constraint.in.h'
        out_file = acados_ocp.con_h.name + '_h_constraint.h'
        # output file
        os_cmd = tera_path + 't_renderer ' + "\"" + template_glob + "\"" + ' ' + "\"" \
                + template_file + "\"" + ' ' + "\"" + '../../' + json_file + \
                "\"" + ' ' + "\"" + out_file + "\""

        status = os.system(os_cmd)
        if (status != 0):
            raise Exception('Rendering of {} failed! Exiting.\n'.format(template_file))

        os.chdir('../..')

    if acados_ocp.constraints.constr_type_e == 'BGH' and acados_ocp.dims.nh_e > 0:
        # create folder
        if not os.path.exists('c_generated_code/' + acados_ocp.con_h_e.name + '_h_e_constraint/'):
            os.mkdir('c_generated_code/' + acados_ocp.con_h_e.name + '_h_e_constraint/')
        os.chdir('c_generated_code/' + acados_ocp.con_h_e.name + '_h_e_constraint/')
        # render source template
        template_file = 'h_e_constraint.in.h'
        out_file = acados_ocp.con_h_e.name + '_h_e_constraint.h'
        # output file
        os_cmd = tera_path + 't_renderer ' + "\"" + template_glob + "\"" + ' ' + "\"" \
                + template_file + "\"" + ' ' + "\"" + '../../' + json_file + \
                "\"" + ' ' + "\"" + out_file + "\""

        status = os.system(os_cmd)
        if (status != 0):
            raise Exception('Rendering of {} failed! Exiting.\n'.format(template_file))

        os.chdir('../..')

    if acados_ocp.cost.cost_type == 'NONLINEAR_LS':
        # create folder
        if not os.path.exists('c_generated_code/' + acados_ocp.cost_r.name + '_r_cost/'):
            os.mkdir('c_generated_code/' + acados_ocp.cost_r.name + '_r_cost/')
        os.chdir('c_generated_code/' + acados_ocp.cost_r.name + '_r_cost/')
        # render source template
        template_file = 'r_cost.in.h'
        out_file = acados_ocp.cost_r.name + '_r_cost.h'
        # output file
        os_cmd = tera_path + 't_renderer ' + "\"" + template_glob + "\"" + ' ' + "\"" \
                + template_file + "\"" + ' ' + "\"" + '../../' + json_file + \
                "\"" + ' ' + "\"" + out_file + "\""

        status = os.system(os_cmd)
        if (status != 0):
            raise Exception('Rendering of {} failed! Exiting.\n'.format(template_file))

        os.chdir('../..')

    if acados_ocp.cost.cost_type_e == 'NONLINEAR_LS':
        # create folder
        if not os.path.exists('c_generated_code/' + acados_ocp.cost_r_e.name + '_r_e_cost/'):
            os.mkdir('c_generated_code/' + acados_ocp.cost_r_e.name + '_r_e_cost/')
        os.chdir('c_generated_code/' + acados_ocp.cost_r_e.name + '_r_e_cost/')
        # render source template
        template_file = 'r_e_cost.in.h'
        out_file = acados_ocp.cost_r_e.name + '_r_e_cost.h'
        # output file
        os_cmd = tera_path + 't_renderer ' + "\"" + template_glob + "\"" + ' ' + "\"" \
                + template_file + "\"" + ' ' + "\"" + '../../' + json_file + \
                "\"" + ' ' + "\"" + out_file + "\""

        status = os.system(os_cmd)
        if (status != 0):
            raise Exception('Rendering of {} failed! Exiting.\n'.format(template_file))

        os.chdir('../..')

    os.chdir('c_generated_code/') 
    # render source template
    template_file = 'Makefile.in'
    out_file = 'Makefile'
    # output file
    os_cmd = tera_path + 't_renderer ' + "\"" + template_glob + "\"" + ' ' + "\"" \
            + template_file + "\"" + ' ' + "\"" + '../' + json_file + \
            "\"" + ' ' + "\"" + out_file + "\""

    status = os.system(os_cmd)
    if (status != 0):
        raise Exception('Rendering of {} failed! Exiting.\n'.format(template_file))

    os.chdir('..')

    os.chdir('c_generated_code/') 
    # render source template
    template_file = 'acados_solver_sfun.in.c'
    out_file = 'acados_solver_sfunction_'  + model.name + '.c'
    # output file
    os_cmd = tera_path + 't_renderer ' + "\"" + template_glob + "\"" + ' ' + "\"" \
            + template_file + "\"" + ' ' + "\"" + '../' + json_file + \
            "\"" + ' ' + "\"" + out_file + "\""

    status = os.system(os_cmd)
    if (status != 0):
        raise Exception('Rendering of {} failed! Exiting.\n'.format(template_file))

    os.chdir('..')

    os.chdir('c_generated_code/') 
    # render source template
    template_file = 'make_sfun.in.m'
    out_file = 'make_sfun.m'
    # output file
    os_cmd = tera_path + 't_renderer ' + "\"" + template_glob + "\"" + ' ' + "\"" \
            + template_file + "\"" + ' ' + "\"" + '../' + json_file + \
            "\"" + ' ' + "\"" + out_file + "\""

    status = os.system(os_cmd)
    if (status != 0):
        raise Exception('Rendering of {} failed! Exiting.\n'.format(template_file))

    os.chdir('..')

    # make 
    os.chdir('c_generated_code')
    os.system('make')
    os.system('make shared_lib')
    os.chdir('..')

    solver = acados_solver(acados_ocp, 'c_generated_code/libacados_solver_' + model.name + '.so')
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
                    \n Possible values are {}. Exiting.".format(out_fields, fields))

        self.shared_lib.ocp_nlp_dims_get_from_attr.argtypes = [c_void_p, c_void_p, c_void_p, c_int, c_char_p]
        self.shared_lib.ocp_nlp_dims_get_from_attr.restype = c_int

        dims = self.shared_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, self.nlp_dims, self.nlp_out, stage_, field)

        out = np.ascontiguousarray(np.zeros((dims,)), dtype=np.float64)
        out_data = cast(out.ctypes.data, POINTER(c_double))

        self.shared_lib.ocp_nlp_out_get.argtypes = [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]
        self.shared_lib.ocp_nlp_out_get(self.nlp_config, self.nlp_dims, self.nlp_out, stage_, field, out_data);

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
            if (field_ not in constraints_fields) and (field_ not in cost_fields) and (field_ not in out_fields):
                raise Exception("acados_solver: {} is not a valid key for method `set(value)`.\
                        \nPossible values are {} and {}. Exiting.".format(field, cost_fields, constraints_fields, out_fields))

            self.shared_lib.ocp_nlp_dims_get_from_attr.argtypes = [c_void_p, c_void_p, c_void_p, c_int, c_char_p]
            self.shared_lib.ocp_nlp_dims_get_from_attr.restype = c_int

            dims = self.shared_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, self.nlp_dims, self.nlp_out, stage_, field)
             
            if value_.shape[0] != dims: 
                raise Exception('acados_solver.set(): mismatching dimension for field "{}" with dimension {} (you have {})'.format(field_,dims, value_.shape[0]))

            value_data = cast(value_.ctypes.data, POINTER(c_double))
            value_data_p = cast((value_data), c_void_p)

            if field_ in constraints_fields:
                self.shared_lib.ocp_nlp_constraints_model_set.argtypes = [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]
                self.shared_lib.ocp_nlp_constraints_model_set(self.nlp_config, self.nlp_dims, self.nlp_in, stage, field, value_data_p);
            elif field_ in cost_fields:
                self.shared_lib.ocp_nlp_cost_model_set.argtypes = [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]
                self.shared_lib.ocp_nlp_cost_model_set(self.nlp_config, self.nlp_dims, self.nlp_in, stage, field, value_data_p);
            elif field_ in out_fields:
                self.shared_lib.ocp_nlp_out_set.argtypes = [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]
                self.shared_lib.ocp_nlp_out_set(self.nlp_config, self.nlp_dims, self.nlp_out, stage, field, value_data_p);

        return

    def __del__(self):
        self.shared_lib.acados_free()


