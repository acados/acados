from jinja2 import Environment, FileSystemLoader
from .generate_c_code_explicit_ode import *
from .generate_c_code_implicit_ode import *
from .generate_c_code_constraint import *
from .acados_ocp_nlp import *
from ctypes import *

def generate_solver(acados_ocp, json_file='acados_ocp_nlp.json'):
    USE_TERA = 1 # EXPERIMENTAL: use Tera standalone parser instead of Jinja2
    
    model = acados_ocp.model
    if acados_ocp.solver_config.integrator_type == 'ERK':
        # explicit model -- generate C code
        generate_c_code_explicit_ode(model)
    else:
        # implicit model -- generate C code
        opts = dict(generate_hess=1)
        generate_c_code_implicit_ode(model, opts)
    
    if acados_ocp.con_p.name is not None and acados_ocp.con_h.name is None:
        raise Exception('h constraint is missing!')

    if acados_ocp.con_h.name is not None:
        # nonlinear part of nonlinear constraints 
        generate_c_code_constraint(acados_ocp.con_h, '_h_constraint')

    if acados_ocp.con_h_e.name is not None:
        # nonlinear part of nonlinear constraints 
        generate_c_code_constraint(acados_ocp.con_h_e, '_h_e_constraint')
    
    if acados_ocp.con_p.name is not None:
        # convex part of nonlinear constraints 
        generate_c_code_constraint(acados_ocp.con_p, '_p_constraint')

    ocp_nlp = acados_ocp
    ocp_nlp.cost = acados_ocp.cost.__dict__
    ocp_nlp.constraints = acados_ocp.constraints.__dict__
    ocp_nlp.solver_config = acados_ocp.solver_config.__dict__
    ocp_nlp.dims = acados_ocp.dims.__dict__
    ocp_nlp.con_h = acados_ocp.con_h.__dict__
    ocp_nlp.con_h_e = acados_ocp.con_h_e.__dict__
    ocp_nlp.con_p = acados_ocp.con_p.__dict__
    ocp_nlp.con_p_e = acados_ocp.con_p_e.__dict__
    ocp_nlp.model = acados_ocp.model.__dict__
    ocp_nlp = ocp_nlp.__dict__

    # need to strip non-numerical stuff from expressions for now
    ocp_nlp['con_h'] = acados_constraint_strip_non_num(ocp_nlp['con_h'])
    ocp_nlp['con_p'] = acados_constraint_strip_non_num(ocp_nlp['con_p'])
    ocp_nlp['con_h_e'] = acados_constraint_strip_non_num(ocp_nlp['con_h_e'])
    ocp_nlp['con_p_e'] = acados_constraint_strip_non_num(ocp_nlp['con_p_e'])

    ocp_nlp['model'] = acados_dae_strip_non_num(ocp_nlp['model'])

    ocp_nlp = dict2json(ocp_nlp)

    with open(json_file, 'w') as f:
        json.dump(ocp_nlp, f, default=np_array_to_list)

    with open(json_file, 'r') as f:
        ocp_nlp_json = json.load(f)

    ocp_nlp_dict = json2dict(ocp_nlp_json, ocp_nlp_json['dims'])
    acados_ocp = ocp_nlp_as_object(ocp_nlp_dict)
    acados_ocp.cost = ocp_nlp_as_object(acados_ocp.cost)
    acados_ocp.constraints = ocp_nlp_as_object(acados_ocp.constraints)
    acados_ocp.solver_config = ocp_nlp_as_object(acados_ocp.solver_config)
    acados_ocp.dims = ocp_nlp_as_object(acados_ocp.dims)

    acados_ocp.con_h = ocp_nlp_as_object(acados_ocp.con_h)
    acados_ocp.con_h_e = ocp_nlp_as_object(acados_ocp.con_h_e)
    acados_ocp.con_p = ocp_nlp_as_object(acados_ocp.con_p)
    acados_ocp.con_p_e = ocp_nlp_as_object(acados_ocp.con_p_e)

    # setting up loader and environment
    acados_path = os.path.dirname(os.path.abspath(__file__))
    if USE_TERA == 0:
        file_loader = FileSystemLoader(acados_path + '/c_templates')
        env = Environment(loader = file_loader)
    else:
        template_glob = acados_path + '/c_templates_tera/*'
        acados_template_path = acados_path + '/c_templates_tera'


    # check render arguments
    check_ra(acados_ocp)

    # create c_generated_code folder
    if not os.path.exists('c_generated_code'):
        os.mkdir('c_generated_code')

    if USE_TERA == 0:
        # render source template
        template = env.get_template('main.in.c')
        output = template.render(ocp=acados_ocp)
        # output file
        out_file = open('./c_generated_code/main_' + model.name + '.c', 'w+')
        out_file.write(output)
    else:
        os.chdir('c_generated_code')
        # render source template
        template_file = 'main.in.c'
        out_file = 'main_' + model.name + '.c'
        # output file
        os_cmd = 't_renderer ' + "\"" + template_glob + "\"" + ' ' + "\"" \
                + template_file + "\"" + ' ' + "\"" + '../' + json_file + \
                "\"" + ' ' + "\"" + out_file + "\""

        os.system(os_cmd)
        os.chdir('..')
        
    if USE_TERA == 0:
        # render source template
        template = env.get_template('acados_solver.in.c')
        output = template.render(ocp=acados_ocp)
        # output file
        out_file = open('./c_generated_code/acados_solver_' + model.name + '.c', 'w+')
        out_file.write(output)
    else:
        os.chdir('c_generated_code')
        # render source template
        template_file = 'acados_solver.in.c'
        out_file = 'acados_solver_' + model.name + '.c'
        # output file
        os_cmd = 't_renderer ' + "\"" + template_glob + "\"" + ' ' + "\"" \
                + template_file + "\"" + ' ' + "\"" + '../' + json_file + \
                "\"" + ' ' + "\"" + out_file + "\""

        os.system(os_cmd)
        os.chdir('..')

    if USE_TERA == 0:
        # render source template
        template = env.get_template('acados_solver.in.h')
        output = template.render(ocp=acados_ocp)
        # output file
        out_file = open('./c_generated_code/acados_solver_' + model.name + '.h', 'w+')
        out_file.write(output)
    else:
        os.chdir('c_generated_code')
        # render source template
        template_file = 'acados_solver.in.h'
        out_file = 'acados_solver_' + model.name + '.h'
        # output file
        os_cmd = 't_renderer ' + "\"" + template_glob + "\"" + ' ' + "\"" \
                + template_file + "\"" + ' ' + "\"" + '../' + json_file + \
                "\"" + ' ' + "\"" + out_file + "\""

        os.system(os_cmd)
        os.chdir('..')

    if USE_TERA == 0:
        # render header templates
        template = env.get_template('model.in.h')
        output = template.render(ocp=acados_ocp)
        # output file
        out_file = open('./c_generated_code/' + model.name + '_model/' + model.name + '_model.h', 'w+')
        out_file.write(output)
    else:
        os.chdir('c_generated_code/' + model.name + '_model/')
        # render source template
        template_file = 'model.in.h'
        out_file = model.name + '_model.h'
        # output file
        os_cmd = 't_renderer ' + "\"" + template_glob + "\"" + ' ' + "\"" \
                + template_file + "\"" + ' ' + "\"" + '../../' + json_file + \
                "\"" + ' ' + "\"" + out_file + "\""

        os.system(os_cmd)
        os.chdir('../..')

    if acados_ocp.dims.npd > 0:
        # create folder
        if not os.path.exists('c_generated_code/' + acados_ocp.con_p.name + '_p_constraint/'):
            os.mkdir('c_generated_code/' + acados_ocp.con_p.name + '_p_constraint/')
        if USE_TERA == 0:
            # render header templates
            template = env.get_template('p_constraint.in.h')
            output = template.render(ocp=acados_ocp)
            # output file
            out_file = open('./c_generated_code/' + acados_ocp.con_p.name + '_p_constraint/' + acados_ocp.con_p.name + '_p_constraint.h', 'w+')
            out_file.write(output)
        else:
            os.chdir('c_generated_code/' + acados_ocp.con_p.name + '_p_constraint/')
            # render source template
            template_file = 'p_constraint.in.h'
            out_file = acados_ocp.con_p.name + '_p_constraint.h'
            # output file
            os_cmd = 't_renderer ' + "\"" + template_glob + "\"" + ' ' + "\"" \
                    + template_file + "\"" + ' ' + "\"" + '../../' + json_file + \
                    "\"" + ' ' + "\"" + out_file + "\""

            os.system(os_cmd)
            os.chdir('../..')

    if acados_ocp.dims.nh > 0:
        # create folder
        if not os.path.exists('c_generated_code/' + acados_ocp.con_h.name + '_h_constraint/'):
            os.mkdir('c_generated_code/' + acados_ocp.con_h.name + '_h_constraint/')
        if USE_TERA == 0:
            # render header templates
            template = env.get_template('h_constraint.in.h')
            output = template.render(ocp=acados_ocp)
            # output file
            out_file = open('./c_generated_code/' + acados_ocp.con_h.name + '_h_constraint/' + acados_ocp.con_h.name + '_h_constraint.h', 'w+')
            out_file.write(output)
        else:
            os.chdir('c_generated_code/' + acados_ocp.con_h.name + '_h_constraint/')
            # render source template
            template_file = 'h_constraint.in.h'
            out_file = acados_ocp.con_h.name + '_h_constraint.h'
            # output file
            os_cmd = 't_renderer ' + "\"" + template_glob + "\"" + ' ' + "\"" \
                    + template_file + "\"" + ' ' + "\"" + '../../' + json_file + \
                    "\"" + ' ' + "\"" + out_file + "\""

            os.system(os_cmd)
            os.chdir('../..')

    if acados_ocp.dims.nh_e > 0:
        # create folder
        if not os.path.exists('c_generated_code/' + acados_ocp.con_h_e.name + '_h_e_constraint/'):
            os.mkdir('c_generated_code/' + acados_ocp.con_h_e.name + '_h_e_constraint/')
        if USE_TERA == 0:
            # render header templates
            template = env.get_template('h_e_constraint.in.h')
            output = template.render(ocp=acados_ocp)
            # output file
            out_file = open('./c_generated_code/' + acados_ocp.con_h_e.name + '_h_e_constraint/' + acados_ocp.con_h_e.name + '_h_e_constraint.h', 'w+')
            out_file.write(output)
        else:
            os.chdir('c_generated_code/' + acados_ocp.con_h_e.name + '_h_e_constraint/')
            # render source template
            template_file = 'h_e_constraint.in.h'
            out_file = acados_ocp.con_h_e.name + '_h_e_constraint.h'
            # output file
            os_cmd = 't_renderer ' + "\"" + template_glob + "\"" + ' ' + "\"" \
                    + template_file + "\"" + ' ' + "\"" + '../../' + json_file + \
                    "\"" + ' ' + "\"" + out_file + "\""

            os.system(os_cmd)
            os.chdir('../..')

    if USE_TERA == 0:
        # render Makefile template
        template = env.get_template('Makefile.in')
        output = template.render(ocp=acados_ocp)

        # output file
        out_file = open('./c_generated_code/Makefile', 'w+')
        out_file.write(output)
    else:
        os.chdir('c_generated_code/') 
        # render source template
        template_file = 'Makefile.in'
        out_file = 'Makefile'
        # output file
        os_cmd = 't_renderer ' + "\"" + template_glob + "\"" + ' ' + "\"" \
                + template_file + "\"" + ' ' + "\"" + '../' + json_file + \
                "\"" + ' ' + "\"" + out_file + "\""

        os.system(os_cmd)
        os.chdir('..')

    if USE_TERA == 0:
        # render S-Function template
        template = env.get_template('acados_solver_sfun.in.c')
        output = template.render(ocp=acados_ocp)

        # output file
        out_file = open('./c_generated_code/acados_solver_sfunction_'  + model.name + '.c', 'w+')
        out_file.write(output)
    else:
        os.chdir('c_generated_code/') 
        # render source template
        template_file = 'acados_solver_sfun.in.c'
        out_file = 'acados_solver_sfunction_'  + model.name + '.c'
        # output file
        os_cmd = 't_renderer ' + "\"" + template_glob + "\"" + ' ' + "\"" \
                + template_file + "\"" + ' ' + "\"" + '../' + json_file + \
                "\"" + ' ' + "\"" + out_file + "\""

        os.system(os_cmd)
        os.chdir('..')

    if USE_TERA == 0:
        # render MATLAB make script
        template = env.get_template('make_sfun.in.m')
        output = template.render(ocp=acados_ocp)

        # output file
        out_file = open('./c_generated_code/make_sfun.m', 'w+')
        out_file.write(output)
    else:
        os.chdir('c_generated_code/') 
        # render source template
        template_file = 'make_sfun.in.m'
        out_file = 'acados_solver_sfun.in.c'
        # output file
        os_cmd = 't_renderer ' + "\"" + template_glob + "\"" + ' ' + "\"" \
                + template_file + "\"" + ' ' + "\"" + '../' + json_file + \
                "\"" + ' ' + "\"" + out_file + "\""

        os.system(os_cmd)
        os.chdir('..')

    # make 
    os.chdir('c_generated_code')
    os.system('make')
    os.system('make shared_lib')
    os.chdir('..')

    solver = acados_solver(acados_ocp, 'c_generated_code/acados_solver_' + model.name + '.so')
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

        self.acados_ocp = acados_ocp

    def solve(self):
        status = self.shared_lib.acados_solve()
        return status

    def get(self, stage_, field_):

        field = field_
        field = field.encode('utf-8')

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
        
        cost = ['y_ref', 'yref']
        constraints = ['lbx', 'ubx', 'lbu', 'ubu']

        # cast value_ to avoid conversion issues
        value_ = value_.astype(float)

        field = field_
        field = field.encode('utf-8')

        self.shared_lib.ocp_nlp_dims_get_from_attr.argtypes = [c_void_p, c_void_p, c_void_p, c_int, c_char_p]
        self.shared_lib.ocp_nlp_dims_get_from_attr.restype = c_int

        dims = self.shared_lib.ocp_nlp_dims_get_from_attr(self.nlp_config, self.nlp_dims, self.nlp_out, stage_, field)
         
        if value_.shape[0] != dims: 
            raise Exception('acados_solver.set(): mismatching dimension for field "{}" with dimension {} (you have {})'.format(field_,dims, value_.shape[0]))

        value_data = cast(value_.ctypes.data, POINTER(c_double))
        value_data_p = cast((value_data), c_void_p)

        stage = c_int(stage_)
        if field_ in constraints:
            self.shared_lib.ocp_nlp_constraints_model_set.argtypes = [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]
            self.shared_lib.ocp_nlp_constraints_model_set(self.nlp_config, self.nlp_dims, self.nlp_in, stage, field, value_data_p);
        if field_ in cost:
            self.shared_lib.ocp_nlp_cost_model_set.argtypes = [c_void_p, c_void_p, c_void_p, c_int, c_char_p, c_void_p]
            self.shared_lib.ocp_nlp_cost_model_set(self.nlp_config, self.nlp_dims, self.nlp_in, stage, field, value_data_p);

        return

    def __del__(self):
        self.shared_lib.acados_free()


