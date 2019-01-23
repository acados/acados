from jinja2 import Environment, FileSystemLoader
from .generate_c_code_explicit_ode import *
from .generate_c_code_implicit_ode import *
from .ocp_nlp_render_arguments import *

def generate_solver(model, ra):
    # setting up loader and environment
    acados_path = os.path.dirname(os.path.abspath(__file__))
    file_loader = FileSystemLoader(acados_path + '/c_templates')
    env = Environment(loader = file_loader)

    # explicit model -- generate C code
    generate_c_code_explicit_ode(model);

    # implicit model -- generate C code
    opts = dict(generate_hess=1)
    generate_c_code_implicit_ode(model, opts);

    # check render arguments
    check_ra(ra)

    # render source template
    template = env.get_template('main.in.c')
    output = template.render(ra=ra)
    # output file
    out_file = open('./c_generated_code/main_' + model.name + '.c', 'w+')
    out_file.write(output)

    # render source template
    template = env.get_template('acados_solver.in.c')
    output = template.render(ra=ra)
    # output file
    out_file = open('./c_generated_code/acados_solver_' + model.name + '.c', 'w+')
    out_file.write(output)

    # render source template
    template = env.get_template('acados_solver.in.h')
    output = template.render(ra=ra)
    # output file
    out_file = open('./c_generated_code/acados_solver_' + model.name + '.h', 'w+')
    out_file.write(output)

    # render header template
    template = env.get_template('model.in.h')
    output = template.render(ra=ra)
    # output file
    out_file = open('./c_generated_code/' + model.name + '_model/' + model.name + '_model.h', 'w+')
    out_file.write(output)

    # render Makefile template
    template = env.get_template('Makefile.in')
    output = template.render(ra=ra)

    # output file
    out_file = open('./c_generated_code/Makefile', 'w+')
    out_file.write(output)

    # render S-Function template
    template = env.get_template('acados_solver_sfun.in.c')
    output = template.render(ra=ra)

    # output file
    out_file = open('./c_generated_code/acados_solver_sfunction_'  + model.name + '.c', 'w+')
    out_file.write(output)

    # render MATLAB make script
    template = env.get_template('make_sfun.in.m')
    output = template.render(ra=ra)

    # output file
    out_file = open('./c_generated_code/make_sfun.m', 'w+')
    out_file.write(output)
