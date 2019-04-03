from acados_template import *
import acados_template
from jinja2 import Environment, FileSystemLoader
import argparse
import json
import os 

def generate_solver_matlab(acados_ocp_nlp_json_file):

    acados_path = os.path.dirname(acados_template.__file__)
    # load MATLAB JSON file instead
    
    with open(acados_ocp_nlp_json_file[0], 'r') as f:
        ocp_nlp_json = json.load(f)

    model_name = ocp_nlp_json['model_name']

    # load JSON layout
    with open(acados_path + '/../acados_layout.json', 'r') as f:
        ocp_nlp_layout = json.load(f)

    ocp_nlp_dict = json2dict(ocp_nlp_json, ocp_nlp_json['dims'])

    ra = ocp_nlp_as_object(ocp_nlp_dict)
    ra.cost = ocp_nlp_as_object(ra.cost)
    ra.constraints = ocp_nlp_as_object(ra.constraints)
    ra.solver_config = ocp_nlp_as_object(ra.solver_config)
    ra.dims = ocp_nlp_as_object(ra.dims)

    # setting up loader and environment
    file_loader = FileSystemLoader(acados_path + '/c_templates')
    env = Environment(loader = file_loader)

    # render source template
    template = env.get_template('main.in.c')
    output = template.render(ra=ra)
    # output file
    out_file = open('./c_generated_code/main_' + model_name + '.c', 'w+')
    out_file.write(output)

    # render source template
    template = env.get_template('acados_solver.in.c')
    output = template.render(ra=ra)
    # output file
    out_file = open('./c_generated_code/acados_solver_' + model_name + '.c', 'w+')
    out_file.write(output)

    # render source template
    template = env.get_template('acados_solver.in.h')
    output = template.render(ra=ra)
    # output file
    out_file = open('./c_generated_code/acados_solver_' + model_name + '.h', 'w+')
    out_file.write(output)

    # render header templates
    template = env.get_template('model.in.h')
    output = template.render(ra=ra)
    # output file
    out_file = open('./c_generated_code/' + model_name + '_model/' + model_name + '_model.h', 'w+')
    out_file.write(output)

    if ra.dims.npd > 0:
        # render header templates
        template = env.get_template('p_constraint.in.h')
        output = template.render(ra=ra)
        # output file
        out_file = open('./c_generated_code/' + ra.con_p_name + '_p_constraint/' + ra.con_p_name + '_p_constraint.h', 'w+')
        out_file.write(output)

    if ra.dims.nh > 0:
        # render header templates
        template = env.get_template('h_constraint.in.h')
        output = template.render(ra=ra)
        # output file
        out_file = open('./c_generated_code/' + ra.con_h_name + '_h_constraint/' + ra.con_h_name + '_h_constraint.h', 'w+')
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
    out_file = open('./c_generated_code/acados_solver_sfunction_'  + model_name + '.c', 'w+')
    out_file.write(output)

    # render MATLAB make script
    template = env.get_template('make_sfun.in.m')
    output = template.render(ra=ra)

    # output file
    out_file = open('./c_generated_code/make_sfun.m', 'w+')
    out_file.write(output)
    
    print('Successfully generated acados solver!\n')

    # build generated code
    os.chdir('c_generated_code')
    os.system('make')
    os.system('make shared_lib')
    os.chdir('..')

    print('Successfully built generated code!')


parser = argparse.ArgumentParser(description='Generate acados solver.')
parser.add_argument('--json_file_name', dest='json_file', nargs='+',
                    help='name of the JSON file containing the acados ocp_nlp description')

args = parser.parse_args()
generate_solver_matlab(args.json_file)

