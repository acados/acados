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

from acados_template import *
import acados_template
from jinja2 import Environment, FileSystemLoader
import argparse
import json
import os 

def generate_solver_matlab(acados_ocp_nlp_json_file):
    USE_TERA = 1 # EXPERIMENTAL: use Tera standalone parser instead of Jinja2

    acados_path = os.path.dirname(acados_template.__file__)
    # load MATLAB JSON file instead
    
    with open(acados_ocp_nlp_json_file[0], 'r') as f:
        ocp_nlp_json = json.load(f)

    model_name = ocp_nlp_json['model_name']

    # load JSON layout
    with open(acados_path + '/acados_layout.json', 'r') as f:
        ocp_nlp_layout = json.load(f)

    ocp_nlp_dict = json2dict(ocp_nlp_json, ocp_nlp_json['dims'])

    ra = ocp_nlp_as_object(ocp_nlp_dict)
    ra.cost = ocp_nlp_as_object(ra.cost)
    ra.constraints = ocp_nlp_as_object(ra.constraints)
    ra.solver_config = ocp_nlp_as_object(ra.solver_config)
    ra.dims = ocp_nlp_as_object(ra.dims)
    acados_ocp = ra

    # setting up loader and environment
    if USE_TERA == 0:
        file_loader = FileSystemLoader(acados_path + '/c_templates')
        env = Environment(loader = file_loader)
    else:
        template_glob = acados_path + '/c_templates_tera/*'
        acados_template_path = acados_path + '/c_templates_tera'

    # render source template
    if USE_TERA == 0:
        # render source template
        template = env.get_template('main.in.c')
        output = template.render(ocp=acados_ocp)
        # output file
        out_file = open('./c_generated_code/main_' + model_name + '.c', 'w+')
        out_file.write(output)
    else:
        os.chdir('c_generated_code')
        # render source template
        template_file = 'main.in.c'
        out_file = 'main_' + model_name + '.c'
        # output file
        os_cmd = 't_renderer ' + "\"" + template_glob + "\"" + ' ' + "\"" \
                + template_file + "\"" + ' ' + "\"" + '../' + json_file + \
                "\"" + ' ' + "\"" + out_file + "\""

        os.system(os_cmd)
        os.chdir('..')


    # render source template
    if USE_TERA == 0:
        # render source template
        template = env.get_template('acados_solver.in.c')
        output = template.render(ocp=acados_ocp)
        # output file
        out_file = open('./c_generated_code/acados_solver_' + model_name + '.c', 'w+')
        out_file.write(output)
    else:
        os.chdir('c_generated_code')
        # render source template
        template_file = 'acados_solver.in.c'
        out_file = 'acados_solver_' + model_name + '.c'
        # output file
        os_cmd = 't_renderer ' + "\"" + template_glob + "\"" + ' ' + "\"" \
                + template_file + "\"" + ' ' + "\"" + '../' + json_file + \
                "\"" + ' ' + "\"" + out_file + "\""

        os.system(os_cmd)
        os.chdir('..')

    # render source template
    if USE_TERA == 0:
        # render source template
        template = env.get_template('acados_solver.in.h')
        output = template.render(ocp=acados_ocp)
        # output file
        out_file = open('./c_generated_code/acados_solver_' + model_name + '.h', 'w+')
        out_file.write(output)
    else:
        os.chdir('c_generated_code')
        # render source template
        template_file = 'acados_solver.in.h'
        out_file = 'acados_solver_' + model_name + '.h'
        # output file
        os_cmd = 't_renderer ' + "\"" + template_glob + "\"" + ' ' + "\"" \
                + template_file + "\"" + ' ' + "\"" + '../' + json_file + \
                "\"" + ' ' + "\"" + out_file + "\""

        os.system(os_cmd)
        os.chdir('..')

    # render header templates
    if USE_TERA == 0:
        # render header templates
        template = env.get_template('model.in.h')
        output = template.render(ocp=acados_ocp)
        # output file
        out_file = open('./c_generated_code/' + model.name + '_model/' + model_name + '_model.h', 'w+')
        out_file.write(output)
    else:
        os.chdir('c_generated_code/' + model_name + '_model/')
        # render source template
        template_file = 'model.in.h'
        out_file = model.name + '_model.h'
        # output file
        os_cmd = 't_renderer ' + "\"" + template_glob + "\"" + ' ' + "\"" \
                + template_file + "\"" + ' ' + "\"" + '../../' + json_file + \
                "\"" + ' ' + "\"" + out_file + "\""

        os.system(os_cmd)
        os.chdir('../..')

    if ra.dims.npd > 0:
        # render header templates
        if USE_TERA == 0:
            # render header templates
            template = env.get_template('p_constraint.in.h')
            output = template.render(ocp=acados_ocp)
            # output file
            out_file = open('./c_generated_code/' + acados_ocp.con_p_name + '_p_constraint/' + acados_ocp.con_p_name + '_p_constraint.h', 'w+')
            out_file.write(output)
        else:
            os.chdir('c_generated_code/' + acados_ocp.con_p_name + '_p_constraint/')
            # render source template
            template_file = 'p_constraint.in.h'
            out_file = acados_ocp.con_p_name + '_p_constraint.h'
            # output file
            os_cmd = 't_renderer ' + "\"" + template_glob + "\"" + ' ' + "\"" \
                    + template_file + "\"" + ' ' + "\"" + '../../' + json_file + \
                    "\"" + ' ' + "\"" + out_file + "\""

            os.system(os_cmd)
            os.chdir('../..')

    if ra.dims.nh > 0:
        if USE_TERA == 0:
            # render header templates
            template = env.get_template('h_constraint.in.h')
            output = template.render(ocp=acados_ocp)
            # output file
            out_file = open('./c_generated_code/' + acados_ocp.con_h_name + '_h_constraint/' + acados_ocp.con_h_name + '_h_constraint.h', 'w+')
            out_file.write(output)
        else:
            os.chdir('c_generated_code/' + acados_ocp.con_h_name + '_h_constraint/')
            # render source template
            template_file = 'h_constraint.in.h'
            out_file = acados_ocp.con_h_name + '_h_constraint.h'
            # output file
            os_cmd = 't_renderer ' + "\"" + template_glob + "\"" + ' ' + "\"" \
                    + template_file + "\"" + ' ' + "\"" + '../../' + json_file + \
                    "\"" + ' ' + "\"" + out_file + "\""

            os.system(os_cmd)
            os.chdir('../..')

    # render Makefile template
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
        out_file = open('./c_generated_code/acados_solver_sfunction_'  + model_name + '.c', 'w+')
        out_file.write(output)
    else:
        os.chdir('c_generated_code/') 
        # render source template
        template_file = 'acados_solver_sfun.in.c'
        out_file = 'acados_solver_sfunction_'  + model_name + '.c'
        # output file
        os_cmd = 't_renderer ' + "\"" + template_glob + "\"" + ' ' + "\"" \
                + template_file + "\"" + ' ' + "\"" + '../' + json_file + \
                "\"" + ' ' + "\"" + out_file + "\""

        os.system(os_cmd)

    # render MATLAB make script
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

