#
# Copyright (c) The acados authors.
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

import itertools as it
import os, json
from random import sample

TEST_SAMPLE = True # only test random sample instead of all possible combinations
SAMPLE_SIZE = 30
GENERATE_SAMPLE = False # generate a sample of combinations

COST_MODULE_values = ['EXTERNAL', 'LS', 'NLS']
COST_MODULE_N_values = ['EXTERNAL', 'LS', 'NLS']
QP_SOLVER_values = ['PARTIAL_CONDENSING_HPIPM', 'FULL_CONDENSING_HPIPM', 'FULL_CONDENSING_QPOASES']
INTEGRATOR_TYPE_values = ['ERK', 'IRK', 'GNSF']
SOLVER_TYPE_values = ['SQP', 'SQP_RTI']
HESS_APPROX_values = ['GAUSS_NEWTON', 'EXACT']
REGULARIZATION_values = ['NO_REGULARIZE', 'MIRROR', 'PROJECT', 'PROJECT_REDUC_HESS', 'CONVEXIFY']


test_parameters = { 'COST_MODULE_values': COST_MODULE_values,
                    'COST_MODULE_N_values': COST_MODULE_N_values,
                    'HESS_APPROX_values': HESS_APPROX_values,
                    'INTEGRATOR_TYPE_values': INTEGRATOR_TYPE_values,
                    'QP_SOLVER_values': QP_SOLVER_values,
                    'REGULARIZATION_values': REGULARIZATION_values,
                    'SOLVER_TYPE_values': SOLVER_TYPE_values}

all_parameter_names = sorted(test_parameters)

# TEST GAUSS_NEWTON
test_parameters_gn = test_parameters
test_parameters_gn['HESS_APPROX_values'] = ['GAUSS_NEWTON']
test_parameters_gn['REGULARIZATION_values'] = ['NO_REGULARIZE']

combinations = list(it.product(*(test_parameters_gn[Name] for Name in all_parameter_names)))

json_file='test_data/test_combinations_pendulum_GAUSS_NEWTON.json'
if GENERATE_SAMPLE:
    combinations = sample(combinations, SAMPLE_SIZE)
    with open(json_file, 'w') as f:
        json.dump(combinations, f, indent=4, sort_keys=True)

if TEST_SAMPLE:
    with open(json_file, 'r') as f:
        combinations = json.load(f)


for parameters in combinations:
    os_cmd = ("python test_ocp_setting.py" +
        " --COST_MODULE_N {}".format(parameters[0]) +
        " --COST_MODULE {}".format(parameters[1]) +
        " --HESS_APPROX {}".format(parameters[2]) +
        " --INTEGRATOR_TYPE {}".format(parameters[3]) +
        " --QP_SOLVER {}".format(parameters[4]) +
        " --REGULARIZATION {}".format(parameters[5]) +
        " --SOLVER_TYPE {}".format(parameters[6])
        )
    status = os.system(os_cmd)
    if status != 0:
        raise Exception("acados status = {} on test {}. Exiting\n".format(status, parameters))


# TEST EXACT HESSIAN
test_parameters_exact = test_parameters
test_parameters_exact['HESS_APPROX_values'] = ['EXACT']
test_parameters_exact['REGULARIZATION_values'] = ['MIRROR', 'PROJECT', 'CONVEXIFY'] #, 'CONVEXIFY', 'PROJECT_REDUC_HESS']
test_parameters_exact['INTEGRATOR_TYPE_values'] = ['ERK', 'IRK']
# test_parameters_exact['COST_MODULE_N_values'] = ['LS', 'NLS', 'EXTERNAL']
# test_parameters_exact['COST_MODULE_values'] = ['LS', 'NLS'] # EXTERNAL

combinations = list(it.product(*(test_parameters_exact[Name] for Name in all_parameter_names)))

json_file='test_data/test_combinations_pendulum_EXACT.json'
if GENERATE_SAMPLE:
    combinations = sample(combinations, SAMPLE_SIZE)
    with open(json_file, 'w') as f:
        json.dump(combinations, f, indent=4, sort_keys=True)

if TEST_SAMPLE:
    with open(json_file, 'r') as f:
        combinations = json.load(f)

# combinations = [('LS', 'LS', 'EXACT', 'ERK', 'PARTIAL_CONDENSING_HPIPM', 'MIRROR', 'SQP')]

for parameters in combinations:
    os_cmd = ("python test_ocp_setting.py" +
        " --COST_MODULE_N {}".format(parameters[0]) +
        " --COST_MODULE {}".format(parameters[1]) +
        " --HESS_APPROX {}".format(parameters[2]) +
        " --INTEGRATOR_TYPE {}".format(parameters[3]) +
        " --QP_SOLVER {}".format(parameters[4]) +
        " --REGULARIZATION {}".format(parameters[5]) +
        " --SOLVER_TYPE {}".format(parameters[6])
        )
    status = os.system(os_cmd)
    if status != 0:
        raise Exception("acados status = {} on test {}. Exiting\n".format(status, parameters))

