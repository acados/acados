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

import itertools as it
import os

COST_MODULE_values = ['EXTERNAL', 'LS', 'NLS']
COST_MODULE_N_values = ['EXTERNAL', 'LS', 'NLS']
QP_SOLVER_values = ['PARTIAL_CONDENSING_HPIPM', 'FULL_CONDENSING_HPIPM', 'FULL_CONDENSING_QPOASES']
INTEGRATOR_TYPE_values = ['ERK', 'IRK']
SOLVER_TYPE_values = ['SQP', 'SQP_RTI']

test_parameters = { 'COST_MODULE_values': COST_MODULE_values,
                    'COST_MODULE_N_values': COST_MODULE_N_values,
                    'QP_SOLVER_values': QP_SOLVER_values,
                    'INTEGRATOR_TYPE_values': INTEGRATOR_TYPE_values,
                    'SOLVER_TYPE_values': SOLVER_TYPE_values}

all_test_parameters = sorted(test_parameters)
combinations = list(it.product(*(test_parameters[Name] for Name in all_test_parameters)))

for parameters in combinations:
    os_cmd = ("python test_ocp_setting.py" +
        " --COST_MODULE {}".format(parameters[0]) +
        " --COST_MODULE_N {}".format(parameters[1]) +
        " --INTEGRATOR_TYPE {}".format(parameters[2]) +
        " --QP_SOLVER {}".format(parameters[3]) +
        " --SOLVER_TYPE {}".format(parameters[4]))
    status = os.system(os_cmd)
    if status != 0:
        raise Exception("acados status  = {} on test {}. Exiting\n".format(status, parameters))

