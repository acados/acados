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

import os
from casadi import *
from .utils import ALLOWED_CASADI_VERSIONS

def generate_c_code_nls_cost( cost ):

    suffix_name = '_r_cost'
    casadi_version = CasadiMeta.version()
    casadi_opts = dict(mex=False, casadi_int='int', casadi_real='double')

    if  casadi_version not in (ALLOWED_CASADI_VERSIONS):
        msg =  'Please download and install CasADi {} '.format(" or ".join(ALLOWED_CASADI_VERSIONS))
        msg += 'to ensure compatibility with acados.\n'
        msg += 'Version {} currently in use.'.format(casadi_version)
        raise Exception(msg)

    # load cost variables and expression
    x = cost.x
    u = cost.u
    p = cost.p
    cost_exp = cost.expr
    cost_name = cost.name

    # get dimensions
    nx = x.size()[0]
    nu = u.size()[0]

    # set up functions to be exported
    fun_name = cost_name + suffix_name

    cost_jac_exp = transpose(jacobian(cost_exp, vertcat(u, x)))

    nls_cost_fun = Function( fun_name, [x, u, p], \
            [ cost_exp, cost_jac_exp ])

    # generate C code
    if not os.path.exists('c_generated_code'):
        os.mkdir('c_generated_code')

    os.chdir('c_generated_code')
    gen_dir = cost_name + suffix_name
    if not os.path.exists(gen_dir):
        os.mkdir(gen_dir)
    gen_dir_location = './' + gen_dir
    os.chdir(gen_dir_location)
    file_name = cost_name + suffix_name

    nls_cost_fun.generate( file_name, casadi_opts )

    os.chdir('../..')
    return

