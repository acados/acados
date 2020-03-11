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
from .utils import ALLOWED_CASADI_VERSIONS, casadi_length

def generate_c_code_nls_cost( model, cost_name, is_terminal ):

    casadi_version = CasadiMeta.version()
    casadi_opts = dict(mex=False, casadi_int='int', casadi_real='double')

    if casadi_version not in (ALLOWED_CASADI_VERSIONS):
        msg =  'Please download and install CasADi {} '.format(" or ".join(ALLOWED_CASADI_VERSIONS))
        msg += 'to ensure compatibility with acados.\n'
        msg += 'Version {} currently in use.'.format(casadi_version)
        raise Exception(msg)

    if is_terminal:
        middle_name = '_cost_y_e'
        u = SX.sym('u', 0, 0)
        cost_expr = model.cost_y_expr_e

    else:
        middle_name = '_cost_y'
        u = model.u
        cost_expr = model.cost_y_expr

    x = model.x
    p = model.p

    # set up directory
    if not os.path.exists('c_generated_code'):
        os.mkdir('c_generated_code')

    os.chdir('c_generated_code')
    gen_dir = cost_name + '_cost'
    if not os.path.exists(gen_dir):
        os.mkdir(gen_dir)
    gen_dir_location = './' + gen_dir
    os.chdir(gen_dir_location)

    # set up expressions
    cost_jac_expr = transpose(jacobian(cost_expr, vertcat(u, x)))

    jac_x = jacobian(cost_expr, x)
    jac_u = jacobian(cost_expr, u)

    ny = casadi_length(cost_expr)

    if isinstance(cost_expr, casadi.SX):
        y = SX.sym('y', ny, 1)
    else:
        y = MX.sym('y', ny, 1)

    y_adj = jtimes(cost_expr, vertcat(u, x), y, True)
    y_hess = jacobian(y_adj, vertcat(u, x))

    ## generate C code
    suffix_name = '_fun'
    fun_name = cost_name + middle_name + suffix_name
    y_fun = Function( fun_name, [x, u, p], \
            [ cost_expr ])
    y_fun.generate( fun_name, casadi_opts )

    suffix_name = '_fun_jac_ut_xt'
    fun_name = cost_name + middle_name + suffix_name
    y_fun_jac_ut_xt = Function(fun_name, [x, u, p], \
            [ cost_expr, cost_jac_expr ])
    y_fun_jac_ut_xt.generate( fun_name, casadi_opts )

    suffix_name = '_hess'
    fun_name = cost_name + middle_name + suffix_name
    y_hess = Function(fun_name, [x, u, p], [ y_hess ])
    y_hess.generate( fun_name, casadi_opts )

    os.chdir('../..')

    return

