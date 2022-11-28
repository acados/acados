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
from .utils import casadi_length, check_casadi_version

def generate_c_code_nls_cost( model, cost_name, stage_type, opts ):

    check_casadi_version()

    casadi_codegen_opts = dict(mex=False, casadi_int='int', casadi_real='double')

    x = model.x
    z = model.z
    p = model.p
    u = model.u

    if isinstance(x, casadi.MX):
        symbol = MX.sym
    else:
        symbol = SX.sym

    if stage_type == 'terminal':
        middle_name = '_cost_y_e'
        u = symbol('u', 0, 0)
        y_expr = model.cost_y_expr_e

    elif stage_type == 'initial':
        middle_name = '_cost_y_0'
        y_expr = model.cost_y_expr_0

    elif stage_type == 'path':
        middle_name = '_cost_y'
        y_expr = model.cost_y_expr

    # change directory
    cwd = os.getcwd()
    cost_dir = os.path.abspath(os.path.join(opts["code_export_directory"], f'{model.name}_cost'))
    if not os.path.exists(cost_dir):
        os.makedirs(cost_dir)
    os.chdir(cost_dir)

    # set up expressions
    cost_jac_expr = transpose(jacobian(y_expr, vertcat(u, x)))
    dy_dz = jacobian(y_expr, z)
    ny = casadi_length(y_expr)

    y = symbol('y', ny, 1)

    y_adj = jtimes(y_expr, vertcat(u, x), y, True)
    y_hess = jacobian(y_adj, vertcat(u, x))

    ## generate C code
    suffix_name = '_fun'
    fun_name = cost_name + middle_name + suffix_name
    y_fun = Function( fun_name, [x, u, z, p], [ y_expr ])
    y_fun.generate( fun_name, casadi_codegen_opts )

    suffix_name = '_fun_jac_ut_xt'
    fun_name = cost_name + middle_name + suffix_name
    y_fun_jac_ut_xt = Function(fun_name, [x, u, z, p], [ y_expr, cost_jac_expr, dy_dz ])
    y_fun_jac_ut_xt.generate( fun_name, casadi_codegen_opts )

    suffix_name = '_hess'
    fun_name = cost_name + middle_name + suffix_name
    y_hess = Function(fun_name, [x, u, z, y, p], [ y_hess ])
    y_hess.generate( fun_name, casadi_codegen_opts )

    os.chdir(cwd)

    return

