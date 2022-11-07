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
import casadi
from .utils import ALLOWED_CASADI_VERSIONS, casadi_version_warning

def generate_c_code_conl_cost(model, cost_name, stage_type, opts):

    casadi_version = casadi.CasadiMeta.version()
    casadi_opts = dict(mex=False, casadi_int='int', casadi_real='double')

    if casadi_version not in (ALLOWED_CASADI_VERSIONS):
        casadi_version_warning(casadi_version)

    x = model.x
    p = model.p

    if isinstance(x, casadi.MX):
        symbol = casadi.MX.sym
    else:
        symbol = casadi.SX.sym

    if stage_type == 'terminal':
        u = symbol('u', 0, 0)

        yref = model.cost_r_in_psi_expr_e
        inner_expr = model.cost_y_expr_e - yref
        outer_expr = model.cost_psi_expr_e
        res_expr = model.cost_r_in_psi_expr_e

        suffix_name_fun = '_conl_cost_e_fun'
        suffix_name_fun_jac_hess = '_conl_cost_e_fun_jac_hess'

    elif stage_type == 'initial':
        u = model.u

        yref = model.cost_r_in_psi_expr_0
        inner_expr = model.cost_y_expr_0 - yref
        outer_expr = model.cost_psi_expr_0
        res_expr = model.cost_r_in_psi_expr_0

        suffix_name_fun = '_conl_cost_0_fun'
        suffix_name_fun_jac_hess = '_conl_cost_0_fun_jac_hess'

    elif stage_type == 'path':
        u = model.u

        yref = model.cost_r_in_psi_expr
        inner_expr = model.cost_y_expr - yref
        outer_expr = model.cost_psi_expr
        res_expr = model.cost_r_in_psi_expr

        suffix_name_fun = '_conl_cost_fun'
        suffix_name_fun_jac_hess = '_conl_cost_fun_jac_hess'

    # set up function names
    fun_name_cost_fun = model.name + suffix_name_fun
    fun_name_cost_fun_jac_hess = model.name + suffix_name_fun_jac_hess

    # set up functions to be exported
    outer_loss_fun = casadi.Function('psi', [res_expr, p], [outer_expr])
    cost_expr = outer_loss_fun(inner_expr, p)

    cost_jac_expr = casadi.jacobian(cost_expr, casadi.vertcat(u, x))

    outer_hess_fun = casadi.Function('inner_hess', [res_expr, p], [casadi.hessian(outer_loss_fun(res_expr, p), res_expr)[0]])
    J_expr = casadi.jacobian(inner_expr, casadi.vertcat(u, x))

    cost_fun = casadi.Function(
        fun_name_cost_fun,
        [x, u, yref, p],
        [cost_expr])

    cost_fun_jac_hess = casadi.Function(
        fun_name_cost_fun_jac_hess,
        [x, u, yref, p],
        [cost_expr, cost_jac_expr.T, J_expr.T @ outer_hess_fun(inner_expr, p) @ J_expr]
    )
    # set up directory
    code_export_dir = opts["code_export_directory"]
    if not os.path.exists(code_export_dir):
        os.makedirs(code_export_dir)

    cwd = os.getcwd()
    os.chdir(code_export_dir)
    gen_dir = cost_name + '_cost'
    if not os.path.exists(gen_dir):
        os.mkdir(gen_dir)
    gen_dir_location = os.path.join('.', gen_dir)
    os.chdir(gen_dir_location)

    # generate C code
    cost_fun.generate(fun_name_cost_fun, casadi_opts)
    cost_fun_jac_hess.generate(fun_name_cost_fun_jac_hess, casadi_opts)

    os.chdir(cwd)

    return

