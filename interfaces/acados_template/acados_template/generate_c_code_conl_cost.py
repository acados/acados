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
from .utils import check_casadi_version

def generate_c_code_conl_cost(model, cost_name, stage_type, opts):

    check_casadi_version()

    casadi_codegen_opts = dict(mex=False, casadi_int='int', casadi_real='double')

    x = model.x
    z = model.z
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

        custom_hess = model.cost_conl_custom_outer_hess_e

    elif stage_type == 'initial':
        u = model.u

        yref = model.cost_r_in_psi_expr_0
        inner_expr = model.cost_y_expr_0 - yref
        outer_expr = model.cost_psi_expr_0
        res_expr = model.cost_r_in_psi_expr_0

        suffix_name_fun = '_conl_cost_0_fun'
        suffix_name_fun_jac_hess = '_conl_cost_0_fun_jac_hess'

        custom_hess = model.cost_conl_custom_outer_hess_0

    elif stage_type == 'path':
        u = model.u

        yref = model.cost_r_in_psi_expr
        inner_expr = model.cost_y_expr - yref
        outer_expr = model.cost_psi_expr
        res_expr = model.cost_r_in_psi_expr

        suffix_name_fun = '_conl_cost_fun'
        suffix_name_fun_jac_hess = '_conl_cost_fun_jac_hess'

        custom_hess = model.cost_conl_custom_outer_hess

    # set up function names
    fun_name_cost_fun = model.name + suffix_name_fun
    fun_name_cost_fun_jac_hess = model.name + suffix_name_fun_jac_hess

    # set up functions to be exported
    outer_loss_fun = casadi.Function('psi', [res_expr, p], [outer_expr])
    cost_expr = outer_loss_fun(inner_expr, p)

    outer_loss_grad_fun = casadi.Function('outer_loss_grad', [res_expr, p], [casadi.jacobian(outer_expr, res_expr).T])

    if custom_hess is None:
        outer_hess_fun = casadi.Function('inner_hess', [res_expr, p], [casadi.hessian(outer_loss_fun(res_expr, p), res_expr)[0]])
    else:
        outer_hess_fun = casadi.Function('inner_hess', [res_expr, p], [custom_hess])

    Jt_ux_expr = casadi.jacobian(inner_expr, casadi.vertcat(u, x)).T
    Jt_z_expr = casadi.jacobian(inner_expr, z).T

    cost_fun = casadi.Function(
        fun_name_cost_fun,
        [x, u, z, yref, p],
        [cost_expr])

    cost_fun_jac_hess = casadi.Function(
        fun_name_cost_fun_jac_hess,
        [x, u, z, yref, p],
        [cost_expr, outer_loss_grad_fun(inner_expr, p), Jt_ux_expr, Jt_z_expr, outer_hess_fun(inner_expr, p)]
    )
    # change directory
    cwd = os.getcwd()
    cost_dir = os.path.abspath(os.path.join(opts["code_export_directory"], f'{model.name}_cost'))
    if not os.path.exists(cost_dir):
        os.makedirs(cost_dir)
    os.chdir(cost_dir)

    # generate C code
    cost_fun.generate(fun_name_cost_fun, casadi_codegen_opts)
    cost_fun_jac_hess.generate(fun_name_cost_fun_jac_hess, casadi_codegen_opts)

    os.chdir(cwd)

    return

