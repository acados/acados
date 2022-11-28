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
# SUBSTITUTE GOODS OR SERVICES LOSS OF USE, DATA, OR PROFITS OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#

import os
import casadi as ca
from .utils import check_casadi_version, casadi_length

def generate_c_code_discrete_dynamics( model, opts ):

    check_casadi_version()

    casadi_codegen_opts = dict(mex=False, casadi_int='int', casadi_real='double')

    # load model
    x = model.x
    u = model.u
    p = model.p
    phi = model.disc_dyn_expr
    model_name = model.name
    nx = casadi_length(x)

    if isinstance(phi, ca.MX):
        symbol = ca.MX.sym
    elif isinstance(phi, ca.SX):
        symbol = ca.SX.sym
    else:
        raise Exception("generate_c_code_disc_dyn: disc_dyn_expr must be a CasADi expression, you have type: {}".format(type(phi)))

    # assume nx1 = nx !!!
    lam = symbol('lam', nx, 1)

    # generate jacobians
    ux = ca.vertcat(u,x)
    jac_ux = ca.jacobian(phi, ux)
    # generate adjoint
    adj_ux = ca.jtimes(phi, ux, lam, True)
    # generate hessian
    hess_ux = ca.jacobian(adj_ux, ux)

    # change directory
    cwd = os.getcwd()
    model_dir = os.path.abspath(os.path.join(opts["code_export_directory"], f'{model_name}_model'))
    if not os.path.exists(model_dir):
        os.makedirs(model_dir)
    os.chdir(model_dir)

    # set up & generate Functions
    fun_name = model_name + '_dyn_disc_phi_fun'
    phi_fun = ca.Function(fun_name, [x, u, p], [phi])
    phi_fun.generate(fun_name, casadi_codegen_opts)

    fun_name = model_name + '_dyn_disc_phi_fun_jac'
    phi_fun_jac_ut_xt = ca.Function(fun_name, [x, u, p], [phi, jac_ux.T])
    phi_fun_jac_ut_xt.generate(fun_name, casadi_codegen_opts)

    fun_name = model_name + '_dyn_disc_phi_fun_jac_hess'
    phi_fun_jac_ut_xt_hess = ca.Function(fun_name, [x, u, lam, p], [phi, jac_ux.T, hess_ux])
    phi_fun_jac_ut_xt_hess.generate(fun_name, casadi_codegen_opts)

    os.chdir(cwd)
