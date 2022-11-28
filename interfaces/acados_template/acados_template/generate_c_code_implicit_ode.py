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
from casadi import SX, MX, Function, jacobian, vertcat, jtimes
from .utils import casadi_length, check_casadi_version

def generate_c_code_implicit_ode( model, opts ):

    check_casadi_version()

    casadi_codegen_opts = dict(mex=False, casadi_int='int', casadi_real='double')

    # load model
    x = model.x
    xdot = model.xdot
    u = model.u
    z = model.z
    p = model.p
    f_impl = model.f_impl_expr
    model_name = model.name

    # get model dimensions
    nx = casadi_length(x)
    nz = casadi_length(z)

    # generate jacobians
    jac_x       = jacobian(f_impl, x)
    jac_xdot    = jacobian(f_impl, xdot)
    jac_u       = jacobian(f_impl, u)
    jac_z       = jacobian(f_impl, z)

    # Set up functions
    p = model.p
    fun_name = model_name + '_impl_dae_fun'
    impl_dae_fun = Function(fun_name, [x, xdot, u, z, p], [f_impl])

    fun_name = model_name + '_impl_dae_fun_jac_x_xdot_z'
    impl_dae_fun_jac_x_xdot_z = Function(fun_name, [x, xdot, u, z, p], [f_impl, jac_x, jac_xdot, jac_z])

    fun_name = model_name + '_impl_dae_fun_jac_x_xdot_u_z'
    impl_dae_fun_jac_x_xdot_u_z = Function(fun_name, [x, xdot, u, z, p], [f_impl, jac_x, jac_xdot, jac_u, jac_z])

    fun_name = model_name + '_impl_dae_fun_jac_x_xdot_u'
    impl_dae_fun_jac_x_xdot_u = Function(fun_name, [x, xdot, u, z, p], [f_impl, jac_x, jac_xdot, jac_u])

    fun_name = model_name + '_impl_dae_jac_x_xdot_u_z'
    impl_dae_jac_x_xdot_u_z = Function(fun_name, [x, xdot, u, z, p], [jac_x, jac_xdot, jac_u, jac_z])

    if opts["generate_hess"]:
        x_xdot_z_u = vertcat(x, xdot, z, u)
        if isinstance(x, MX):
            symbol = MX.sym
        else:
            symbol = SX.sym
        multiplier = symbol('multiplier', nx + nz)
        ADJ = jtimes(f_impl, x_xdot_z_u, multiplier, True)
        HESS = jacobian(ADJ, x_xdot_z_u)
        fun_name = model_name + '_impl_dae_hess'
        impl_dae_hess = Function(fun_name, [x, xdot, u, z, multiplier, p], [HESS])

    # change directory
    cwd = os.getcwd()
    model_dir = os.path.abspath(os.path.join(opts["code_export_directory"], f'{model_name}_model'))
    if not os.path.exists(model_dir):
        os.makedirs(model_dir)
    os.chdir(model_dir)

    # generate C code
    fun_name = model_name + '_impl_dae_fun'
    impl_dae_fun.generate(fun_name, casadi_codegen_opts)

    fun_name = model_name + '_impl_dae_fun_jac_x_xdot_z'
    impl_dae_fun_jac_x_xdot_z.generate(fun_name, casadi_codegen_opts)

    fun_name = model_name + '_impl_dae_jac_x_xdot_u_z'
    impl_dae_jac_x_xdot_u_z.generate(fun_name, casadi_codegen_opts)

    fun_name = model_name + '_impl_dae_fun_jac_x_xdot_u_z'
    impl_dae_fun_jac_x_xdot_u_z.generate(fun_name, casadi_codegen_opts)

    fun_name = model_name + '_impl_dae_fun_jac_x_xdot_u'
    impl_dae_fun_jac_x_xdot_u.generate(fun_name, casadi_codegen_opts)

    if opts["generate_hess"]:
        fun_name = model_name + '_impl_dae_hess'
        impl_dae_hess.generate(fun_name, casadi_codegen_opts)

    os.chdir(cwd)
