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


from casadi import *
import os

def generate_c_code_constraint( constraint ):

    casadi_version = CasadiMeta.version()
    casadi_opts = dict(mex=False, casadi_int='int', casadi_real='double')

    if  casadi_version not in ('3.4.5', '3.4.0'):
        # old casadi versions
        raise Exception('Please download and install Casadi 3.4.0 to ensure compatibility with acados. Version ' + casadi_version + ' currently in use.')

    # load constraint variables and expression
    x = constraint.x
    u = constraint.u
    r = constraint.r
    z = constraint.z
    p = constraint.p
    # nc = nh or np 
    nh = constraint.nh 
    nr = constraint.nr
    con_h_expr = constraint.con_h_expr
    con_r_expr = constraint.con_r_expr
    con_name = constraint.name

    # get dimensions
    nx = x.size()[0]
    nu = u.size()[0]
    if r is not None:
        nr = r.size()[0]
    else:
        nr = 0

    if type(p) is list:
        # check that z is empty
        if len(p) == 0:
            np = 0
            p = SX.sym('p', 0, 0)
        else:
            raise Exception('p is a non-empty list. It should be either an empty list or an SX object.')
    else:
        np = p.size()[0]

    if type(z) is list:
        # check that z is empty
        if len(z) == 0:
            nz = 0
            z = SX.sym('z', 0, 0)
        else:
            raise Exception('z is a non-empty list. It should be either an empty list or an SX object.')
    else:
        nz = z.size()[0]

    # set up functions to be exported
    fun_name = con_name + '_h_constraint'
    if nr == 0: # BGH constraint
        jac_x = jacobian(con_h_expr, x);
        jac_u = jacobian(con_h_expr, u);
        jac_z = jacobian(con_h_expr, z);
        constraint_fun_jac_tran = Function(fun_name, [x, u, z, p], [con_h_expr, vertcat(transpose(jac_u), transpose(jac_x)), transpose(jac_z)])

        # generate C code
        if not os.path.exists('c_generated_code'):
            os.mkdir('c_generated_code')

        os.chdir('c_generated_code')
        gen_dir = con_name + '_h_constraint'
        if not os.path.exists(gen_dir):
            os.mkdir(gen_dir)
        gen_dir_location = './' + gen_dir
        os.chdir(gen_dir_location)
        file_name = con_name + '_h_constraint'
        constraint_fun_jac_tran.generate(file_name, casadi_opts)
        os.chdir('../..')
    else: # BGHP constraint
        con_h_expr_x_u = substitute(con_h_expr, r, con_r_expr)
        jac_u = jacobian(con_h_expr_x_u, u)
        jac_x = jacobian(con_h_expr_x_u, x)

        hess = hessian(con_h_expr[0], r)[0]
        for i in range(1, nh):
            hess = vertcat(hess, hessian(con_h_expr[i], r)[0])

        constraint_fun_jac_tran_hess = Function(fun_name, [x, u, z, p], [con_h_expr_x_u, vertcat(transpose(jac_u), transpose(jac_x)), hess])

        # generate C code
        if not os.path.exists('c_generated_code'):
            os.mkdir('c_generated_code')

        os.chdir('c_generated_code')
        gen_dir = con_name + '_h_constraint'
        if not os.path.exists(gen_dir):
            os.mkdir(gen_dir)
        gen_dir_location = './' + gen_dir
        os.chdir(gen_dir_location)
        file_name = con_name + '_h_constraint'
        constraint_fun_jac_tran_hess.generate(file_name, casadi_opts)
        os.chdir('..')


        jac_u = jacobian(con_r_expr, u);
        jac_x = jacobian(con_r_expr, x);
        fun_name = con_name + '_p_constraint'
        constraint_residual_fun_jac_tran = Function(fun_name, [x, u, z, p], [con_r_expr, vertcat(transpose(jac_u), transpose(jac_x))])

        gen_dir = con_name + '_p_constraint'
        if not os.path.exists(gen_dir):
            os.mkdir(gen_dir)
        gen_dir_location = './' + gen_dir
        os.chdir(gen_dir_location)
        file_name = con_name + '_p_constraint'
        constraint_residual_fun_jac_tran.generate(file_name, casadi_opts)

        os.chdir('../..')

    return
