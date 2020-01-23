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
from .utils import ALLOWED_CASADI_VERSIONS, is_empty, casadi_length

def generate_c_code_constraint( model, con_name ):

    casadi_version = CasadiMeta.version()
    casadi_opts = dict(mex=False, casadi_int='int', casadi_real='double')

    if casadi_version not in (ALLOWED_CASADI_VERSIONS):
        msg =  'Please download and install CasADi {} '.format(" or ".join(ALLOWED_CASADI_VERSIONS))
        msg += 'to ensure compatibility with acados.\n'
        msg += 'Version {} currently in use.'.format(casadi_version)
        raise Exception(msg)

    # load constraint variables and expression
    x = model.x
    u = model.u
    r = model.r
    z = model.z
    p = model.p

    con_h_expr = model.con_h_expr
    con_phi_expr = model.con_phi_expr

    if (not is_empty(con_h_expr)) and (not is_empty(con_phi_expr)):
        raise Exception("acados: you can either have constraint_h, or constraint_phi, not both.")

    if not (is_empty(con_h_expr) and is_empty(con_phi_expr)):
        if is_empty(con_h_expr):
            constr_type = 'BGP'
        else:
            constr_type = 'BGH'

        # get dimensions
        if x is not None:
            nx = x.size()[0]
        else:
            nx = 0

        if u is not None:
            nu = u.size()[0]
        else:
            nu = 0

        if r is not None:
            nr = r.size()[0]
        else:
            nr = 0

        if is_empty(p):
            p = SX.sym('p', 0, 0)

        if is_empty(z):
            z = SX.sym('z', 0, 0)

        # set up & change directory
        if not os.path.exists('c_generated_code'):
            os.mkdir('c_generated_code')
        os.chdir('c_generated_code')
        gen_dir = con_name + '_constraints'
        if not os.path.exists(gen_dir):
            os.mkdir(gen_dir)
        gen_dir_location = './' + gen_dir
        os.chdir(gen_dir_location)

        # export casadi functions
        if constr_type == 'BGH':
            fun_name = con_name + '_constr_h_fun_jac_uxt_zt'
            file_name = con_name + '_constr_h_fun_jac_uxt_zt'
            jac_x = jacobian(con_h_expr, x)
            jac_u = jacobian(con_h_expr, u)
            jac_z = jacobian(con_h_expr, z)
            constraint_fun_jac_tran = \
                Function(fun_name, [x, u, z, p], \
                [con_h_expr, vertcat(transpose(jac_u), \
                transpose(jac_x)), transpose(jac_z)])

            constraint_fun_jac_tran.generate(file_name, casadi_opts)

        else: # BGP constraint
            nphi = casadi_length(con_phi_expr)
            con_r_expr = model.con_r_expr
            fun_name = con_name + '_phi_constraint'
            con_phi_expr_x_u_z = substitute(con_phi_expr, r, con_r_expr)
            phi_jac_u = jacobian(con_phi_expr_x_u_z, u)
            phi_jac_x = jacobian(con_phi_expr_x_u_z, x)
            phi_jac_z = jacobian(con_phi_expr_x_u_z, z)

            hess = hessian(con_phi_expr[0], r)[0]
            for i in range(1, nphi):
                hess = vertcat(hess, hessian(con_phi_expr[i], r)[0])

            r_jac_u = jacobian(con_r_expr, u)
            r_jac_x = jacobian(con_r_expr, x)

            constraint_phi = \
                Function(fun_name, [x, u, z, p], \
                [con_phi_expr_x_u_z, \
                vertcat(transpose(phi_jac_u), \
                transpose(phi_jac_x)), \
                transpose(phi_jac_z), \
                hess, vertcat(transpose(r_jac_u), \
                transpose(r_jac_x))])

            file_name = con_name + '_phi_constraint'
            constraint_phi.generate(file_name, casadi_opts)

        # change directory back
        os.chdir('../..')

    return
