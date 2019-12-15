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

def generate_c_code_constraint_e( constraint, con_name ):

    casadi_version = CasadiMeta.version()
    casadi_opts = dict(mex=False, casadi_int='int', casadi_real='double')

    if casadi_version not in (ALLOWED_CASADI_VERSIONS):
        msg =  'Please download and install CasADi {} '.format(" or ".join(ALLOWED_CASADI_VERSIONS))
        msg += 'to ensure compatibility with acados.\n'
        msg += 'Version {} currently in use.'.format(casadi_version)
        raise Exception(msg)

    # load constraint variables and expression
    x = constraint.x
    r = constraint.r
    p = constraint.p
    nh = constraint.nh
    nphi = constraint.nphi
    if nh > 0 and nphi > 0:
        raise Exception('cannot have both nh_e and phi_e > 0.')

    if nh > 0 or nphi > 0:

        # get dimensions
        nx = x.size()[0]
        if r is not None:
            nr = r.size()[0]
        else:
            nr = 0

        if x is not None:
            nx = x.size()[0]
        else:
            nx = 0

        if type(p) is list:
            # check that p is empty
            if len(p) == 0:
                np = 0
                p = SX.sym('p', 0, 0)
            else:
                raise Exception('p is a non-empty list. It should be either an empty list or an SX object.')
        else:
            np = p.size()[0]

        # create dummy u
        u = SX.sym('u', 0, 0)
        # create dummy r
        z = SX.sym('z', 0, 0)
        # set up functions to be exported
        fun_name = con_name + '_constr_h_e'
        if nr == 0: # BGH constraint
            con_h_expr = constraint.con_h_expr
            jac_x = jacobian(con_h_expr, x)
            jac_z = jacobian(con_h_expr, z)
            constraint_fun_jac_tran = Function(fun_name, [x, u, z, p], [con_h_expr, transpose(jac_x), transpose(jac_z)])

            # generate C code
            if not os.path.exists('c_generated_code'):
                os.mkdir('c_generated_code')

            os.chdir('c_generated_code')
            gen_dir = con_name + '_h_e_constraint'
            if not os.path.exists(gen_dir):
                os.mkdir(gen_dir)
            gen_dir_location = './' + gen_dir
            os.chdir(gen_dir_location)
            file_name = con_name + '_h_e_constraint'
            constraint_fun_jac_tran.generate(file_name, casadi_opts)
            os.chdir('../..')
        else: # BGP constraint
            con_phi_expr = constraint.con_phi_expr
            con_r_expr = constraint.con_r_expr
            gen_dir = con_name + '_phi_e_constraint'
            fun_name = con_name + '_phi_e_constraint'
            con_phi_expr_x = substitute(con_phi_expr, r, con_r_expr)
            phi_jac_x = jacobian(con_phi_expr_x, x)
            phi_jac_u = jacobian(con_phi_expr_x, u)
            phi_jac_z = jacobian(con_phi_expr_x, z)

            r_jac_x = jacobian(con_r_expr, x)
            r_jac_u = jacobian(con_r_expr, u)

            hess = hessian(con_phi_expr[0], r)[0]
            for i in range(1, nh):
                hess = vertcat(hess, hessian(con_phi_expr[i], r)[0])

            constraint_phi = Function(fun_name, [x, u, z, p], [con_phi_expr_x, transpose(phi_jac_x), transpose(phi_jac_z), hess, transpose(r_jac_x)])

            # generate C code
            if not os.path.exists('c_generated_code'):
                os.mkdir('c_generated_code')

            os.chdir('c_generated_code')
            gen_dir = con_name + '_phi_e_constraint'
            if not os.path.exists(gen_dir):
                os.mkdir(gen_dir)
            gen_dir_location = './' + gen_dir
            os.chdir(gen_dir_location)
            file_name = con_name + '_phi_e_constraint'
            constraint_phi.generate(file_name, casadi_opts)
            os.chdir('../..')

            # jac_x = jacobian(con_r_expr, x);
            # fun_name = con_name + '_r_e_constraint'
            # constraint_residual_fun_jac_tran = Function(fun_name, [x, u, z, p], [con_r_expr, transpose(jac_x)])

            # gen_dir = con_name + '_r_e_constraint'
            # if not os.path.exists(gen_dir):
            #     os.mkdir(gen_dir)
            # gen_dir_location = './' + gen_dir
            # os.chdir(gen_dir_location)
            # file_name = con_name + '_r_e_constraint'
            # constraint_residual_fun_jac_tran.generate(file_name, casadi_opts)

            # Jr_tran_eval = (constraint_residual_fun_jac_tran([0,0], [])[1]).full()
            # hess_eval = (constraint_fun_jac_tran_hess([0,0], [])[2]).full()
            # tmp1 = mtimes(Jr_tran_eval, hess_eval[0:3,0:3])
            # HESS1 = mtimes(mtimes(Jr_tran_eval, hess_eval[0:3,0:3]), transpose(Jr_tran_eval))
            # tmp2 = mtimes(Jr_tran_eval, hess_eval[3:6,0:3])
            # HESS2= HESS1+ mtimes(mtimes(Jr_tran_eval, hess_eval[3:6, 0:3]), transpose(Jr_tran_eval))
            # import pdb; pdb.set_trace()
            # os.chdir('../..')

    return
