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

def generate_c_code_constraint_e( constraint, suffix_name ):

    casadi_version = CasadiMeta.version()
    casadi_opts = dict(mex=False, casadi_int='int', casadi_real='double')

    if  casadi_version not in ('3.4.5', '3.4.0'):
        # old casadi versions
        raise Exception('Please download and install Casadi 3.4.0 to ensure compatibility with acados. Version ' + casadi_version + ' currently in use.')

    # load constraint variables and expression
    x = constraint.x
    u = constraint.u
    # nc = nh or np 
    nc = constraint.nc 
    con_exp = constraint.expr
    con_name = constraint.name

    # get dimensions
    nx = x.size()[0]
    nu = u.size()[0]

    # set up functions to be exported
    fun_name = con_name + suffix_name
    # TODO(andrea): first output seems to be ignored in the C code
    jac_x = jacobian(con_exp, x);
    # jac_u = jacobian(con_exp, u);
    constraint_fun_jac_tran = Function(fun_name, [x, u], [con_exp, transpose(jac_x)])

    # generate C code
    if not os.path.exists('c_generated_code'):
        os.mkdir('c_generated_code')

    os.chdir('c_generated_code')
    gen_dir = con_name + suffix_name 
    if not os.path.exists(gen_dir):
        os.mkdir(gen_dir)
    gen_dir_location = './' + gen_dir
    os.chdir(gen_dir_location)
    file_name = con_name + suffix_name
    constraint_fun_jac_tran.generate(file_name, casadi_opts)
    os.chdir('../..')

    return
