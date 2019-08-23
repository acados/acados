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

from setuptools import setup, find_packages

import sys
print(sys.version_info)
if sys.version_info < (3,5) or sys.version_info > (3,7):
    sys.exit('3.5 <= Python version < 3.7 required. Exiting.')


setup(name='acados_template',
   version='0.1',
   python_requires='>=3.5, <3.7',
   description='A templating framework for acados',
   url='http://github.com/zanellia/acados',
   author='Andrea Zanelli',
   license='LGPL',
   packages = find_packages(),
   include_package_data = True,
   install_requires=[
      'jinja2',
      'numpy',
      'scipy',
      'casadi==3.4.0'
   ],
   package_data={'': [
       'c_templates/main.in.c',
       'c_templates/Makefile.in',
       'c_templates/model.in.h',
       'c_templates/main.in.h',
       'c_templates/acados_solver.in.c',
       'c_templates/acados_solver.in.h',
       'c_templates/acados_sim_solver.in.c',
       'c_templates/acados_sim_solver.in.h',
       'c_templates/acados_solver_sfun.in.c',
       'c_templates/p_constraint.in.h',
       'c_templates/h_constraint.in.h',
       'c_templates/make_sfun.in.m',
       'acados_layout.json'
       ]},
   zip_safe=False
)
