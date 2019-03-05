from setuptools import setup, find_packages

import acados_template

setup(name='acados_template',
   version='0.1',
   description='a templating framework for acados',
   url='http://github.com/zanellia/acados',
   author='Andrea Zanelli',
   license='LGPL',
   packages = find_packages(),
   include_package_data = True,
   package_data={'': ['c_templates/main.in.c', 
       'c_templates/Makefile.in', 
       'c_templates/model.in.h',
       'c_templates/acados_solver.in.c', 
       'c_templates/acados_solver_sfun.in.c', 
       'c_templates/acados_solver.in.h',
       'c_templates/p_constraint.in.h',
       'c_templates/h_constraint.in.h',
       'c_templates/make_sfun.in.m', 
       ]},
   zip_safe=False)
