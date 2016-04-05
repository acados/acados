#!/usr/bin/env python

"""
setup.py file for SWIG example
"""

from distutils.core import setup, Extension


erk_integrator_module = Extension('_erk_integrator',
                           sources=['erk_integrator_wrap.c', 'erk_integrator.c', 'auxiliary_functions.c', 'model.c', 'timing_functions.c'],
                           )

setup (name = 'erk_integrator',
       version = '0.1',
       author      = "SWIG Docs",
       description = """Simple swig example from docs""",
       ext_modules = [erk_integrator_module],
       py_modules = ["erk_integrator"],
       )
