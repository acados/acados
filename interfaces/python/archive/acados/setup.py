from setuptools import setup, find_packages

import acados

setup(name='acados',
   version='0.1',
   description='Python interface to acados',
   url='http://github.com/giaf/hpipm',
   author='Gianluca Frison - Andrea Zanelli',
   license='GPL+CE',
   packages = find_packages(),
   zip_safe=False)
