from ctypes import *
import ctypes.util 
import numpy as np
from casadi import *
from os import system

from generate_wrapper import set_function_pointers

#import faulthandler

#faulthandler.enable()



#class acados_integrator_model:
#	def __init__(self):
#		
#		self.type = 'explicit'
#	
#	def set(self, field, value):
#		if field=='ode_expr':
#			# call casadi codgen
#		if field=='x':
#			# set casadi variable x
#		if field=='u':
#			# set casadi variable u
#		if field=='xdot':
#			# set casadi variable xdot



#class acados_integrator_opts:
#	def __init__(self, model):
#		
#		# default
#		self.ns = 4
#		self.num_steps = 1
#		
#		if model.type=='explicit':
#			self.type = 'erk'
#	
#	def set(self, field, value):
#		if field=='ns':
#			self.ns = value
#		if field=='num_steps':
#			self.num_steps = value



class acados_integrator:
#	def __init__(self, opts, model):
	def __init__(self):
		
#		print(CasadiMeta.version())

		# load acados library
		__acados = CDLL('libacados_c.so')
		self.__acados = __acados


		# nx
		nx = 4
		nu = 1

		x = SX.sym('x', nx, 1)
		casadi_ode_expr = -2*x

		# Form a function and generate C code
		user_fun_name = 'ode_expr'
		python_ode_expr = Function(user_fun_name, [x], [casadi_ode_expr], ['x'], ['ode_expr'])
		cname = python_ode_expr.generate()

		model_name = 'model.so'
		system('gcc -fPIC -shared '+user_fun_name+'.c -o ' + model_name)


		## load model library
		__model = CDLL(model_name)
		self.__model = __model



		## external function
		ext_fun_struct_size = __acados.external_function_casadi_struct_size()
		ext_fun_struct = cast(create_string_buffer(ext_fun_struct_size), c_void_p)
		self.ext_fun = ext_fun_struct

		# set function pointers
		set_function_pointers(__acados, model_name, user_fun_name, self.ext_fun)

		# create external function
		__acados.external_function_casadi_create(self.ext_fun)



		## config
		self.config = cast(__acados.sim_config_create( 0 ), c_void_p)
		print(self.config)



		## dims
		self.dims = cast(__acados.sim_dims_create(self.config), c_void_p)
		print(self.dims)
		__acados.sim_config_set_nx(self.config, self.dims, nx)
		__acados.sim_config_set_nu(self.config, self.dims, nu)



		## opts
		self.opts = cast(__acados.sim_opts_create(self.config), c_void_p)
		print(self.opts)


	

	# TODO free stuff !!!!!!!!
	def __del__(self):
		
		self.__acados.external_function_casadi_free(self.ext_fun)
		self.__acados.sim_config_free(self.config)
		self.__acados.sim_dims_free(self.dims)
#		self.__acados.sim_opts_free(self.opts)



