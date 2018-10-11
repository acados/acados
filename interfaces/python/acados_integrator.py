from ctypes import *
import ctypes.util 
import numpy as np
from casadi import *
from os import system

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
		ode_expr = -2*x

		# Form a function and generate C code
		name = 'ode_expr'
		ode_expr = Function(name, [x], [ode_expr], ['x'], ['ode_expr'])
		cname = ode_expr.generate()

		oname = 'model.so'
		system('gcc -fPIC -shared ' + cname + ' -o ' + oname)


		## load model library
		__model = CDLL('model.so')
		self.__model = __model

		fun_ptr = cast(addressof(__model.ode_expr), c_void_p)
		print(fun_ptr)




		## config
		self.config = cast(__acados.sim_config_create( 0 ), c_void_p)
		print(self.config)


		## dims
		self.dims = cast(__acados.sim_dims_create(self.config), c_void_p)
		print(self.dims)
		__acados.sim_config_set_nx(self.config, self.dims, nx)
		__acados.sim_config_set_nu(self.config, self.dims, nu)


	

	# TODO free stuff !!!!!!!!
	def __del__(self):
		
		self.__acados.sim_config_free(self.config)
		self.__acados.sim_dims_free(self.dims)



