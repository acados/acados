from ctypes import *
import ctypes.util 
import numpy as np
from casadi import *
from os import system

from generate_wrapper import set_function_pointers

#import faulthandler

#faulthandler.enable()





class acados_integrator_model:



	def __init__(self):
		
		self.type = 'explicit'
		self.model_name = 'model'
		self.x = SX.sym('x', 0, 1)
		self.u = SX.sym('u', 0, 1)
		self.xdot = SX.sym('xdot', 0, 1)
	


	def set(self, field, value):

		if field=='type':
			self.type = value

		elif field=='ode_expr':
			self.ode_expr = value

		elif field=='x':
			self.x = value

		elif field=='u':
			self.u = value

		elif field=='xdot':
			self.xdot = value

		elif field=='model_name':
			self.model_name = value

		else:
			disp('acados_integrator_model.set(): wrong field')





class acados_integrator_opts:



	def __init__(self):

		self.scheme = 'erk'
		self.sens_forw = 'false'



	def set(self, field, value):
		
		if field=='scheme':
			self.scheme = value
		
		if field=='sens_forw':
			self.sens_forw = value





class acados_integrator:



	def __init__(self, model, opts):
		
		# checks

		if not (model.type=='explicit'):
			print('\nwrong model type value!\n')
			exit()

		if not (opts.scheme=='erk'):
			print('\nwrong opts scheme value!\n')
			exit()

		if not (opts.sens_forw=='false'):
			print('\nwrong opts sens_forw value!\n')
			exit()


		
#		print(CasadiMeta.version())

		# load acados library
		__acados = CDLL('libacados_c.so')
		self.__acados = __acados


		# sizes
		self.nx = model.x.size()[0]
		self.nu = model.u.size()[0]


		# define functions & generate C code
#		if opts.scheme=='erk'
		casadi_opts = dict(casadi_int='int', casadi_real='double')
		c_sources = ' '

		fun_name = 'expl_ode_fun'
		fun = Function(fun_name, [model.x], [model.ode_expr])
		fun.generate(casadi_opts)
		c_sources = c_sources + fun_name + '.c '

#		jac_x_name = 'expl_ode_jac_x'
#		jac_x_expr = jacobian(model.ode_expr, model.x)
#		jac_x = Function(jac_x_name, [model.x], [jac_x_expr])
#		jac_x.generate(casadi_opts)
#		c_sources = c_sources + jac_x_name + '.c '

#		jac_u_name = 'expl_ode_jac_u'
#		jac_u_expr = jacobian(model.ode_expr, model.u)
#		jac_u = Function(jac_u_name, [model.u], [jac_u_expr])
#		jac_u.generate(casadi_opts)
#		c_sources = c_sources + jac_u_name + '.c '

		# create model library
		lib_name = model.model_name + '.so'
		system('gcc -fPIC -shared ' + c_sources + ' -o ' + lib_name)

		## load model library
		self.__model = CDLL(lib_name)

#		self.__model = model.model


		## external function
		ext_fun_struct_size = __acados.external_function_casadi_struct_size()
		ext_fun_struct = cast(create_string_buffer(ext_fun_struct_size), c_void_p)
		self.ext_fun = ext_fun_struct

		# set function pointers
		set_function_pointers(__acados, lib_name, fun_name, self.ext_fun)

		# create external function
		__acados.external_function_casadi_create(self.ext_fun)


		## config
		if opts.scheme=='erk':
			self.config = cast(__acados.sim_config_create( 0 ), c_void_p)

		## dims
		self.dims = cast(__acados.sim_dims_create(self.config), c_void_p)
		__acados.sim_dims_set_nx(self.config, self.dims, self.nx)
		__acados.sim_dims_set_nu(self.config, self.dims, self.nu)

		## opts
		self.opts = cast(__acados.sim_opts_create(self.config, self.dims), c_void_p)
		if opts.sens_forw=='false':
			__acados.sim_opts_set_sens_forw(self.opts, 0)
		else:
			__acados.sim_opts_set_sens_forw(self.opts, 1)

		## sim_in
		self.sim_in = cast(__acados.sim_in_create(self.config, self.dims), c_void_p)
		if opts.scheme=='erk':
			__acados.sim_set_model(self.config, self.sim_in, "expl_ode_fun", self.ext_fun)

		## sim_out
		self.sim_out = cast(__acados.sim_out_create(self.config, self.dims), c_void_p)

		## sim solver
		self.solver = cast(__acados.sim_create(self.config, self.dims, self.opts), c_void_p)



	def set(self, field, value):

		if field=='x':
			tmp = np.ascontiguousarray(value, dtype=np.float64)
			tmp = cast(tmp.ctypes.data, POINTER(c_double))
			self.__acados.sim_in_set_x(self.config, self.dims, tmp, self.sim_in)

		if field=='t':
			self.__acados.sim_in_set_T(self.config, c_double(value), self.sim_in)
			


	def solve(self):

		# solve
		flag = self.__acados.sim_solve(self.solver, self.sim_in, self.sim_out)

		return flag
	


	def get(self, field):
		
		if field=='xn':
			value = np.zeros((self.nx, 1))
			tmp = cast(value.ctypes.data, POINTER(c_double))
			self.__acados.sim_out_get_xn(self.config, self.dims, self.sim_out, tmp)

		return value



	def __del__(self):
		
		self.__acados.external_function_casadi_free(self.ext_fun)
		self.__acados.sim_config_free(self.config)
		self.__acados.sim_dims_free(self.dims) # double free ???
		self.__acados.sim_opts_free(self.opts)
		self.__acados.sim_out_free(self.sim_in)
		self.__acados.sim_in_free(self.sim_out)
		self.__acados.sim_free(self.solver)



