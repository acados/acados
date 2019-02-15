from ctypes import *
import ctypes.util 
import _ctypes
import numpy as np
from os import system
import sys

from acados.sim.generate_wrapper import set_function_pointers

#import faulthandler

#faulthandler.enable()





class acados_integrator_model:



	def __init__(self):
		
		self.model_name = 'model'
		self.type = 'explicit'
		self.ode_expr = None
		self.x = None
		self.u = None
		self.xdot = None
		self.z = None
		self.nx = 0
		self.nu = 0
		self.nz = 0
		self.ode_expr_hash = None
	


	def set(self, field, value):

		if field=='model_name':
			self.model_name = value

		elif field=='type':
			self.type = value

		elif field=='ode_expr':
			self.ode_expr = value
			self.ode_expr_hash = hash(str(self.ode_expr))

		elif field=='ode_expr_hash':
			self.ode_expr_hash = value

		elif field=='x':
			self.x = value
			self.nx = self.x.size()[0]

		elif field=='u':
			self.u = value
			self.nu = self.u.size()[0]

		elif field=='xdot':
			self.xdot = value

		elif field=='z':
			self.z = value
			self.nz = self.z.size()[0]

		elif field=='nx':
			self.nx = value

		elif field=='nu':
			self.nu = value

		elif field=='nz':
			self.nz = value

		else:
			print('acados_integrator_model.set(): wrong field')





class acados_integrator_opts:



	def __init__(self):

		self.scheme = 'erk'
		self.sens_forw = 'false'
		self.codgen_model = 'true'



	def set(self, field, value):
		
		if field=='scheme':
			self.scheme = value
		
		if field=='sens_forw':
			self.sens_forw = value

		if field=='codgen_model':
			self.codgen_model = value





class acados_integrator:



	def __init__(self, model, opts):
		
#		print(id(self))
		
#		print(model.ode_expr)
#		print(str(model.ode_expr))
#		print(hash(str(model.ode_expr)))

		# checks

		if not ((model.type=='explicit') | (model.type=='implicit')):
			print('\nwrong model type value!\n')
			exit()

		if not ((opts.scheme=='erk') | (opts.scheme=='irk')):
			print('\nwrong opts scheme value!\n')
			exit()

		if not ((opts.sens_forw=='true') | (opts.sens_forw=='false')):
			print('\nwrong opts sens_forw value!\n')
			exit()

		if ((model.type=='explicit') & (opts.scheme!='erk')):
			print('\nwrong opts scheme value for explicit model!\n')
			exit()

		if ((model.type=='implicit') & (opts.scheme!='irk')):
			print('\nwrong opts scheme value for implicit model!\n')
			exit()


		self.nx = model.nx
		self.nu = model.nu
		self.nz = model.nz
		self.scheme = opts.scheme
		self.sens_forw = opts.sens_forw

		
#		print(CasadiMeta.version())

		# load acados library
		__acados = CDLL('libacados_c.so')
		self.__acados = __acados


		# codgen model library
		if opts.codgen_model=='true':
			self.codgen_model(model, opts)


		## load model library
		lib_name = model.model_name

#		lib_name = lib_name + '_' + str(id(self))

		if opts.scheme=='erk':
			lib_name = lib_name + '_erk'
		elif opts.scheme=='irk':
			lib_name = lib_name + '_irk'

		if opts.sens_forw=='false':
			lib_name = lib_name + '_0'
		else:
			lib_name = lib_name + '_1'

		lib_name = lib_name + '_' + str(model.ode_expr_hash)

		lib_name = lib_name + '.so'

#		print(self.__model._handle)
		self.__model = CDLL(lib_name)
#		print(self.__model._handle)



#		self.__model = model.model

		## external function
		ext_fun_struct_size = __acados.external_function_casadi_struct_size()

		if opts.scheme=='erk':

			if opts.sens_forw=='false':

				fun_name = 'expl_ode_fun'
				ext_fun_struct = cast(create_string_buffer(ext_fun_struct_size), c_void_p)
				self.expl_ode_fun = ext_fun_struct
				# set function pointers
				set_function_pointers(__acados, lib_name, fun_name, self.expl_ode_fun)
				# create external function
				__acados.external_function_casadi_create(self.expl_ode_fun)

			else:

				fun_name = 'expl_vde_for'
				ext_fun_struct = cast(create_string_buffer(ext_fun_struct_size), c_void_p)
				self.expl_vde_for = ext_fun_struct
				# set function pointers
				set_function_pointers(__acados, lib_name, fun_name, self.expl_vde_for)
				# create external function
				__acados.external_function_casadi_create(self.expl_vde_for)

		elif opts.scheme=='irk':

			fun_name = 'impl_ode_fun'
			ext_fun_struct = cast(create_string_buffer(ext_fun_struct_size), c_void_p)
			self.impl_ode_fun = ext_fun_struct
			# set function pointers
			set_function_pointers(__acados, lib_name, fun_name, self.impl_ode_fun)
			# create external function
			__acados.external_function_casadi_create(self.impl_ode_fun)

			fun_name = 'impl_ode_fun_jac_x_xdot_z'
			ext_fun_struct = cast(create_string_buffer(ext_fun_struct_size), c_void_p)
			self.impl_ode_fun_jac_x_xdot_z = ext_fun_struct
			# set function pointers
			set_function_pointers(__acados, lib_name, fun_name, self.impl_ode_fun_jac_x_xdot_z)
			# create external function
			__acados.external_function_casadi_create(self.impl_ode_fun_jac_x_xdot_z)

			if opts.sens_forw=='true':

				fun_name = 'impl_ode_jac_x_xdot_u_z'
				ext_fun_struct = cast(create_string_buffer(ext_fun_struct_size), c_void_p)
				self.impl_ode_jac_x_xdot_u_z = ext_fun_struct
				# set function pointers
				set_function_pointers(__acados, lib_name, fun_name, self.impl_ode_jac_x_xdot_u_z)
				# create external function
				__acados.external_function_casadi_create(self.impl_ode_jac_x_xdot_u_z)


		## config
		__acados.sim_config_create.restype = c_void_p
		if opts.scheme=='erk':
			self.config = cast(__acados.sim_config_create( 0 ), c_void_p)
		elif opts.scheme=='irk':
			self.config = cast(__acados.sim_config_create( 2 ), c_void_p)

		## dims
		__acados.sim_dims_create.restype = c_void_p
		self.dims = cast(__acados.sim_dims_create(self.config), c_void_p)
		__acados.sim_dims_set_nx(self.config, self.dims, self.nx)
		__acados.sim_dims_set_nu(self.config, self.dims, self.nu)

		## opts
		__acados.sim_opts_create.restype = c_void_p
		self.opts = cast(__acados.sim_opts_create(self.config, self.dims), c_void_p)
		if opts.sens_forw=='false':
			__acados.sim_opts_set_sens_forw(self.opts, 0)
		else:
			__acados.sim_opts_set_sens_forw(self.opts, 1)

		## sim_in
		__acados.sim_in_create.restype = c_void_p
		self.sim_in = cast(__acados.sim_in_create(self.config, self.dims), c_void_p)
		flag = 0

		if opts.scheme=='erk':

			if opts.sens_forw=='false':

				ext_fun = 'expl_ode_fun'
				ext_fun_b = ext_fun.encode('utf-8')
				flag = flag + __acados.sim_set_model(self.config, self.sim_in, c_char_p(ext_fun_b), self.expl_ode_fun)

			else:

				ext_fun = 'expl_vde_for'
				ext_fun_b = ext_fun.encode('utf-8')
				flag = flag + __acados.sim_set_model(self.config, self.sim_in, c_char_p(ext_fun_b), self.expl_vde_for)

		if opts.scheme=='irk':

			ext_fun = 'impl_ode_fun'
			ext_fun_b = ext_fun.encode('utf-8')
			flag = flag + __acados.sim_set_model(self.config, self.sim_in, c_char_p(ext_fun_b), self.impl_ode_fun)

			ext_fun = 'impl_ode_fun_jac_x_xdot' #_z'
			ext_fun_b = ext_fun.encode('utf-8')
			flag = flag + __acados.sim_set_model(self.config, self.sim_in, c_char_p(ext_fun_b), self.impl_ode_fun_jac_x_xdot_z)

			if opts.sens_forw=='true':

				ext_fun = 'impl_ode_jac_x_xdot_u' #_z'
				ext_fun_b = ext_fun.encode('utf-8')
				flag = flag + __acados.sim_set_model(self.config, self.sim_in, c_char_p(ext_fun_b), self.impl_ode_jac_x_xdot_u_z)

		if flag!=0:
			print("\nwrong set model name in acados\n")
			exit()

		if opts.sens_forw=='true':
			# set Sx
			Sx0 = np.zeros((self.nx, self.nx), order='F')
			for ii in range(self.nx):
				Sx0[ii][ii] = 1.0
			tmp = np.ascontiguousarray(Sx0, dtype=np.float64)
			tmp = cast(tmp.ctypes.data, POINTER(c_double))
			self.__acados.sim_in_set_Sx(self.config, self.dims, tmp, self.sim_in)
			# set Su
			Su0 = np.zeros((self.nx, self.nu), order='F')
			tmp = np.ascontiguousarray(Su0, dtype=np.float64)
			tmp = cast(tmp.ctypes.data, POINTER(c_double))
			self.__acados.sim_in_set_Su(self.config, self.dims, tmp, self.sim_in)

		## sim_out
		__acados.sim_out_create.restype = c_void_p
		self.sim_out = cast(__acados.sim_out_create(self.config, self.dims), c_void_p)

		## sim solver
		__acados.sim_create.restype = c_void_p
		self.solver = cast(__acados.sim_solver_create(self.config, self.dims, self.opts), c_void_p)

		## unload model library
#		_ctypes.dlclose(self.__model._handle)



	def codgen_model(self, model, opts):

#		from casadi import * # syntax valid only for the entire module
		import casadi

		# x
		if model.x==None:
			x = casadi.SX.sym('x', 0, 1)
		else:
			x = model.x
		# xdot
		if model.xdot==None:
			xdot = casadi.SX.sym('xdot', 0, 1)
		else:
			xdot = model.xdot
		# u
		if model.u==None:
			u = casadi.SX.sym('u', 0, 1)
		else:
			u = model.u
		# z
		if model.z==None:
			z = casadi.SX.sym('z', 0, 1)
		else:
			z = model.z

		# fun
		fun = model.ode_expr

		# sizes
		nx = model.nx
		nu = model.nu
		nz = model.nz

		# define functions & generate C code
		casadi_opts = dict(casadi_int='int', casadi_real='double')
		c_sources = ' '

		if opts.scheme=='erk':

			if opts.sens_forw=='false':

				fun_name = 'expl_ode_fun'
				casadi_fun = casadi.Function(fun_name, [x, u], [fun])
				casadi_fun.generate(casadi_opts)
				c_sources = c_sources + fun_name + '.c '

			else:

				fun_name = 'expl_vde_for'
				Sx = casadi.SX.sym('Sx', nx, nx)
				Su = casadi.SX.sym('Su', nx, nu)
				vde_x = casadi.jtimes(fun, x, Sx)
				vde_u = casadi.jacobian(fun, u) + casadi.jtimes(fun, x, Su)
				casadi_fun = casadi.Function(fun_name, [x, Sx, Su, u], [fun, vde_x, vde_u])
				casadi_fun.generate(casadi_opts)
				c_sources = c_sources + fun_name + '.c '

		elif opts.scheme=='irk':

			fun_name = 'impl_ode_fun'
			casadi_fun = casadi.Function(fun_name, [x, xdot, u, z], [fun])
			casadi_fun.generate(casadi_opts)
			c_sources = c_sources + fun_name + '.c '

			fun_name = 'impl_ode_fun_jac_x_xdot_z'
			jac_x = casadi.jacobian(fun, x)
			jac_xdot = casadi.jacobian(fun, xdot)
			jac_z = casadi.jacobian(fun, z)
			casadi_fun = casadi.Function(fun_name, [x, xdot, u, z], [fun, jac_x, jac_xdot, jac_z])
			casadi_fun.generate(casadi_opts)
			c_sources = c_sources + fun_name + '.c '

			if opts.sens_forw=='true':

				fun_name = 'impl_ode_jac_x_xdot_u_z'
				jac_x = casadi.jacobian(fun, x)
				jac_xdot = casadi.jacobian(fun, xdot)
				jac_u = casadi.jacobian(fun, u)
				jac_z = casadi.jacobian(fun, z)
				casadi_fun = casadi.Function(fun_name, [x, xdot, u, z], [jac_x, jac_xdot, jac_u, jac_z])
				casadi_fun.generate(casadi_opts)
				c_sources = c_sources + fun_name + '.c '

		# create model library
		lib_name = model.model_name

#		lib_name = lib_name + '_' + str(id(self))

		if opts.scheme=='erk':
			lib_name = lib_name + '_erk'
		elif opts.scheme=='irk':
			lib_name = lib_name + '_irk'

		if opts.sens_forw=='false':
			lib_name = lib_name + '_0'
		else:
			lib_name = lib_name + '_1'

		lib_name = lib_name + '_' + str(model.ode_expr_hash)

		lib_name = lib_name + '.so'

		system('gcc -fPIC -shared ' + c_sources + ' -o ' + lib_name)



	def set(self, field, value):

		if field=='x':
			tmp = np.ascontiguousarray(value, dtype=np.float64)
			tmp = cast(tmp.ctypes.data, POINTER(c_double))
			self.__acados.sim_in_set_x(self.config, self.dims, tmp, self.sim_in)

		if field=='xdot':
			tmp = np.ascontiguousarray(value, dtype=np.float64)
			tmp = cast(tmp.ctypes.data, POINTER(c_double))
			self.__acados.sim_in_set_xdot(self.config, self.dims, tmp, self.sim_in)

		if field=='t':
			self.__acados.sim_in_set_T(self.config, c_double(value), self.sim_in)
			


	def solve(self):

		# solve
		flag = self.__acados.sim_solve(self.solver, self.sim_in, self.sim_out)

		## unload model library
#		_ctypes.dlclose(self.__model._handle)

		return flag
	


	def get(self, field):
		
		if field=='xn':
			value = np.zeros((self.nx, 1))
			tmp = cast(value.ctypes.data, POINTER(c_double))
			self.__acados.sim_out_get_xn(self.config, self.dims, self.sim_out, tmp)

		if field=='Sxn':
			value = np.zeros((self.nx, self.nx), order='F')
			tmp = cast(value.ctypes.data, POINTER(c_double))
			self.__acados.sim_out_get_Sxn(self.config, self.dims, self.sim_out, tmp)

		return value



#	def unload_model(self):
#		_ctypes.dlclose(self.__model._handle)



	def __del__(self):
		
#		print('\nin descruction\n')
		if self.scheme=='erk':
			if self.sens_forw=='false':
				self.__acados.external_function_casadi_free(self.expl_ode_fun)
			else:
				self.__acados.external_function_casadi_free(self.expl_vde_for)
		if self.scheme=='irk':
			if self.sens_forw=='false':
				self.__acados.external_function_casadi_free(self.impl_ode_fun)
				self.__acados.external_function_casadi_free(self.impl_ode_fun_jac_x_xdot_z)
		self.__acados.sim_config_destroy(self.config)
		self.__acados.sim_dims_destroy(self.dims)
		self.__acados.sim_opts_destroy(self.opts)
		self.__acados.sim_out_destroy(self.sim_in)
		self.__acados.sim_in_destroy(self.sim_out)
		self.__acados.sim_solver_destroy(self.solver)
#		print(self.__model)
#		print(self.__model._handle)
		# on POSIX systems; on windows call FreeLibrary instead
#		print(self.__model._handle)
		_ctypes.dlclose(self.__model._handle)



