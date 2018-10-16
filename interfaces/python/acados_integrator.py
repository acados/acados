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
		self.z = SX.sym('z', 0, 1)
	


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

		elif field=='z':
			self.z = value

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


		self.scheme = opts.scheme
		self.sens_forw = opts.sens_forw

		
#		print(CasadiMeta.version())

		# load acados library
		__acados = CDLL('libacados_c.so')
		self.__acados = __acados


		# sizes
		self.nx = model.x.size()[0]
		self.nu = model.u.size()[0]
		self.nz = model.z.size()[0]

		nx = self.nx
		nu = self.nu
		nz = self.nz


		# define functions & generate C code
		casadi_opts = dict(casadi_int='int', casadi_real='double')
		c_sources = ' '

		if opts.scheme=='erk':

			if opts.sens_forw=='false':

				fun_name = 'expl_ode_fun'
				fun = Function(fun_name, [model.x, model.u], [model.ode_expr])
				fun.generate(casadi_opts)
				c_sources = c_sources + fun_name + '.c '

			else:

				fun_name = 'expl_vde_for'
				Sx = SX.sym('Sx', nx, nx)
				Su = SX.sym('Su', nx, nu)
				vde_x = SX.zeros(nx, nx) + jtimes(model.ode_expr, model.x, Sx) # TODO try sparse !!!
				vde_u = SX.zeros(nx, nu) + jacobian(model.ode_expr, model.u) + jtimes(model.ode_expr, model.x, Su) # TODO try sparse !!!
				fun = Function(fun_name, [model.x, Sx, Su, model.u], [model.ode_expr, vde_x, vde_u])
				fun.generate(casadi_opts)
				c_sources = c_sources + fun_name + '.c '

		elif opts.scheme=='irk':

			fun_name = 'impl_ode_fun'
			fun = Function(fun_name, [model.x, model.xdot, model.u, model.z], [model.ode_expr])
			fun.generate(casadi_opts)
			c_sources = c_sources + fun_name + '.c '

			fun_name = 'impl_ode_fun_jac_x_xdot_z'
			jac_x = SX.zeros(nx, nx) + jacobian(model.ode_expr, model.x) # TODO try sparse !!!
			jac_xdot = SX.zeros(nx, nx) + jacobian(model.ode_expr, model.xdot) # TODO try sparse !!!
			jac_z = SX.zeros(nx, nz) + jacobian(model.ode_expr, model.z) # TODO try sparse !!!
			fun = Function(fun_name, [model.x, model.xdot, model.u, model.z], [model.ode_expr, jac_x, jac_xdot, jac_z])
			fun.generate(casadi_opts)
			c_sources = c_sources + fun_name + '.c '

			if opts.sens_forw=='true':

				fun_name = 'impl_ode_jac_x_xdot_u_z'
				jac_x = SX.zeros(nx, nx) + jacobian(model.ode_expr, model.x) # TODO try sparse !!!
				jac_xdot = SX.zeros(nx, nx) + jacobian(model.ode_expr, model.xdot) # TODO try sparse !!!
				jac_u = SX.zeros(nx, nu) + jacobian(model.ode_expr, model.u) # TODO try sparse !!!
				jac_z = SX.zeros(nx, nz) + jacobian(model.ode_expr, model.z) # TODO try sparse !!!
				fun = Function(fun_name, [model.x, model.xdot, model.u, model.z], [jac_x, jac_xdot, jac_u, jac_z])
				fun.generate(casadi_opts)
				c_sources = c_sources + fun_name + '.c '


		# create model library
		lib_name = model.model_name + '.so'
		system('gcc -fPIC -shared ' + c_sources + ' -o ' + lib_name)

		## load model library
		self.__model = CDLL(lib_name)

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
		if opts.scheme=='erk':
			self.config = cast(__acados.sim_config_create( 0 ), c_void_p)
		elif opts.scheme=='irk':
			self.config = cast(__acados.sim_config_create( 2 ), c_void_p)

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
			Sx0 = np.zeros((nx, nx), order='F')
			for ii in range(nx):
				Sx0[ii][ii] = 1.0
			tmp = np.ascontiguousarray(Sx0, dtype=np.float64)
			tmp = cast(tmp.ctypes.data, POINTER(c_double))
			self.__acados.sim_in_set_Sx(self.config, self.dims, tmp, self.sim_in)
			# set Su
			Su0 = np.zeros((nx, nu), order='F')
			tmp = np.ascontiguousarray(Su0, dtype=np.float64)
			tmp = cast(tmp.ctypes.data, POINTER(c_double))
			self.__acados.sim_in_set_Su(self.config, self.dims, tmp, self.sim_in)

		## sim_out
		self.sim_out = cast(__acados.sim_out_create(self.config, self.dims), c_void_p)

		## sim solver
		self.solver = cast(__acados.sim_create(self.config, self.dims, self.opts), c_void_p)



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



	def __del__(self):
		
		if self.scheme=='erk':
			if self.sens_forw=='false':
				self.__acados.external_function_casadi_free(self.expl_ode_fun)
			else:
				self.__acados.external_function_casadi_free(self.expl_vde_for)
		if self.scheme=='irk':
			if self.sens_forw=='false':
				self.__acados.external_function_casadi_free(self.impl_ode_fun)
				self.__acados.external_function_casadi_free(self.impl_ode_fun_jac_x_xdot_z)
		self.__acados.sim_config_free(self.config)
		self.__acados.sim_dims_free(self.dims)
		self.__acados.sim_opts_free(self.opts)
		self.__acados.sim_out_free(self.sim_in)
		self.__acados.sim_in_free(self.sim_out)
		self.__acados.sim_free(self.solver)



