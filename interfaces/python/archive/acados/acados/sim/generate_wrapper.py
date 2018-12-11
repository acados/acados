from ctypes import *
import ctypes.util 
import _ctypes
from os import system



def generate_getter(user_fun_name):

	file = open("casadi_fun_ptr_getter.c", "w")

	file.write("int "+user_fun_name+"(const double** arg, double** res, int* iw, double* w, void* mem);\n")
	file.write("int "+user_fun_name+"_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);\n")
	file.write("const int* "+user_fun_name+"_sparsity_in(int i);\n")
	file.write("const int* "+user_fun_name+"_sparsity_out(int i);\n")
	file.write("int "+user_fun_name+"_n_in(void);\n")
	file.write("int "+user_fun_name+"_n_out(void);\n")
	file.write("\n")
#	file.write("void *get_fun_fun() { return &"+user_fun_name+"; }\n")
	file.write("void *get_fun_fun() { return &"+user_fun_name+"; }\n")
	file.write("void *get_fun_work() { return &"+user_fun_name+"_work; }\n")
	file.write("void *get_fun_sparsity_in() { return &"+user_fun_name+"_sparsity_in; }\n")
	file.write("void *get_fun_sparsity_out() { return &"+user_fun_name+"_sparsity_out; }\n")
	file.write("void *get_fun_n_in() { return &"+user_fun_name+"_n_in; }\n")
	file.write("void *get_fun_n_out() { return &"+user_fun_name+"_n_out; }\n")


	file.close()



def set_function_pointers(acados, model_name, user_fun_name, acados_ext_fun):
	
#	print(user_fun_name)

	# generate function pointers getter
	generate_getter(user_fun_name)

	# compile and load shared library
	model_getter_name = user_fun_name + '_getter.so'
#	model_getter_name = 'fun_getter.so' # XXX it needs an unique name !!!
	system('gcc -fPIC -shared casadi_fun_ptr_getter.c -o ' + model_getter_name + ' -L. ' + model_name)
#	print(model_name)
#	system('nm ' + model_name)
	model_getter = CDLL(model_getter_name)

	# set up function pointers
	model_getter.get_fun_fun.restype = c_void_p
	tmp_ptr = model_getter.get_fun_fun()
#	print(tmp_ptr)
	acados.external_function_casadi_set_fun.argtypes = [c_void_p, c_void_p]
	acados.external_function_casadi_set_fun(acados_ext_fun, tmp_ptr)

	model_getter.get_fun_work.restype = c_void_p
	tmp_ptr = model_getter.get_fun_work()
	acados.external_function_casadi_set_work.argtypes = [c_void_p, c_void_p]
	acados.external_function_casadi_set_work(acados_ext_fun, tmp_ptr)

	model_getter.get_fun_sparsity_in.restype = c_void_p
	tmp_ptr = model_getter.get_fun_sparsity_in()
	acados.external_function_casadi_set_sparsity_in.argtypes = [c_void_p, c_void_p]
	acados.external_function_casadi_set_sparsity_in(acados_ext_fun, tmp_ptr)

	model_getter.get_fun_sparsity_out.restype = c_void_p
	tmp_ptr = model_getter.get_fun_sparsity_out()
	acados.external_function_casadi_set_sparsity_out.argtypes = [c_void_p, c_void_p]
	acados.external_function_casadi_set_sparsity_out(acados_ext_fun, tmp_ptr)

	model_getter.get_fun_n_in.restype = c_void_p
	tmp_ptr = model_getter.get_fun_n_in()
	acados.external_function_casadi_set_n_in.argtypes = [c_void_p, c_void_p]
	acados.external_function_casadi_set_n_in(acados_ext_fun, tmp_ptr)

	model_getter.get_fun_n_out.restype = c_void_p
	tmp_ptr = model_getter.get_fun_n_out()
	acados.external_function_casadi_set_n_out.argtypes = [c_void_p, c_void_p]
	acados.external_function_casadi_set_n_out(acados_ext_fun, tmp_ptr)

	_ctypes.dlclose(model_getter._handle)
	system('rm casadi_fun_ptr_getter.c')
	system('rm ' + model_getter_name)




