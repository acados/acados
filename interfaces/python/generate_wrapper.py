def generate_wrapper(user_name, acados_name):

#	print(user_name)
#	print(acados_name)

	file = open("wrapper_"+user_name+".c", "w")

	file.write("int "+user_name+"(const double** arg, double** res, int* iw, double* w, void* mem);\n")
	file.write("int "+user_name+"_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);\n")
	file.write("const int* "+user_name+"_sparsity_in(int i);\n")
	file.write("const int* "+user_name+"_sparsity_out(int i);\n")
	file.write("int "+user_name+"_n_in(void);\n")
	file.write("int "+user_name+"_n_out(void);\n")
	file.write("\n")
	file.write("void *get_"+acados_name+"_fun() { return &"+user_name+"; }\n")
	file.write("void *get_"+acados_name+"_work() { return &"+user_name+"_work; }\n")
	file.write("void *get_"+acados_name+"_sparsity_in() { return &"+user_name+"_sparsity_in; }\n")
	file.write("void *get_"+acados_name+"_sparsity_out() { return &"+user_name+"_sparsity_out; }\n")
	file.write("void *get_"+acados_name+"_n_in() { return &"+user_name+"_n_in; }\n")
	file.write("void *get_"+acados_name+"_n_out() { return &"+user_name+"_n_out; }\n")


	file.close()


