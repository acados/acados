def generate_wrapper(name):
	print(name)

	file = open("wrapper_"+name+".c", "w")

	file.write("int "+name+"(const double** arg, double** res, int* iw, double* w, void* mem);\n")
	file.write("int "+name+"_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);\n")
	file.write("const int* "+name+"_sparsity_in(int i);\n")
	file.write("const int* "+name+"_sparsity_out(int i);\n")
	file.write("int "+name+"_n_in(void);\n")
	file.write("int "+name+"_n_out(void);\n")
	file.write("\n")
	file.write("void *get_"+name+"_fun() { return &"+name+"; }\n")
	file.write("void *get_"+name+"_work() { return &"+name+"_work; }\n")
	file.write("void *get_"+name+"_sparsity_in() { return &"+name+"_sparsity_in; }\n")
	file.write("void *get_"+name+"_sparsity_out() { return &"+name+"_sparsity_out; }\n")
	file.write("void *get_"+name+"_n_in() { return &"+name+"_n_in; }\n")
	file.write("void *get_"+name+"_n_out() { return &"+name+"_n_out; }\n")


	file.close()


