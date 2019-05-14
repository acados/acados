function compile_casadi_functions(model_struct, opts_struct)

c_sources = ' ';
if (strcmp(opts_struct.method, 'erk'))
	% generate c for function and derivatives using casadi
	generate_c_code_explicit_ode(model_struct, opts_struct);
	% compile the code in a shared library
	c_sources = [c_sources, 'sim_model_dyn_expl_ode_fun.c '];
	c_sources = [c_sources, 'sim_model_dyn_expl_vde_for.c '];
	c_sources = [c_sources, 'sim_model_dyn_expl_vde_adj.c '];
	c_sources = [c_sources, 'sim_model_dyn_expl_ode_hes.c '];
elseif (strcmp(opts_struct.method, 'irk'))
	% generate c for function and derivatives using casadi
	generate_c_code_implicit_ode(model_struct, opts_struct);
	% compile the code in a shared library
	c_sources = [c_sources, 'sim_model_dyn_impl_ode_fun.c '];
	c_sources = [c_sources, 'sim_model_dyn_impl_ode_fun_jac_x_xdot.c '];
	c_sources = [c_sources, 'sim_model_dyn_impl_ode_fun_jac_x_xdot_u.c '];
	c_sources = [c_sources, 'sim_model_dyn_impl_ode_jac_x_xdot_u.c '];
	c_sources = [c_sources, 'sim_model_dyn_impl_ode_hess.c '];
else
	fprintf('\ncodegen_model: method not supported: %s\n', opts_struct.method);
	return;
end
lib_name = ['libsim_model.so'];
system(['gcc -O2 -fPIC -shared ', c_sources, ' -o ', lib_name]);

