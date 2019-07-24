function sim_generate_casadi_ext_fun(model_struct, opts_struct)

model_name = model_struct.name;

c_sources = ' ';
if (strcmp(opts_struct.method, 'erk'))
	% generate c for function and derivatives using casadi
	generate_c_code_explicit_ode(model_struct, opts_struct);
	% compile the code in a shared library
	c_sources = [c_sources, model_name, '_dyn_expl_ode_fun.c '];
	c_sources = [c_sources, model_name, '_dyn_expl_vde_for.c '];
	c_sources = [c_sources, model_name, '_dyn_expl_vde_adj.c '];
	c_sources = [c_sources, model_name, '_dyn_expl_ode_hes.c '];
elseif (strcmp(opts_struct.method, 'irk'))
	% generate c for function and derivatives using casadi
	generate_c_code_implicit_ode(model_struct, opts_struct);
	% compile the code in a shared library
	c_sources = [c_sources, model_name, '_dyn_impl_ode_fun.c '];
	c_sources = [c_sources, model_name, '_dyn_impl_ode_fun_jac_x_xdot.c '];
	c_sources = [c_sources, model_name, '_dyn_impl_ode_fun_jac_x_xdot_u.c '];
	c_sources = [c_sources, model_name, '_dyn_impl_ode_jac_x_xdot_u.c '];
	c_sources = [c_sources, model_name, '_dyn_impl_ode_hess.c '];
elseif (strcmp(opts_struct.method, 'irk_gnsf'))
	% generate c for function and derivatives using casadi
	generate_c_code_gnsf(model_struct); %, opts_struct);
	% compile the code in a shared library
	c_sources = [c_sources, model_name, '_dyn_gnsf_f_lo_fun_jac_x1k1uz.c '];
	c_sources = [c_sources, model_name, '_dyn_gnsf_get_matrices_fun.c '];
	c_sources = [c_sources, model_name, '_dyn_gnsf_phi_fun.c '];
	c_sources = [c_sources, model_name, '_dyn_gnsf_phi_fun_jac_y.c '];
	c_sources = [c_sources, model_name, '_dyn_gnsf_phi_jac_y_uhat.c '];
else
	fprintf('\nsim_generate_casadi_ext_fun: method not supported: %s\n', opts_struct.method);
	return;
end

if ispc
  lib_name = ['lib', model_name, '.lib'];
else
  lib_name = ['lib', model_name, '.so'];
end

system(['gcc -O2 -fPIC -shared ', c_sources, ' -o ', lib_name]);

c_files = split(c_sources);
for k=1:length(c_files)
  if ~isempty(c_files{k})
    movefile(c_files{k}, 'build')
  end
end

movefile(lib_name, 'build')
