function ocp_generate_casadi_ext_fun(model_struct, opts_struct)

model_name = model_struct.name;

% select files to compile
c_sources = ' ';
% dynamics
if (strcmp(model_struct.dyn_type, 'explicit'))
	% generate c for function and derivatives using casadi
	generate_c_code_explicit_ode(model_struct, opts_struct);
	% sources list
	c_sources = [c_sources, model_name, '_dyn_expl_ode_fun.c '];
	c_sources = [c_sources, model_name, '_dyn_expl_vde_for.c '];
	c_sources = [c_sources, model_name, '_dyn_expl_vde_adj.c '];
	c_sources = [c_sources, model_name, '_dyn_expl_ode_hes.c '];
elseif (strcmp(model_struct.dyn_type, 'implicit'))
	if (strcmp(opts_struct.sim_method, 'irk'))
		% generate c for function and derivatives using casadi
		generate_c_code_implicit_ode(model_struct, opts_struct);
		% sources list
		c_sources = [c_sources, model_name, '_dyn_impl_ode_fun.c '];
		c_sources = [c_sources, model_name, '_dyn_impl_ode_fun_jac_x_xdot.c '];
		c_sources = [c_sources, model_name, '_dyn_impl_ode_fun_jac_x_xdot_u.c '];
		c_sources = [c_sources, model_name, '_dyn_impl_ode_jac_x_xdot_u.c '];
		c_sources = [c_sources, model_name, '_dyn_impl_ode_hess.c '];
	elseif (strcmp(opts_struct.sim_method, 'irk_gnsf'))
		% generate c for function and derivatives using casadi
		generate_c_code_gnsf(model_struct); %, opts_struct);
		% compile the code in a shared library
		c_sources = [c_sources, model_name, '_dyn_gnsf_f_lo_fun_jac_x1k1uz.c '];
		c_sources = [c_sources, model_name, '_dyn_gnsf_get_matrices_fun.c '];
		c_sources = [c_sources, model_name, '_dyn_gnsf_phi_fun.c '];
		c_sources = [c_sources, model_name, '_dyn_gnsf_phi_fun_jac_y.c '];
		c_sources = [c_sources, model_name, '_dyn_gnsf_phi_jac_y_uhat.c '];
	else
		fprintf('\nocp_generate_casadi_ext_fun: sim_method not supported: %s\n', opts_struct.sim_method);
		return;
	end
elseif (strcmp(model_struct.dyn_type, 'discrete'))
	% generate c for function and derivatives using casadi
	generate_c_code_disc_dyn(model_struct, opts_struct);
	% sources list
	c_sources = [c_sources, model_name, '_dyn_disc_phi_fun_jac.c '];
	c_sources = [c_sources, model_name, '_dyn_disc_phi_fun_jac_hess.c '];
else
	fprintf('\ncodegen_model: dyn_type not supported: %s\n', model_struct.dyn_type);
	return;
end
% nonlinear constraints
if (strcmp(model_struct.constr_type, 'bgh') && (isfield(model_struct, 'constr_expr_h') || isfield(model_struct, 'constr_expr_h_e')))
	% generate c for function and derivatives using casadi
	generate_c_code_nonlinear_constr(model_struct, opts_struct);
	% sources list
	if isfield(model_struct, 'constr_expr_h')
		c_sources = [c_sources, model_name, '_constr_h_fun_jac_ut_xt.c '];
		c_sources = [c_sources, model_name, '_constr_h_fun_jac_ut_xt_hess.c '];
	end
	if isfield(model_struct, 'constr_expr_h_e')
		c_sources = [c_sources, model_name, '_constr_h_e_fun_jac_ut_xt.c '];
		c_sources = [c_sources, model_name, '_constr_h_e_fun_jac_ut_xt_hess.c '];
	end
end
% nonlinear least squares
if (strcmp(model_struct.cost_type, 'nonlinear_ls') || strcmp(model_struct.cost_type_e, 'nonlinear_ls'))
	% generate c for function and derivatives using casadi
	generate_c_code_nonlinear_least_squares(model_struct, opts_struct);
	% sources list
	if isfield(model_struct, 'cost_expr_y')
		c_sources = [c_sources, model_name, '_cost_y_fun_jac_ut_xt.c '];
		c_sources = [c_sources, model_name, '_cost_y_hess.c '];
	end
	if isfield(model_struct, 'cost_expr_y_e')
		c_sources = [c_sources, model_name, '_cost_y_e_fun_jac_ut_xt.c '];
		c_sources = [c_sources, model_name, '_cost_y_e_hess.c '];
	end
end
% external cost
if (strcmp(model_struct.cost_type, 'ext_cost') || strcmp(model_struct.cost_type_e, 'ext_cost'))
	% generate c for function and derivatives using casadi
	generate_c_code_ext_cost(model_struct, opts_struct);
	% sources list
	if isfield(model_struct, 'cost_expr_ext_cost')
		c_sources = [c_sources, model_name, '_cost_ext_cost_jac_hes.c '];
	end
	if isfield(model_struct, 'cost_expr_ext_cost_e')
		c_sources = [c_sources, model_name, '_cost_ext_cost_e_jac_hes.c '];
	end
end

c_files = strsplit(c_sources);
c_files = c_files(~cellfun(@isempty, c_files)); % remove empty cells

if ispc
  ldext = '.lib';
else
  ldext = '.so';
end

lib_name = 'ocp_model';

if ispc
  mbuild(c_files{:}, '-output', lib_name, 'CFLAGS="-fPIC $CFLAGS"', 'LDTYPE="-shared"', ['LDEXT=', ldext]);
else
  system(['gcc -O2 -fPIC -shared ', c_sources, ' -o ', [lib_name, ldext]]);
end

movefile([lib_name, ldext], fullfile(opts_struct.output_dir, [lib_name, ldext]));

for k=1:length(c_files)
  movefile(c_files{k}, opts_struct.output_dir);
end


