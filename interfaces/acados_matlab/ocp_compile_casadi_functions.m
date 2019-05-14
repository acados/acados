function ocp_compile_casadi_functions(model_struct, opts_struct)

% select files to compile
c_sources = ' ';
% dynamics
if (strcmp(model_struct.dyn_type, 'explicit'))
	% generate c for function and derivatives using casadi
	generate_c_code_explicit_ode(model_struct, opts_struct);
	% sources list
	c_sources = [c_sources, 'ocp_model_dyn_expl_ode_fun.c '];
	c_sources = [c_sources, 'ocp_model_dyn_expl_vde_for.c '];
	c_sources = [c_sources, 'ocp_model_dyn_expl_vde_adj.c '];
	c_sources = [c_sources, 'ocp_model_dyn_expl_ode_hes.c '];
elseif (strcmp(model_struct.dyn_type, 'implicit'))
	% generate c for function and derivatives using casadi
	generate_c_code_implicit_ode(model_struct, opts_struct);
	% sources list
	c_sources = [c_sources, 'ocp_model_dyn_impl_ode_fun.c '];
	c_sources = [c_sources, 'ocp_model_dyn_impl_ode_fun_jac_x_xdot.c '];
	c_sources = [c_sources, 'ocp_model_dyn_impl_ode_fun_jac_x_xdot_u.c '];
	c_sources = [c_sources, 'ocp_model_dyn_impl_ode_jac_x_xdot_u.c '];
	c_sources = [c_sources, 'ocp_model_dyn_impl_ode_hess.c '];
elseif (strcmp(model_struct.dyn_type, 'discrete'))
	% generate c for function and derivatives using casadi
	generate_c_code_disc_dyn(model_struct, opts_struct);
	% sources list
	c_sources = [c_sources, 'ocp_model_dyn_disc_phi_fun_jac.c '];
	c_sources = [c_sources, 'ocp_model_dyn_disc_phi_fun_jac_hess.c '];
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
		c_sources = [c_sources, 'ocp_model_constr_h_fun_jac_ut_xt.c '];
		c_sources = [c_sources, 'ocp_model_constr_h_fun_jac_ut_xt_hess.c '];
	end
	if isfield(model_struct, 'constr_expr_h_e')
		c_sources = [c_sources, 'ocp_model_constr_h_e_fun_jac_ut_xt.c '];
		c_sources = [c_sources, 'ocp_model_constr_h_e_fun_jac_ut_xt_hess.c '];
	end
end
% nonlinear least squares
if (strcmp(model_struct.cost_type, 'nonlinear_ls') || strcmp(model_struct.cost_type_e, 'nonlinear_ls'))
	% generate c for function and derivatives using casadi
	generate_c_code_nonlinear_least_squares(model_struct, opts_struct);
	% sources list
	if isfield(model_struct, 'cost_expr_y')
		c_sources = [c_sources, 'ocp_model_cost_y_fun_jac_ut_xt.c '];
		c_sources = [c_sources, 'ocp_model_cost_y_hess.c '];
	end
	if isfield(model_struct, 'cost_expr_y_e')
		c_sources = [c_sources, 'ocp_model_cost_y_e_fun_jac_ut_xt.c '];
		c_sources = [c_sources, 'ocp_model_cost_y_e_hess.c '];
	end
end
% external cost
if (strcmp(model_struct.cost_type, 'ext_cost') || strcmp(model_struct.cost_type_e, 'ext_cost'))
	% generate c for function and derivatives using casadi
	generate_c_code_ext_cost(model_struct, opts_struct);
	% sources list
	if isfield(model_struct, 'cost_expr_ext_cost')
		c_sources = [c_sources, 'ocp_model_cost_ext_cost_jac_hes.c '];
	end
	if isfield(model_struct, 'cost_expr_ext_cost_e')
		c_sources = [c_sources, 'ocp_model_cost_ext_cost_e_jac_hes.c '];
	end
end
lib_name = ['libocp_model.so'];
system(['gcc -O2 -fPIC -shared ', c_sources, ' -o ', lib_name]);

