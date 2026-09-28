%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.


function model = get_linear_mass_spring_model()

	import casadi.*

	%% dims
	num_mass = 4;

	nx = 2*num_mass;
	nu = num_mass-1;

	%% symbolic variables
	sym_x = SX.sym('x', nx, 1); % states
	sym_u = SX.sym('u', nu, 1); % controls
	sym_xdot = SX.sym('xdot',size(sym_x)); %state derivatives

	%% dynamics
	% continuous time
	Ac = zeros(nx, nx);
	for ii=1:num_mass
		Ac(ii,num_mass+ii) = 1.0;
		Ac(num_mass+ii,ii) = -2.0;
	end
	for ii=1:num_mass-1
		Ac(num_mass+ii,ii+1) = 1.0;
		Ac(num_mass+ii+1,ii) = 1.0;
	end

	Bc = zeros(nx, nu);
	for ii=1:nu
		Bc(num_mass+ii, ii) = 1.0;
	end

	c_const = zeros(nx, 1);

	% discrete time
	Ts = 0.5; % sampling time
	M = expm([Ts*Ac, Ts*Bc; zeros(nu, 2*nx/2+nu)]);
	A = M(1:nx,1:nx);
	B = M(1:nx,nx+1:end);

	dyn_expr_f_expl = Ac*sym_x + Bc*sym_u + c_const;
	dyn_expr_f_impl = dyn_expr_f_expl - sym_xdot;
	dyn_expr_phi = A*sym_x + B*sym_u;

	%% populate structure
	model = AcadosModel();
	model.x = sym_x;
	model.xdot = sym_xdot;
	model.u = sym_u;
	model.f_expl_expr = dyn_expr_f_expl;
	model.f_impl_expr = dyn_expr_f_impl;
	model.disc_dyn_expr = dyn_expr_phi;
end
