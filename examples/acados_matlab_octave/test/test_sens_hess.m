%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% The 2-Clause BSD License
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.;

%

%% test of native matlab interface

addpath('../pendulum_on_cart_model/');

% TODO: include irk_gnsf, as soon as hessians are implemented
for integrator = {'erk', 'irk'} %, 'irk_gnsf'}
	%% arguments
	compile_interface = 'auto';
	codgen_model = 'true';
	method = integrator{1}; %'irk'; 'irk_gnsf'; 'erk';
	sens_forw = 'true';
	sens_adj = 'true';
	sens_hess = 'true';
	num_stages = 4;
	num_steps = 3;

	Ts = 0.1;
	x0 = [1e-1; 1e0; 2e-1; 2e0];
	u = 0;
	FD_epsilon = 1e-6;

	%% model
	model = pendulum_on_cart_model();

	model_name = ['pendulum_' method];
	nx = model.nx;
	nu = model.nu;

	%% acados sim model
	sim_model = acados_sim_model();
	sim_model.set('T', Ts);
	sim_model.set('name', model_name);

    sim_model.set('sym_x', model.sym_x);
    if isfield(model, 'sym_u')
        sim_model.set('sym_u', model.sym_u);
    end
    if isfield(model, 'sym_p')
        sim_model.set('sym_p', model.sym_p);
    end

    if (strcmp(method, 'erk'))
        sim_model.set('dyn_type', 'explicit');
        sim_model.set('dyn_expr_f', model.expr_f_expl);
    else % irk irk_gnsf
        sim_model.set('dyn_type', 'implicit');
        sim_model.set('dyn_expr_f', model.expr_f_impl);
        sim_model.set('sym_xdot', model.sym_xdot);
    %	if isfield(model, 'sym_z')
    %		sim_model.set('sym_z', model.sym_z);
    %	end
    end


	%% acados sim opts
	sim_opts = acados_sim_opts();
	sim_opts.set('compile_interface', compile_interface);
	sim_opts.set('codgen_model', codgen_model);
	sim_opts.set('num_stages', num_stages);
	sim_opts.set('num_steps', num_steps);
	sim_opts.set('method', method);
	sim_opts.set('sens_forw', sens_forw);
	sim_opts.set('sens_adj', sens_adj);
	sim_opts.set('sens_hess', sens_hess);
	%sim_opts.opts_struct

	%% acados sim
	% create sim
	sim = acados_sim(sim_model, sim_opts);

	% compute hessian sensitivities using internal numerical differentiation
	S_hess_ind = zeros(nx+nu, nx+nu, nx);

	% compute hessian sensitivities using finite differences
	S_hess_fd = zeros(nx+nu, nx+nu, nx);

	for jj=1:nx % loop over unit seeds
		% set initial state
		sim.set('x', x0);
		sim.set('u', u);

		% internal numerical differentiation seed
		lambda = zeros(nx, 1);
		lambda(jj) = 1.0;
		sim.set('seed_adj', lambda);

		% solve
		sim.solve();

		% S_hess
		S_hess = sim.get('S_hess');
		S_hess_ind(:, :, jj) = S_hess;

		% S_adj
		S_adj = sim.get('S_adj');

		%% asymmetric finite differences
		for ii=1:nx
			dx = zeros(nx, 1);
			dx(ii) = 1.0;

			sim.set('x', x0+FD_epsilon*dx);
			sim.set('u', u);

			sim.solve();
			S_adj_tmp = sim.get('S_adj');
			S_hess_fd(:, ii, jj) = (S_adj_tmp - S_adj) / FD_epsilon;
		
		end

		for ii=1:nu
			du = zeros(nu, 1);
			du(ii) = 1.0;

			sim.set('x', x0);
			sim.set('u', u+FD_epsilon*du);

			sim.solve();
			S_adj_tmp = sim.get('S_adj');
			S_hess_fd(:, nx+ii, jj) = (S_adj_tmp - S_adj) / FD_epsilon;
        end
	end

	%% compute & check error
	error_abs = max(max(max(abs(S_hess_fd - S_hess_ind))));
	disp(' ')
	disp(['integrator:  ' method]);
	disp(['error hessian (wrt finite differences):   ' num2str(error_abs)])
	disp(' ')
    if error_abs > 1e-6
        disp(['hessian error too large'])
        error(strcat('test_sens_hess FAIL: second order sensitivity error too large: \n',...
			'for integrator:\t', method));
	end
end

fprintf('\nTEST_HESSIANS: success!\n\n');

return;
