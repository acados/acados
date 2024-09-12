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



% NOTE: `acados` currently supports both an old MATLAB/Octave interface (< v0.4.0)
% as well as a new interface (>= v0.4.0).

% THIS EXAMPLE still uses the OLD interface. If you are new to `acados` please start
% with the examples that have been ported to the new interface already.
% see https://github.com/acados/acados/issues/1196#issuecomment-2311822122)

clear all; clc;

% check that env.sh has been run
env_run = getenv('ENV_RUN');
if (~strcmp(env_run, 'true'))
	error('env.sh has not been sourced! Before executing this example, run: source env.sh');
end

%% options
compile_interface = 'auto'; % true, false
gnsf_detect_struct = 'true'; % true, false
method = 'irk'; % irk, irk_gnsf, [erk]
sens_forw = 'true'; % true, false
jac_reuse = 'false'; % true, false
num_stages = 3;
num_steps = 3;
newton_iter = 3;
model_name = 'pend_dae';

% x = [xpos, ypos, alpha, vx, vy, valpha]
% x0 = [1; -5; 1; 0.1; -0.5; 0.1];
% x0 = zeros(nx, 1); %
length_pendulum = 5;

alpha0 = 0.1;
xp0 = length_pendulum * sin(alpha0);
yp0 = - length_pendulum * cos(alpha0);
x0 = [ xp0; yp0; alpha0; 0; 0; 0];

u = 0;

% steady state
% x0 = [ 0; -length_pendulum; 0; 0; 0; 0];

%% model
model = pendulum_dae_model();
disp('state')
disp(model.sym_x)

nx = length(model.sym_x);
nu = length(model.sym_u);
nz = length(model.sym_z);

%% acados sim model
sim_model = acados_sim_model();
sim_model.set('name', model_name);
sim_model.set('T', 0.1); % simulation time

sim_model.set('sym_x', model.sym_x);
if isfield(model, 'sym_u')
    sim_model.set('sym_u', model.sym_u);
end
if isfield(model, 'sym_p')
    sim_model.set('sym_p', model.sym_p);
end

% Note: DAEs can only be used with implicit integrator
sim_model.set('dyn_type', 'implicit');
sim_model.set('dyn_expr_f', model.expr_f_impl);
sim_model.set('sym_xdot', model.sym_xdot);
if isfield(model, 'sym_z')
	sim_model.set('sym_z', model.sym_z);
end

%% acados sim opts
sim_opts = acados_sim_opts();
sim_opts.set('compile_interface', compile_interface);
sim_opts.set('num_stages', num_stages);
sim_opts.set('num_steps', num_steps);
sim_opts.set('newton_iter', newton_iter);
sim_opts.set('method', method);
sim_opts.set('sens_forw', sens_forw);
sim_opts.set('sens_adj', 'true');
sim_opts.set('sens_algebraic', 'true');
sim_opts.set('output_z', 'true');
sim_opts.set('sens_hess', 'false');
sim_opts.set('jac_reuse', jac_reuse);
if (strcmp(method, 'irk_gnsf'))
	sim_opts.set('gnsf_detect_struct', gnsf_detect_struct);
end


%% acados sim
% create integrator
sim_solver = acados_sim(sim_model, sim_opts);

N_sim = 100;

x_sim = zeros(nx, N_sim+1);
x_sim(:,1) = x0;

% initialization
xdot0 = zeros(nx, 1);
z0 = zeros(nz, 1);

tic
for ii=1:N_sim

	% set initial state
	sim_solver.set('x', x_sim(:,ii));
	sim_solver.set('u', u);

	% set adjoint seed
	sim_solver.set('seed_adj', ones(nx,1));

    % initialize implicit integrator
    if (strcmp(method, 'irk'))
        sim_solver.set('xdot', xdot0);
        sim_solver.set('z', z0);
    elseif (strcmp(method, 'irk_gnsf'))
        y_in = sim_solver.sim.model.dyn_gnsf_L_x * x0 ...
                + sim_solver.sim.model.dyn_gnsf_L_xdot * xdot0 ...
                + sim_solver.sim.model.dyn_gnsf_L_z * z0;
        u_hat = sim_solver.sim.model.dyn_gnsf_L_u * u;
        phi_fun = Function([model_name,'_gnsf_phi_fun'],...
                        {sim_solver.sim.model.sym_gnsf_y, sim_solver.sim.model.sym_gnsf_uhat},...
                            {sim_solver.sim.model.dyn_gnsf_expr_phi(:)});

        phi_guess = full( phi_fun( y_in, u_hat ) );
        n_out = sim_solver.sim.dims.gnsf_nout;
        sim_solver.set('phi_guess', zeros(n_out,1));
    end

	% solve
	sim_solver.solve();

	% get simulated state
	x_sim(:,ii+1) = sim_solver.get('xn');

	% forward sensitivities
	S_forw = sim_solver.get('S_forw');
	Sx = sim_solver.get('Sx');
	Su = sim_solver.get('Su');

end
format short e
xfinal = sim_solver.get('xn')';
S_adj = sim_solver.get('S_adj')';
z = sim_solver.get('zn')'; % approximate value of algebraic variables at start of simulation
S_alg = sim_solver.get('S_algebraic'); % sensitivities of algebraic variables z

simulation_time = toc


figure;
Nsub = 4;
subplot(Nsub, 1, 1);
plot(1:N_sim+1, x_sim(1,:));
legend('xpos');

subplot(Nsub, 1, 2);
plot(1:N_sim+1, x_sim(2,:));
legend('ypos');

subplot(Nsub, 1, 3);
plot(1:N_sim+1, x_sim(3,:));
legend('alpha')

subplot(Nsub, 1, 4);
plot(1:N_sim+1, x_sim(4:6,:));
legend('vx', 'vy', 'valpha');

if is_octave()
    waitforbuttonpress;
end

% check consistency
xp = x_sim(1,:);
yp = x_sim(2,:);
check = xp.^2 + yp.^2 - length_pendulum^2;
if any( max(abs(check)) > 1e-10 )
    disp('note: check for constant pendulum length failed');
end

