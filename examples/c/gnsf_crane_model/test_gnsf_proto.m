clear all
% close all

import acados.*
import casadi.*

export_crane_model_gnsf;
load crane_controls.mat
controls_MPC = controls_MPC';
X0 = [0 0 0.8 0 0 0 0 0 0]';

stepsize = 0.1;
time_grid = 0:stepsize:stepsize; % dt; %5;
n_simulations = 10;

%% set simulation options

q = 4; % number of stages;
n_steps = 2;
disp(['Number of stages for collocation = ' num2str(q)]);
disp(['Number of steps per simulation =   ' num2str(n_steps)]);
adj_sens_mode = 1;
forw_sens_mode = 1;
opts = set_integrator_opts(q, stepsize, n_steps, forw_sens_mode, adj_sens_mode);
opts.sensi_mode = 'direct';
opts.max_newton = 3;
% opts.

%% Using the NLF_integrator
s = prepare_gnsf(s, opts);
% adj_seed = rand(nx + nu,1);
% adj_seed = zeros(nx + nu,1);
% adj_seed(s.nx+2) = 1;
% adj_seed(s.nx1+s.nx2) = 1;
adj_seed = [ones(nx,1); zeros(nu,1)]; %adj_seed;ones(s.nx+s.nu,1);
opts.adj_seed = adj_seed;
u0 = controls_MPC(:,1);
x0 = zeros(s.nx,1);
x0(3) = 0.8;
tic
output = proto_gnsf(s, opts, x0, u0);
toc
load('gnsf1_output.mat')
sim_error = max(abs(output.xf - output.xf))
sens_error = max(max(abs(output.forw_sensi - output.forw_sensi)))
sens_adj_error = max(max(abs(output.adj_sensi - output.adj_sensi)))
error_z = max(abs(output.Zf - output.Zf))

keyboard
% output = nsf_dae_lo(s, opts, x0, u0);

sim_error = max(abs(output.xf - output2.xf));
sens_error = max(max(abs(output.forw_sensi - output2.forw_sensi)));
sens_adj_error = max(max(abs(output.adj_sensi - output2.adj_sensi)));
error_z = max(abs(output.Zf - output2.Zf));
if max([sens_error, sim_error, sens_adj_error, error_z])> 1E-10
    error('GO BACK TO OLD VERSION')
end

reference_time = toc;
xfref = output.xf;
sensi1 = output.forw_sensi;
adj_sensi_ref = output.adj_sensi;
Zf_ref = output.Zf;
forw_sensi = [output.forw_sensi; zeros(s.nu,s.nu + s.nx)];
forw_sensi(s.nx+1:s.nx+s.nu, s.nx+1:s.nx+s.nu) = eye( s.nu);
forw_times_adseed = forw_sensi' * adj_seed;

% adj_sensi_ref = forw_sensi' * adj_seed;

output_adj = output.adj_sensi;
error_adj = forw_sensi' * adj_seed - output.adj_sensi;
% table(output_adj, forw_times_adseed, error_adj)

table(le_max(output_adj), le_max(forw_times_adseed), le_max(error_adj))

%% set simulation options
imax = 20;
for i = 1:imax
    n_steps = i;
    q = 4; % number of stages;
%     disp(['Number of stages for collocation = ' num2str(q)]);
%     disp(['Number of steps per simulation =   ' num2str(n_steps)]); 
    opts = set_integrator_opts(q, stepsize, n_steps, forw_sens_mode, adj_sens_mode);
    opts.sensi_mode = 'direct';
    opts.max_newton = 10;
    opts.adj_seed = adj_seed;

    %% Using the NLF_integrator
    s = prepare_struct_exploiting_sim(s, opts);
    tic
    [ output ] = proto_gnsf(s, opts, x0, u0);
    time( i ) = output.time_forw + output.time_adj;
    xf2 = output.xf;
    sensi2 = output.forw_sensi;
    adj_sensi = output.adj_sensi;
    Zf = output.Zf;
    error_Zf(i) = max(abs(Zf - Zf_ref));
    error_adj_sensi(n_steps) = max(max(abs( adj_sensi - adj_sensi_ref)));
    error_xf(n_steps) = max(abs(xf2-xfref));
    error_sensi(n_steps) = max(max(abs(sensi1 - sensi2)));
    time_forw(i) = output.time_forw;
    time_adj(i ) = output.time_adj;
end
plot(1:imax, error_xf)
hold on
plot(1:imax, error_sensi, 'r-')
plot(1:imax, error_adj_sensi, 'g-.')
plot(1:imax, error_Zf, 'k:')
legend('error sim', 'error forw sensi', 'error adj sensi', 'error Z');
xlabel('number of IRK steps')
ylabel('error')
set(gca, 'YScale', 'log')

figure
hold on
plot(1:imax, time, '*')
plot(1:imax, time_forw)
plot(1:imax, time_adj)
legend('time', 'time forw', 'time adj')

% disp('===== iterative sensitivity propagation =====')
if 0 
    %% LOOP
    disp('without algebraic structure exploitation')
    simulation_it = x0;
    xk = x0;
    tic
    for k=1:n_simulations
        uk = controls_MPC(:,k);
        [ output, mem ] = nlf_proto_final(s, opts, mem, xk, uk);
        xk = output.xf;
        simulation_it = [simulation_it xk];
    end
    toc

    dffZ_dx1u = output.dffZ_dx1u;
    sens_iterative = output.forw_sensi;
end