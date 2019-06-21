clear all

import casadi.*
import acados.*

%addpath('../../interfaces/matlab/hpipm_matlab')

%CasadiMeta.version()



fprintf('\nfirst integrator\n');

nx = 4;
%nu = 0;

x = SX.sym('x', nx, 1);
xdot = SX.sym('xdot', nx, 1);

expl_ode_expr = -2*x;
impl_ode_expr = xdot - expl_ode_expr;

tic

sim_model = acados_integrator_model();
sim_model.set('model_name', 'expl_model');
%sim_model.set('model_name', 'impl_model');
sim_model.set('type', 'explicit');
%sim_model.set('type', 'implicit');
sim_model.set('ode_expr', expl_ode_expr);
%sim_model.set('ode_expr', impl_ode_expr);
sim_model.set('x', x);
%sim_model.set('xdot', xdot);

sim_model_create_time = toc;
fprintf('sim_model time %e\n', sim_model_create_time);

% sim_model


tic

sim_opts = acados_integrator_opts();
sim_opts.set('scheme', 'erk');
%sim_opts.set('scheme', 'irk');
%sim_opts.set('sens_forw', 'true');
sim_opts.set('sens_forw', 'false');
sim_opts.set('codgen_model', 'true');
%sim_opts.set('codgen_model', 'false');

sim_opts_create_time = toc;
fprintf('sim_opts time %e\n', sim_opts_create_time);

%sim_opts



tic

sim = acados_integrator(sim_model, sim_opts);

sim_create_time = toc;
fprintf('sim create time %e\n', sim_create_time);

%sim



tic

x0 = [1; 0; 2; -1];
xdot0 = x0 + 2*x0;
sim.set('x', x0);
sim.set('xdot', xdot0);
sim.set('t', 0.05);

sim_set_x_time = toc;
fprintf('sim set x time %e\n', sim_set_x_time);



tic

flag = sim.solve();

sim_solve_time = toc;
fprintf('sim solve time %e\n', sim_solve_time);



tic

xn = sim.get('xn');
Sxn = sim.get('Sxn');

sim_get_xn_time = toc;
fprintf('sim get xn time %e\n', sim_get_xn_time);
xn
Sxn





fprintf('\nsecond integrator\n');

% update expression
expl_ode_expr = -1*x;
impl_ode_expr = xdot - expl_ode_expr;

% update model
%sim_model.set('model_name', 'new_expl_model');
%sim_model.set('type', 'explicit');
sim_model.set('ode_expr', expl_ode_expr);

% update opts
%sim_opts.set('scheme', 'erk');
sim_opts.set('sens_forw', 'true');
%sim_opts.set('codgen_model', 'true');

% new sim solver
new_sim = acados_integrator(sim_model, sim_opts);

new_sim.set('x', x0);
new_sim.set('xdot', xdot0);
new_sim.set('t', 0.05);

flag = new_sim.solve();

xn = new_sim.get('xn');
Sxn = new_sim.get('Sxn');

xn
Sxn
