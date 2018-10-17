import casadi.*
import acados.*

%CasadiMeta.version()

nx = 4;
%nu = 0;

x = SX.sym('x', nx, 1);
xdot = SX.sym('xdot', nx, 1);

expl_ode_expr = -2*x;
impl_ode_expr = xdot - expl_ode_expr;

sim_model = acados_integrator_model();
%sim_model = sim_model.set('type', 'explicit');
sim_model = sim_model.set('type', 'implicit');
%sim_model = sim_model.set('ode_expr', expl_ode_expr);
sim_model = sim_model.set('ode_expr', impl_ode_expr);
sim_model = sim_model.set('x', x);
sim_model = sim_model.set('xdot', xdot);
%sim_model = sim_model.set('model_name', 'expl_model');
sim_model = sim_model.set('model_name', 'impl_model');

sim_model



sim_opts = acados_integrator_opts();
%sim_opts = sim_opts.set('scheme', 'erk');
sim_opts = sim_opts.set('scheme', 'irk');
sim_opts = sim_opts.set('sens_forw', 'true');
%sim_opts = sim_opts.set('sens_forw', 'false');
sim_opts = sim_opts.set('codgen_model', 'true');
%sim_opts = sim_opts.set('codgen_model', 'false');

sim_opts



tic
sim = acados_integrator(sim_model, sim_opts);
sim_create_time = toc

sim



tic
flag = sim.solve()
sim_solve_time = toc
