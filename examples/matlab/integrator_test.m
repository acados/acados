import acados.*
import casadi.*


pendulum_ode = pendulum_model();

opts1 = struct('step', 0.1, ...% only mandatory argument
               'sens_forw', 1, ...
               'sens_adj', 1, ...
               'jac_reuse', 0);
sim1 = acados.integrator(pendulum_ode, opts1);

opts2 = struct('step', 0.2, 'use_MX', true);
chen_ode = chen_model();
sim2 = acados.integrator(chen_ode, opts2);


x1 = [ 4.2; 0.42; 0; 0];
u = 1;

x2 = [ 4.2; 0];

t1 = 0.0;
t2 = 0.0;
M = 100;

for k=1:M
    
    tic
    x1 = cell2mat(sim1.integrate(x1, u));
    t1 = t1 + toc;
    
    tic
    x2 = cell2mat(sim2.integrate(x2, u));
    t2 = t2 + toc;
    
    sim1.set_step(sim1.step()*0.99);
    disp(x1)
    disp(x2)
end

disp(['avg time sim1: ', num2str(t1/M)]);
disp(['avg time sim2: ', num2str(t2/M)]);
