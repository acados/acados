minimal_example_ocp;
u_ref = ocp.get('u');
x_ref = ocp.get('x');

ocp.generate_c_code

cd c_generated_code/

% templated MEX
t_ocp = pendulum_mex_solver;
t_ocp.solve
t_ocp.print
u_t = t_ocp.get('u');
x_t = t_ocp.get('x');


% comparison
err_u = max(abs(u_ref - u_t))
err_ = max(max(abs(x_ref - x_t)))
cd ..
keyboard
clear all

