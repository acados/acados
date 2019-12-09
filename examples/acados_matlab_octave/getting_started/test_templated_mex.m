minimal_example_ocp;
u_ref = ocp.get('u');

ocp.generate_c_code

cd c_generated_code/

% templated MEX
t_ocp = pendulum_mex_solver;
t_ocp.solve
t_ocp.print
u_t = t_ocp.get('u');


% comparison
err_u = max(abs(u_ref - u_t))
cd ..
keyboard
clear all

