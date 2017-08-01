import acados.*

qp = ocp_qp(struct('N',5, 'nx',2, 'nu', 1));

% specify OCP
qp.A = [0 0; 1 0];
qp.B = [0; 1];
qp.Q = eye(2);
qp.R = 1;

% specify initial condition
x0 = [1; 1];
qp.lb{1} = x0;
qp.ub{1} = x0;

% solve QP
solver = ocp_qp_solver('qpdunes', qp);
output = solver.evaluate();
assert(abs(-0.5 - output.controls{1}) < 1e-8)
disp(output.states)