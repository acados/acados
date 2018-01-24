import acados.*

qp = ocp_qp(5, 2, 1);

% specify OCP
qp.update('A', [0 1; 0 0]);
qp.update('B', [0; 1]);
qp.update('Q', eye(2));
qp.update('R', 1);

% specify initial condition
x0 = [1; 1];
qp.update('lbx', 0, x0);
qp.update('ubx', 0, x0);

% solve QP
disp(qp);
solver = ocp_qp_solver(FULL_CONDENSING_QPOASES, qp);
output = solver.evaluate(qp);

xopt = output.states();

for i=1:numel(xopt)
    disp(xopt{i})
end
