import acados.*

qp = ocp_qp(5, 2, 1);

% specify OCP
qp.set('A', [0 1; 0 0]);
qp.set('B', [0; 1]);
qp.set('Q', eye(2));
qp.set('R', 1);

% specify initial condition
x0 = [1; 1];
qp.set('lbx', 0, x0);
qp.set('ubx', 0, x0);

% solve QP
output = qp.solve('qpoases');

xopt = output.states();

for i=1:numel(xopt)
    disp(xopt{i})
end
