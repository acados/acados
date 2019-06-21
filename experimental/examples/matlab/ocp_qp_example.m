
import acados.*

qp = ocp_qp(5, 2, 1);

% specify OCP
qp.set_field('A', [1 1; 0 1]);
qp.set_field('B', [0; 1]);
qp.set_field('Q', eye(2));
qp.set_field('R', 1);

% specify bounds
qp.set_field('lbx', [0.5; -inf]);
qp.set_field('ubx', [3.0; +inf]);

% specify initial condition
x0 = [1.1; 1.1];
qp.set_field('lbx', 0, x0);
qp.set_field('ubx', 0, x0);

% solve QP
qp.initialize_solver('qpoases');
output = qp.solve();

xopt = output.states();

for i=1:numel(xopt)
    disp(xopt{i})
end
