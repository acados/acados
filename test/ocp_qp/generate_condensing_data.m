% Condensing routine that outputs data against which acados is tested.

% Let randn always return the same output
randn("state", 0);

% LTV system
N = 20;
nx = 4;
nu = 2;

save('N.dat', 'N', '-ascii', '-double');
save('nx.dat', 'nx', '-ascii', '-double');
save('nu.dat', 'nu', '-ascii', '-double');
x0 = randn(nx, 1);
save('x0.dat', 'x0', '-ascii', '-double');
A = randn(nx, nx);
save('A.dat', 'A', '-ascii', '-double');
B = randn(nx, nu);
save('B.dat', 'B', '-ascii', '-double');
b = randn(nx, 1);
save('b.dat', 'b', '-ascii', '-double');

A_bar = [-eye(nx),zeros(nx,(N-1)*nx)];
for i=1:N-1
    A_bar = [A_bar; zeros(nx, (i-1)*(nx)), A, -eye(nx), zeros(nx, (N-i-1)*nx)];
end
B_bar = kron(eye(N),B);
b_bar = repmat(b, N, 1);
b_bar(1:nx) = b_bar(1:nx) + A*x0;

transition_vector = -A_bar \ b_bar;
save('transition_vector.dat', 'transition_vector', '-ascii', '-double');
transition_matrix = -A_bar \ B_bar;
save('transition_matrix.dat', 'transition_matrix', '-ascii', '-double');
