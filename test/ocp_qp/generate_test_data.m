% Condensing routine that outputs data against which acados is tested.

% Always give the same output
rng(1);

% LTV system
N = 20;
nx = 4;
nu = 2;

x0 = randn(nx, 1);
save('x0.dat', 'x0', '-ascii', '-double')
A = randn(nx, nx);
save('A.dat', 'A', '-ascii', '-double')
B = randn(nx, nu);
save('B.dat', 'B', '-ascii', '-double')
c = randn(nx, 1);
save('c.dat', 'c', '-ascii', '-double')

A_bar = [-eye(nx),zeros(nx,(N-1)*nx)];
for i=1:N-1
    A_bar = [A_bar; zeros(nx, (i-1)*(nx)), A, -eye(nx), zeros(nx, (N-i-1)*nx)];
end
B_bar = kron(eye(N),B);
c_bar = repmat(c, N, 1);
c_bar(1:nx) = c_bar(1:nx) + A*x0;

transition_vector = -A_bar \ c_bar;
save('transition_vector.dat', 'transition_vector', '-ascii', '-double')
transition_matrix = -A_bar \ B_bar;
save('transition_matrix.dat', 'transition_matrix', '-ascii', '-double')

