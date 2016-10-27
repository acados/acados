#!/usr/bin/octave -qf

% Condensing routine that outputs data against which acados is tested.
clear

% Let randn always return the same output
randn('state', 0);

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
save('bv.dat', 'b', '-ascii', '-double');
Q = randn(nx, nx);
save('Q.dat', 'Q', '-ascii', '-double');
S = 0*randn(nu, nx);
save('S.dat', 'S', '-ascii', '-double');
R = 0.01*eye(nu);randn(nu, nu);
save('R.dat', 'R', '-ascii', '-double');
q = 0*randn(nx, 1);
save('qv.dat', 'q', '-ascii', '-double');
r = 0*randn(nu, 1);
save('rv.dat', 'r', '-ascii', '-double');

A_bar = [-eye(nx),zeros(nx,(N-1)*nx)];
for i=1:N-1
    A_bar = [A_bar; zeros(nx, (i-1)*(nx)), A, -eye(nx), zeros(nx, (N-i-1)*nx)];
end
B_bar = kron(eye(N),B);
b_bar = repmat(b, N, 1);
b_bar(1:nx) = b_bar(1:nx) + A*x0;

g = -A_bar \ b_bar;
save('transition_vector.dat', 'g', '-ascii', '-double');
G = -A_bar \ B_bar;
save('transition_matrix.dat', 'G', '-ascii', '-double');

Q_bar = kron(eye(N), Q);
S_bar = kron(eye(N), S.');
R_bar = kron(eye(N), R);
q_bar = repmat(q, N, 1);
r_bar = repmat(r, N, 1);

h_bar = r_bar + transpose(G)*(q_bar + Q_bar*g) + S_bar.'*[g(1:end-nu);zeros(nu,1)];
h_bar = r_bar + G.'*(q_bar + Q_bar*g); + S_bar.'*g;
save('condensed_gradient.dat', 'h_bar', '-ascii', '-double');
H_bar = R_bar + G.'*Q_bar*G + S_bar.'*G + G.'*S_bar;
save('condensed_hessian.dat', 'H_bar', '-ascii', '-double');

% Test condensing
% H = blkdiag(kron(eye(N), [Q, S.'; S, R]), Q);
% h = [repmat([q;r], N, 1);q];
% G_bar = [-eye(nx),zeros(nx,N*(nx+nu))];
% for i=1:N
%     G_bar = [G_bar; zeros(nx, (i-1)*(nx+nu)), A, B, -eye(nx), zeros(nx, (N-i)*(nx+nu))];
% end
% g_bar = [-x0; zeros(N*nx, 1)];
% w_star_ocp = qp([], H, h, G_bar, g_bar);
% XU = reshape([w_star_ocp;zeros(nu,1)], nx+nu, N+1);
% u_star_ocp = XU(nx+1:end, 1:end-1);
% u_star_ocp = u_star_ocp(:);
% w_star_condensed = qp([], H_bar, h_bar);
%
% [u_star_ocp w_star_condensed]
% norm(u_star_ocp - w_star_condensed)
