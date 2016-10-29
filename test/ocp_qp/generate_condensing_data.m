% Condensing routine that outputs data against which acados is tested.
clear
% pkg install -forge struct
% pkg install -forge optim
pkg load optim

% Let randn always return the same output
randn('state', 0);

% LTI system
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
do
    Q = randn(nx, nx);
    Q = Q.'*Q;
    S = randn(nu, nx);
    R = randn(nu, nu);
    R = R.'*R;
until(all(eig([Q, S.'; S, R]) > 1e-2))

save('Q.dat', 'Q', '-ascii', '-double');
save('S.dat', 'S', '-ascii', '-double');
save('R.dat', 'R', '-ascii', '-double');
q = randn(nx, 1);
save('qv.dat', 'q', '-ascii', '-double');
r = randn(nu, 1);
save('rv.dat', 'r', '-ascii', '-double');

A_bar = [-eye(nx),zeros(nx,(N-1)*nx)];
for i=1:N-1
    A_bar = [A_bar; zeros(nx, (i-1)*(nx)), A, -eye(nx), zeros(nx, (N-i-1)*nx)];
end
B_bar = kron(eye(N),B);
b_bar = repmat(b, N, 1);

g = -A_bar \ (b_bar + [A*x0; zeros((N-1)*nx,1)]);
save('transition_vector.dat', 'g', '-ascii', '-double');
G = -A_bar \ B_bar;
save('transition_matrix.dat', 'G', '-ascii', '-double');

Q_bar = kron(eye(N), Q);
S_bar = kron(eye(N), S.');
R_bar = kron(eye(N), R);
q_bar = repmat(q, N, 1);
r_bar = repmat(r, N, 1);

% Test condensing
H = blkdiag(kron(eye(N), [Q, S.'; S, R]), Q);
h = [repmat([q;r], N, 1);q];
G_bar = [-eye(nx),zeros(nx,N*(nx+nu))];
for i=1:N
    G_bar = [G_bar; zeros(nx, (i-1)*(nx+nu)), A, B, -eye(nx), zeros(nx, (N-i)*(nx+nu))];
end
g_bar = [x0; b_bar];

ubx = 5*(x0 + abs(randn(nx,1)));
ubu = 5*abs(randn(nu,1));
lbx = -ubx;
lbu = -ubu;
save('upper_bound_x.dat', 'ubx', '-ascii', '-double');
save('lower_bound_x.dat', 'lbx', '-ascii', '-double');
save('upper_bound_u.dat', 'ubu', '-ascii', '-double');
save('lower_bound_u.dat', 'lbu', '-ascii', '-double');
ubw = [repmat([ubx;ubu], N, 1); ubx];
lbw = [repmat([lbx;lbu], N, 1); lbx];

[w_star_ocp, ~, exit_flag, ~, lambda_star_ocp] = quadprog(H, h, [], [], G_bar, -g_bar, lbw, ubw);
if(~(exit_flag == 1))
    Z = null(G_bar);
    disp(['convex QP? : ', num2str(all(eig(Z.'*H*Z) > 1e-4))])
    error(['QP solution failed with code: ', num2str(exit_flag)]);
end
% Octave follows a different convention from us
lambda_star_ocp = -lambda_star_ocp.eqlin;

KKT_system_ocp = [];
for i=0:N-1
    xk = w_star_ocp(i*(nx+nu)+1:i*(nx+nu)+nx);
    uk = w_star_ocp(i*(nx+nu)+nx+1:(i+1)*(nx+nu));
    lambdak = lambda_star_ocp(i*nx+1:(i+1)*nx);
    lambdak_1 = lambda_star_ocp((i+1)*nx+1:(i+2)*nx);
    KKT_system_ocp = [KKT_system_ocp; Q*xk + q + S.'*uk - lambdak + A.'*lambdak_1];
    KKT_system_ocp = [KKT_system_ocp; R*uk + r + S*xk + B.'*lambdak_1];
end
lambdaN = lambda_star_ocp(N*(nx)+1:end);
xN = w_star_ocp(N*(nx+nu)+1:end);
KKT_system_ocp = [KKT_system_ocp; Q*xN + q - lambdaN];

% There is a bug somewhere in the following lines
% h_bar = r_bar + G.'*(q_bar + Q_bar*g) + S_bar.'*[g(1:end-nx);zeros(nx,1)];
% h_bar = r_bar + G.'*(q_bar + Q_bar*g) + S_bar.'*g;
% save('condensed_gradient.dat', 'h_bar', '-ascii', '-double');
% S_cut = blkdiag(kron(eye(N-1), S.'), zeros(nx, nu));
% G_cut = [G(1:end-nx, :); zeros(nx, N*nu)];
% H_bar = R_bar + G.'*Q_bar*G + S_cut.'*G + G.'*S_cut;
% save('condensed_hessian.dat', 'H_bar', '-ascii', '-double');

% Janicks' way
c = b_bar;
L = [A; zeros((N-1)*nx, nx)];
G2 = [zeros(nx, N*(nu)); -A_bar \ B_bar];
g2 = [zeros(nx, 1); -A_bar \ c];
Ge = [eye(nx); -A_bar \ L];
Q_all = kron(eye(N+1), Q);
S_all = [kron(eye(N), S), zeros(N*nu, nx)];
R_all = kron(eye(N), R);
r_all = repmat(r, N, 1);
q_all = repmat(q, N+1, 1);
Se = G2.'*Q_all*Ge + S_all*Ge;

r_condensed = r_all + G2.'*(q_all + Q_all*g2) + S_all*g2 + Se*x0;
R_condensed = R_all + G2.'*Q_all*G2 + S_all*G2 + G2.'*S_all.';
save('condensed_gradient.dat', 'r_condensed', '-ascii', '-double');
save('condensed_hessian.dat', 'R_condensed', '-ascii', '-double');

Dx = blkdiag(kron(eye(N), [eye(nx); zeros(nu, nx)]), eye(nx));
Du = [kron(eye(N), [zeros(nx, nu); eye(nu)]); zeros(nx, N*nu)];

% Stack bound constraints
d = [ubw; -lbw];
Dx = [Dx; -Dx];
Du = [Du; -Du];
Dex0 = Dx*Ge*x0;

D_bar = Du + Dx*G2;
d_bar = d - Dx*g2 - Dex0;

% w_star_condensed = qp([], H_bar, h_bar);
w_star_condensed_quadprog = quadprog(R_condensed, r_condensed, D_bar, d_bar);

u_lb = repmat(lbu, N, 1);
u_ub = repmat(ubu, N, 1);
A_lb = repmat(lbx, N+1, 1) - g2 - Ge*x0;
A_ub = repmat(ubx, N+1, 1) - g2 - Ge*x0;
% First nx rows of G2 are zero
A_lb = A_lb(nx+1:end);
A_ub = A_ub(nx+1:end);
A_condensed = G2(nx+1:end,:);
w_star_condensed = qp([], R_condensed, r_condensed, [], [], u_lb, u_ub, A_lb, A_condensed, A_ub);

save('u_lower_bound.dat', 'u_lb', '-ascii', '-double');
save('u_upper_bound.dat', 'u_ub', '-ascii', '-double');
save('condensed_lower_bound.dat', 'A_lb', '-ascii', '-double');
save('condensed_upper_bound.dat', 'A_ub', '-ascii', '-double');
save('condensed_constraint_matrix.dat', 'A_condensed', '-ascii', '-double');

% Compare sparse and condensed solution
XU = reshape([w_star_ocp;zeros(nu,1)], nx+nu, N+1);
u_star_ocp = XU(nx+1:end, 1:end-1);
u_star_ocp = u_star_ocp(:);
disp(['Difference between condensed and sparse solution: ', num2str(norm(u_star_ocp - w_star_condensed))])
