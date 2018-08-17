
N = 20;
nx = 4;
nu = 1;


A = [1.00000   0.00000   0.00000   0.00000 ;
  0.05000   1.00000   0.00000   0.00000 ;
  0.00123   0.04933   1.01691   0.67823 ;
  0.00002   0.00123   0.05028   1.01691].';

B = [0.00125   0.05003   0.00157   0.06285 ].';

Q = 1e-10*eye(nx);
R = 1e-4;
q = zeros(nx,1);
r = 0;

x0 = [0; 0; pi; 0];

G = [];
H = [];
f = [];
A_bar = [];
B_bar = [];
c = zeros(N*nx, 1);
for i=0:N-1
    G = [G, zeros(i*nx, nx+nu); zeros(nx, i*(nx+nu)), A, B, -eye(nx)];
    H = blkdiag(H, blkdiag(Q, R));
    f = [f; q; r];
    if i<N-1
        A_bar = [[A_bar; zeros(nx, i*nx)], [zeros(i*nx, nx); eye(nx); -A]];
    end
    B_bar = blkdiag(B_bar, B);
end
H = blkdiag(H, Q);
f = [f; q];
A_bar = [A_bar, [zeros((N-1)*nx, nx); eye(nx)]];
B_bar = [[A; zeros((N-1)*nx, nx)], B_bar];

w = quadprog(H, f, [zeros(1,100), -1.6, 0, 1.28, 0], -0.6, G, zeros(size(G, 1), 1), x0, x0);

G_bar = A_bar \ B_bar;
g_bar = A_bar \ c;

Q_bar = kron(eye(N), Q);
R_bar = kron(eye(N), R);
h = G_bar.' * (repmat(q, N, 1) + Q_bar*g_bar);
H = R_bar + G_bar.' * Q_bar * G_bar;
