import casadi.*
close all
clear variables

% Time horizon
Ts = 0.1;
Q = diag([1, 1]);
R = 0.05;
x0 = [0.5; 0].';

% Declare model variables
x1 = MX.sym('x1');
x2 = MX.sym('x2');
x = [x1; x2];
u1 = MX.sym('u1');
u = [u1];

nx = 2;
nu = 1;

% Model equations
xdot = [];
xdot=[ xdot; x2+u*(0.5+0.5*x1)];
xdot=[xdot;x1+u*(0.5-2*x2)];

% Continuous time dynamics
f = Function('f', {x, u}, {xdot});

% Control discretization
N = 13; % number of control intervals
M = 50; % RK4 steps per interval
T = Ts*N;
DT = Ts/M;
X0 = MX.sym('X0', 2);
U = MX.sym('U',1);
X = X0;
for j=1:M
    [k1] = f(X, U);
    [k2] = f(X + DT/2 * k1, U);
    [k3] = f(X + DT/2 * k2, U);
    [k4] = f(X + DT * k3, U);
    X=X+DT/(6)*(k1 +2*k2 +2*k3 +k4);
end
F = Function('F', {X0, U}, {X});
W = [X0;U];
JF = Function('JF',{X0,U},{jacobian(X,W)});

% Start with an empty NLP
w={};
w0 = [];
lbw = [];
ubw = [];
J = 0;
g={};
lbg = [];
ubg = [];

% "Lift" initial conditions
X0 = MX.sym('X0', 2);
w = {w{:}, X0};
lbw = [lbw; x0(1); x0(2)];
ubw = [ubw; x0(1); x0(2)];
w0 = [w0; x0(1); x0(2)];

% Formulate the NLP
Xk = X0;
for k=0:N-1
    % New NLP variable for the control
    Uk = MX.sym(['U_' num2str(k)],1);
    w = {w{:}, Uk};
    lbw = [lbw; -Inf];
    ubw = [ubw;  Inf];
    w0 = [w0;  0];

    % Integrate till the end of the interval
    [Xk_end] = F(Xk, Uk);

    % New NLP variable for state at end of interval
    Xk = MX.sym(['X_' num2str(k+1)], 2);
    J=J+0.5*(Xk(1)^2 + Xk(2)^2 + 0.05*Uk(1)^2);

    w = {w{:}, Xk};
    lbw = [lbw; -Inf; -inf];
    ubw = [ubw;  inf;  inf];
    w0 = [w0; 0; 0];
        
    % Add equality constraint
    g = {g{:}, Xk_end-Xk};
    lbg = [lbg; 0; 0];
    ubg = [ubg; 0; 0];
end

% Create an NLP solver
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
opts = [];
opts.ipopt.fixed_variable_treatment= 'make_constraint';
solver = nlpsol('solver', 'ipopt', prob);

% Solve the NLP
sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,'lbg', lbg, 'ubg', ubg);
w_opt = full(sol.x);

% Plot the solution
x1_opt = w_opt(1:3:end);
x2_opt = w_opt(2:3:end);
u1_opt = w_opt(3:3:end);

sol_x = [x1_opt x2_opt];
sol_x.'



tgrid = linspace(0, T, N+1);
clf;
subplot(211)
hold on
plot(tgrid, x1_opt, '--')
plot(tgrid, x2_opt, '-')
grid on
subplot(212)
hold on
stairs(tgrid, [u1_opt; nan], '-.')
grid on
xlabel('t')
legend('x1','x2','u')

sol_x = [x1_opt x2_opt]
sol_u = u1_opt

sol_pi  = full(sol.lam_g)

sol_lam_x  = full(sol.lam_x)

sol_pi = reshape(full(sol.lam_g),nx,N);
% Compute residuals
res_stat = Inf*ones((nx+nu),N+1);
J = full(JF(sol_x(1,:),sol_u(1,:)));
B = J(:,nx+1:end);
res_stat(nx+1:end,1) = R*sol_u(1,:) + B.'*sol_pi(:,1);
for i = 2:N
    J = full(JF(sol_x(i,:),sol_u(i,:)));
    A = J(:,1:nx);
    B = J(:,nx+1:end);
    res_stat(1:nx,i) = Q*sol_x(i,:).' +  A.'*sol_pi(:,i) - sol_pi(:,i-1);
    res_stat(nx+1:end,i) = R*sol_u(i,:).' +  B.'*sol_pi(:,i);
end

i = N+1;

J = full(JF(sol_x(i,:),ones(nu,1)));
    A = J(:,1:nx);
    res_stat(1:nx,i) = Q*sol_x(i,:).' - sol_pi(:,i-1);
