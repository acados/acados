clc;
clear all;
close all;

import casadi.*
mkdir('chain');

% Get collocation points
method = '';
Ns = 2; % NUMBER OF INTEGRATION STEPS

resX = []; resU = [];
for d = 0:3
for Nm = 2:4
    disp(['---- Nm value = ' num2str(Nm) '----']);

% Environment
g = 9.81;     % [N/kg]
L = 0.033;
D = 1.0;
m = 0.03;
x0 = zeros(3,1);
xN = [1 0 0].';

wall_pos = -0.01;

T = 3.0;
N = 15;

% Number of variables
nx = (Nm-1)*2*3;
nu = 3;

% Trivial LS cost function
cost_x = SX.sym('x',nx);
cost_u = SX.sym('u',nu);
cost_xu = [cost_x;cost_u];
cost_y = [cost_x;cost_u];
cost_jac_y = jacobian(cost_y,cost_xu);
ls_cost = Function(['ls_cost_nm' num2str(Nm)], {cost_x,cost_u},{cost_y,cost_jac_y});

opts = struct('mex', false);
ls_cost.generate(['ls_cost_nm' num2str(Nm)], opts);

% Trivial terminal LS cost function
cost_yN = cost_x;
cost_jac_yN = jacobian(cost_yN,cost_x);
ls_costN = Function(['ls_costN_nm' num2str(Nm)], {cost_x}, {cost_yN, cost_jac_yN});
ls_costN.generate(['ls_costN_nm' num2str(Nm)], opts);

% Empty path constraints
pathcon_g = SX.zeros(0);
pathcon_jac_g = jacobian(pathcon_g, cost_xu);
pathcon = Function(['pathcon_nm' num2str(Nm)], {cost_x, cost_u}, {pathcon_g, pathcon_jac_g});
pathcon.generate(['pathcon_nm' num2str(Nm)], opts);

% Empty terminal path constraints
pathcon_gN = SX.zeros(0);
pathcon_jac_gN = jacobian(pathcon_gN, cost_x);
pathconN = Function(['pathconN_nm' num2str(Nm)], {cost_x}, {pathcon_gN, pathcon_jac_gN});
pathconN.generate(['pathconN_nm' num2str(Nm)], opts);

% State variables
u = SX.sym('u',3);
dae.p = u;

dae.x = [];
states = [];
for i = 1:Nm-1
    p = SX.sym(['p' num2str(i)],3);
    v = SX.sym(['v' num2str(i)],3);

    x_struct = struct('p',p,'v',v);
    states = [states; x_struct];
    dae.x = [dae.x; casadi_struct2vec(x_struct)];
end

% Compute forces
F = {};
for i = 1:Nm-1
    if i == 1
        dist = states(1).p-x0;
    else
        dist = states(i).p-states(i-1).p;
    end
    tmp = D*(1 - L/sqrt(dist.'*dist));
    F = {F{:}, tmp*dist};
end

% Set up ODE
dae.ode = [];
for i = 1:Nm-2
    f = 1/m*(F{i+1} - F{i}) - [0;0;g];
    dae.ode = [dae.ode; casadi_vec(x_struct,'p',states(i).v,'v',f)];
end
dae.ode = [dae.ode; casadi_vec(x_struct,'p',states(end).v,'v',u)];

if d > 0
tau_root = collocation_points(d,'legendre');

collfun = simpleColl(dae,tau_root,T/(Ns*N));
collfun = collfun.expand();
else
expl_ode = Function('ode',{dae.x,dae.p},{dae.ode});
RK4fun = simpleRK(expl_ode, Ns, 4);
end

%% Find rest position
Xpoints = linspace(0,1,Nm);
x0_guess = [Xpoints(2:end);zeros(5,Nm-1)];
x0_guess = x0_guess(:);
u_guess = zeros(3,1);

odeFun = Function(['ode_chain_nm' num2str(Nm)],{dae.x,dae.p},{dae.ode,jacobian(dae.ode,[dae.x;dae.p])+SX.zeros(nx,nx+nu)});

Sx = SX.sym('Sx',nx,nx);
Sp = SX.sym('Sp',nx,nu);

Sx_dot = SX.sym('Sx_dot',nx,nx);
Sp_dot = SX.sym('Sp_dot',nx,nu);

Sx_full = [Sx;Sx_dot];
Sp_full = [Sp;Sp_dot];

x_dot = SX.sym('x_dot',nx,1);

x_impl = [dae.x;x_dot];

ode_impl = x_dot - dae.ode;

vdeX = Sx_dot;
vdeX = vdeX - jtimes(dae.ode,dae.x,Sx);

vdeP = Sp_dot - jacobian(dae.ode,dae.p);
vdeP = vdeP - jtimes(dae.ode,dae.x,Sp);

vdeFun = Function(['vde_chain_nm' num2str(Nm)],{x_impl,Sx_full,Sp_full,dae.p},{ode_impl,vdeX,vdeP});

jacX = SX.zeros(nx,2*nx) + jacobian(ode_impl,x_impl);
jacFun = Function(['jac_chain_nm' num2str(Nm)],{x_impl,dae.p},{ode_impl,jacX});

opts = struct('mex', false);
vdeFun.generate(['vde_chain_nm' num2str(Nm)], opts);
jacFun.generate(['jac_chain_nm' num2str(Nm)], opts);

val = odeFun(x0_guess,u_guess);
while norm(full(val)) > 1e-10
    [val,jac] = odeFun(x0_guess,u_guess);
    jac = full(jac); val = full(val);
    delta = -jac\val;
    x0_guess = x0_guess + delta(1:nx);
end
% x0_guess
xN_term = x0_guess;
err_rest = norm(full(val))

x0_mat2 = [zeros(6,1) reshape(x0_guess,6,Nm-1)];

Ypoints = linspace(0,1.5,Nm);
Zpoints = linspace(0,0.5,Nm);

x0_init = [zeros(1,Nm-1);Ypoints(2:end);Zpoints(2:end);zeros(3,Nm-1)];
x0_init = x0_init(:);
u_init = zeros(3,1);

val = odeFun(x0_init,u_init);
while norm(full(val)) > 1e-10
    [val,jac] = odeFun(x0_init,u_init);
    jac = full(jac); val = full(val);
    delta = -jac\val;
    x0_init = x0_init + delta(1:nx);
end
% x0_init
err_rest = norm(full(val))

save(['chain/x0_nm' num2str(Nm) '_d' num2str(d) '.dat'], 'x0_init', '-ascii', '-double');
save(['chain/xN_nm' num2str(Nm) '_d' num2str(d) '.dat'], 'xN_term', '-ascii', '-double');

x0_mat = [zeros(6,1) reshape(x0_init,6,Nm-1)];

Fontsize = 20;
set(0,'DefaultAxesFontSize',Fontsize)

 x0 = [];
 for k = 1:N
     x0 = [x0;repmat(x0_guess,Ns*d+1,1);zeros(nu,1)];
 end
 x0 = [x0;x0_guess];

for rho = [0.7 0.2 0.0]
%% Optimal Control Problem
Xs = {};
for i=1:N+1
   Xs{i} = MX.sym(['X_' num2str(i)],nx);
end
XCs = {};
Us = {};
for i=1:N
   XCs{i} = MX.sym(['XC_' num2str(i)],nx,Ns*d);
   Us{i}  = MX.sym(['U_' num2str(i)],nu);
end

V_block = struct();
V_block.X  = Sparsity.dense(nx,1);
V_block.XC  = Sparsity.dense(nx,Ns*d);
V_block.U  = Sparsity.dense(nu,1);

% Simple bounds on states
lbx = {};
ubx = {};

% List of constraints
g = {};

% List of all decision variables (determines ordering)
V = {};
for k=1:N
  % Add decision variables
  V = {V{:} casadi_vec(V_block,'X',Xs{k},'XC',XCs{k},'U',Us{k})};

  if k==1
    % Bounds at t=0
    x_lb = x0_init;
    x_ub = x0_init;
    u_lb = -10*ones(3,1);
    u_ub = 10*ones(3,1);
    lbx = {lbx{:} casadi_vec(V_block,-inf,'X',x_lb,'U',u_lb)};
    ubx = {ubx{:} casadi_vec(V_block,inf, 'X',x_ub,'U',u_ub)};
  else %k < N
    % Bounds for other t
    m_lb  = [-inf;wall_pos;-inf;-inf;-inf;-inf];
    m_ub  = [inf;inf;inf;inf;inf;inf];
    x_lb = repmat(m_lb,Nm-1,1);
    x_ub = repmat(m_ub,Nm-1,1);
    u_lb = -10*ones(3,1);
    u_ub = 10*ones(3,1);
    lbx = {lbx{:} casadi_vec(V_block,-inf,'X',x_lb,'U',u_lb)};
    ubx = {ubx{:} casadi_vec(V_block,inf, 'X',x_ub,'U',u_ub)};
  end
  % Obtain collocation expressions
  Xcur = Xs{k};
  if d > 0
    for i = 1:Ns
		  [Xcur,coll_out] = collfun(Xcur,XCs{k}(:,1+(i-1)*d:i*d),Us{k});
		  g = {g{:} coll_out};         % collocation constraints
    end
  else
	  Xcur = RK4fun(Xcur,Us{k},T/N);
  end
  g = {g{:} Xs{k+1}-Xcur}; % gap closing
end

V = {V{:} Xs{end}};

% Bounds for final t
x_lb = (1-rho)*xN_term+rho*x0_init;
x_ub = (1-rho)*xN_term+rho*x0_init;
lbx = {lbx{:} x_lb};
ubx = {ubx{:} x_ub};

% Objective function
fun_ref = (vertcat(Xs{:})-repmat(xN_term,N+1,1));
controls = vertcat(Us{:});
effort = 1/2*(controls.'*controls + 1e-2*fun_ref.'*fun_ref);

nlp = struct('x',vertcat(V{:}), 'f',effort, 'g', vertcat(g{:}));

nlpfun = Function('nlp',nlp,char('x','p'),char('f','g'));

solver_opts = struct('print_time', false);
solver_opts.ipopt = struct('print_level', 0, 'linear_solver', 'mumps');
%opts.ipopt = struct('linear_solver','ma27','acceptable_tol', 1e-12, 'tol', 1e-12);
solver = nlpsol('solver','ipopt',nlp, solver_opts);

%args = struct;
%args.x0 = x0;
%args.lbx = vertcat(lbx{:});
%args.ubx = vertcat(ubx{:});
%args.lbg = 0;
%args.ubg = 0;

res = solver('x0',x0,'lbx',vertcat(lbx{:}),'ubx',vertcat(ubx{:}),'lbg',0,'ubg',0);

x0 = full(res.x);
end
struct_res = res;

mu = [];
for k = 1:N
   mu = [mu; full(res.lam_g((k-1)*(Ns*nx*d+nx)+1:(k-1)*(Ns*nx*d+nx)+Ns*nx*d))];
end
lam = [];
for k = 1:N
   lam = [lam; full(res.lam_g((k-1)*(Ns*nx*d+nx)+Ns*nx*d+1:k*(Ns*nx*d+nx)))];
end

dim = size(casadi_struct2vec(V_block));
res_split = vertsplit(res.x,dim(1));

res_U = {}; resK = [];
for r=res_split(1:end-1)
    rs = casadi_vec2struct(V_block,r{1});
    res_U = {res_U{:} rs.U};

    k_mat = full(rs.XC);
    for i = 1:Ns
        resK = [resK; k_mat(:,1+(i-1)*d:i*d)];
    end
end

res_X = {};
for r=res_split(1:end-1)
    rs = casadi_vec2struct(V_block,r{1});
    res_X = {res_X{:} rs.X};
end
res_X = {res_X{:} res_split{end}};

resX = vertcat(res_X{:});
resX = full(reshape(resX,nx,N+1));

resU = vertcat(res_U{:});
resU = full(reshape(resU,nu,N));

res = full(res.x);

save(['chain/resX_nm' num2str(Nm) '_d' num2str(d) '.dat'], 'resX', '-ascii', '-double');
save(['chain/resU_nm' num2str(Nm) '_d' num2str(d) '.dat'], 'resU', '-ascii', '-double');

end
end
