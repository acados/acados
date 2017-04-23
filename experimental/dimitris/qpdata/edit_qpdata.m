
clear variables;clc

% TODO(dimitris): S_vertcat
% TODO(dimitris): Test for varying dimensions

%% load raw qp data

load('nx.txt');
load('nu.txt');
load('A_vertcat.txt');
load('B_vertcat.txt');
load('bv_vertcat.txt');

load('Q_vertcat.txt');
load('R_vertcat.txt');
load('qv_vertcat.txt');
load('rv_vertcat.txt');

load('ub_vertcat.txt');
load('lb_vertcat.txt');

N = length(nx)-1;
nvar = length(lb_vertcat);
neq  = sum(nx(2:end));

%% convert to cells

A = {};
B = {};
b = {};

Q = {};
R = {};
q = {};
r = {};

acc_Q = 0;
acc_R = 0;
acc_q = 0;
acc_r = 0;
acc_A = 0;
acc_B = 0;
acc_b = 0;

for ii = 1:N
    
    Q{ii} = reshape(Q_vertcat(acc_Q+1:acc_Q + nx(ii)*nx(ii)), nx(ii), nx(ii));
    R{ii} = reshape(R_vertcat(acc_R+1:acc_R + nu(ii)*nu(ii)), nu(ii), nu(ii));
    q{ii} = qv_vertcat(acc_q+1:acc_q + nx(ii));
    r{ii} = rv_vertcat(acc_r+1:acc_r + nu(ii));
    
    A{ii} = reshape(A_vertcat(acc_A+1:acc_A + nx(ii)*nx(ii)), nx(ii), nx(ii));
    B{ii} = reshape(B_vertcat(acc_B+1:acc_B + nx(ii)*nu(ii)), nx(ii), nu(ii));
    b{ii} = bv_vertcat(acc_b+1:acc_b + nx(ii));
 
    acc_Q = acc_Q + nx(ii)*nx(ii);
    acc_R = acc_R + nu(ii)*nu(ii);
    acc_q = acc_q + nx(ii);
    acc_r = acc_r + nu(ii);
    acc_A = acc_A + nx(ii)*nx(ii);
    acc_B = acc_B + nx(ii)*nu(ii);
    acc_b = acc_b + nx(ii);
end

Q{end+1} = reshape(Q_vertcat(acc_Q+1:acc_Q + nx(end)*nx(end)), nx(end), nx(end));
q{end+1} = qv_vertcat(acc_q+1:acc_q + nx(end));

%% build quadprog data

H = [];
h = [];
for ii = 1:N
    H = blkdiag(H,Q{ii});
    H = blkdiag(H,R{ii});
    h = [h; q{ii}];
    h = [h; r{ii}];
end
H = blkdiag(H,Q{end});
h = [h; q{end}];

Aeq  = sparse(neq,nvar);
beq  = zeros(neq,1);

for ii = 2:N+1
   
    % row index
    rind = (ii-2)*nx(ii)+1:(ii-1)*nx(ii);
    
    % column indicies TODO
    cind_xp = (ii-2)*(nx(ii-1)+nu(ii-1))+1:(ii-2)*(nx(ii-1)+nu(ii-1))+nx(ii-1);
    cind_up = (ii-2)*(nx(ii-1)+nu(ii-1))+nx(ii-1)+1:(ii-1)*(nx(ii-1)+nu(ii-1));
    cind_x  = (ii-1)*(nx(ii-1)+nu(ii-1))+1:(ii-1)*(nx(ii-1)+nu(ii-1))+nx(ii);
    
    Aeq(rind,cind_xp) = A{ii-1};
    Aeq(rind,cind_up) = B{ii-1};
    Aeq(rind,cind_x)  = -eye(nx(ii));

    beq(rind) = -b{ii-1};
    
end

%% solve with quadprog

opts = optimoptions('quadprog','Display','off');
%opts.TolCon = 1e-20;
%opts.TolFun = 1e-20;
%opts.TolX   = 1e-20;

[sol, fval, flag, output, lam_quadprog] = quadprog(H,h,[],[],Aeq,beq,lb_vertcat,ub_vertcat,[], opts);

xopt = {};
uopt = {};

acc = 0;
for ii = 1:N
    xopt{ii} = sol(acc+1:acc+nx(ii));
    uopt{ii} = sol(acc+nx(ii)+1:acc+nx(ii)+nu(ii));
    acc = acc + nx(ii) + nu(ii);
end
xopt{end+1} = sol(end-nx(end)+1:end);

%% print first input

for ii = 1:nu(ii)
    fprintf('u0[0] = %f\n', uopt{1}(ii));
end
