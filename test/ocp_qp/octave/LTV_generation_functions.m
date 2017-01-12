
% Let randn always return the same output
randn('state', 0);

function [N, nx, nu, nb, nc] = generate_dimensions()
    N = 20;
    nx = 4;
    nu = 2;
    nb = nx+nu;
    nc = nx+nu;
endfunction

function [A_all, B_all, b_all, x0] = generate_dynamics()
    [N, nx, nu] = generate_dimensions();
    A_all = [];
    B_all = [];
    b_all = [];
    for k=0:N-1
        do
            A = randn(nx, nx);
        until(all(abs(eig(A)) < 0.7)) % Unstable dynamics leads to unstable condensing
        A_all = [A_all, A];
        B_all = [B_all, randn(nx, nu)];
        b_all = [b_all, randn(nx, 1)];
    end
    x0 = randn(nx, 1);
endfunction

function [Q_all, S_all, R_all, q_all, r_all] = generate_cost_function()
    [N, nx, nu] = generate_dimensions();
    Q_all = [];
    S_all = [];
    R_all = [];
    q_all = [];
    r_all = [];
    for k=0:N-1
        do
            Q = randn(nx, nx);
            Q = Q.'*Q;
            S = randn(nu, nx);
            R = randn(nu, nu);
            R = R.'*R;
        until(all(eig([Q, S.'; S, R]) > 1e-2)) % Implies a convex QP
        Q_all = [Q_all, Q];
        S_all = [S_all, S];
        R_all = [R_all, R];
        q_all = [q_all, randn(nx, 1)];
        r_all = [r_all, randn(nu, 1)];
    end
    Q = randn(nx, nx);
    Q = Q.'*Q;
    Q_all = [Q_all, Q];
    q_all = [q_all, randn(nx, 1)];
endfunction

function [xl_all, xu_all, ul_all, uu_all] = generate_bounds(x0)
    [N, nx, nu, nb] = generate_dimensions();
    xl_all = [];
    xu_all = [];
    ul_all = [];
    uu_all = [];
    for k=0:N-1
        xl_all = [xl_all, -10*abs(randn(nx,1))];
        xu_all = [xu_all, +10*abs(randn(nx,1))];
        ul_all = [ul_all, -20*abs(randn(nu,1))];
        uu_all = [uu_all, +20*abs(randn(nu,1))];
    end
    xl_all = [xl_all, -10*abs(randn(nx,1))];
    xu_all = [xu_all, +10*abs(randn(nx,1))];

endfunction

function [Cx_all, Cu_all, cl_all, cu_all] = generate_general_constraints()
    [N, nx, nu, ~, nc] = generate_dimensions();
    Cx_all = [];
    Cu_all = [];
    cl_all = [];
    cu_all = [];
    for k=1:N
        Cx_all = [Cx_all, randn(nc, nx)];
        Cu_all = [Cu_all, randn(nc, nu)];
        cl_all = [cl_all, -20*abs(randn(nc, 1))];
        cu_all = [cu_all, +20*abs(randn(nc, 1))];
    end
    Cx_all = [Cx_all, randn(nc, nx)];
    cl_all = [cl_all, -20*abs(randn(nc, 1))];
    cu_all = [cu_all, +20*abs(randn(nc, 1))];
endfunction
