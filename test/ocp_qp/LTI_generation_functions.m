
% Let randn always return the same output
randn('state', 0);

function [N, nx, nu, nb, nc] = generate_dimensions()
    N = 20;
    nx = 4;
    nu = 2;
    nb = nx+nu;
    nc = nx+nu;
endfunction

function [A, B, b, x0] = generate_dynamics()
    [~, nx, nu] = generate_dimensions();
    do
        A = randn(nx, nx);
    until(all(abs(eig(A)) < 0.9)) % Unstable dynamics leads to unstable condensing
    B = randn(nx, nu);
    b = randn(nx, 1);
    x0 = randn(nx, 1);
endfunction

function [Q, S, R, q, r] = generate_cost_function()
    [~, nx, nu] = generate_dimensions();
    do
        Q = randn(nx, nx);
        Q = Q.'*Q;
        S = randn(nu, nx);
        R = randn(nu, nu);
        R = R.'*R;
    until(all(eig([Q, S.'; S, R]) > 1e-2)) % Implies a convex QP
    q = randn(nx, 1);
    r = randn(nu, 1);
endfunction

function [xl, xu, ul, uu] = generate_bounds(x0)
    [~, nx, nu, nb] = generate_dimensions();
    xl = -5*(abs(x0) + abs(randn(nx,1)));
    xu = +5*(abs(x0) + abs(randn(nx,1)));
    ul = -10*abs(randn(nu,1));
    uu = +10*abs(randn(nu,1));
endfunction

function [Cx, Cu, cl, cu] = generate_general_constraints()
    [~, nx, nu, ~, nc] = generate_dimensions();
    Cx = randn(nc, nx);
    Cu = randn(nc, nu);
    cl = -10*abs(randn(nc, 1));
    cu = +10*abs(randn(nc, 1));
endfunction
