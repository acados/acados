%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% The 2-Clause BSD License
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.;

%


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
    [N, nx, nu] = generate_dimensions();
    do
        A = randn(nx, nx);
    until(all(abs(eig(A)) < 0.9)) % Unstable dynamics leads to unstable condensing
    A = repmat(A, 1, N);
    B = randn(nx, nu);
    B = repmat(B, 1, N);
    b = randn(nx, 1);
    b = repmat(b, 1, N);
    x0 = randn(nx, 1);
endfunction

function [Q, S, R, q, r] = generate_cost_function()
    [N, nx, nu] = generate_dimensions();
    do
        Q = randn(nx, nx);
        Q = Q.'*Q;
        S = randn(nu, nx);
        R = randn(nu, nu);
        R = R.'*R;
    until(all(eig([Q, S.'; S, R]) > 1e-2)) % Implies a convex QP
    Q = repmat(Q, 1, N+1);
    S = repmat(S, 1, N);
    R = repmat(R, 1, N);
    q = randn(nx, 1);
    q = repmat(q, 1, N+1);
    r = randn(nu, 1);
    r = repmat(r, 1, N);
endfunction

function [xl, xu, ul, uu] = generate_bounds(x0)
    [N, nx, nu, nb] = generate_dimensions();
    xl = -5*(abs(x0) + abs(randn(nx,1)));
    xl = repmat(xl, 1, N+1);
    xu = +5*(abs(x0) + abs(randn(nx,1)));
    xu = repmat(xu, 1, N+1);
    ul = -10*abs(randn(nu,1));
    ul = repmat(ul, 1, N);
    uu = +10*abs(randn(nu,1));
    uu = repmat(uu, 1, N);
endfunction

function [Cx, Cu, cl, cu] = generate_general_constraints()
    [N, nx, nu, ~, nc] = generate_dimensions();
    Cx = randn(nc, nx);
    Cx = repmat(Cx, 1, N+1);
    Cu = randn(nc, nu);
    Cu = repmat(Cu, 1, N);
    cl = -10*abs(randn(nc, 1));
    cl = repmat(cl, 1, N+1);
    cu = +10*abs(randn(nc, 1));
    cu = repmat(cu, 1, N+1);
endfunction
