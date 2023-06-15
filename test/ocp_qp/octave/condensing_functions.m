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


1;

function [H, h] = ocp_cost_function(N, nx, nu, Q, S, R, q, r)
    % cost function: min 0.5 * w^T H w + h^T w
    H = [];
    h = [];
    for k=1:N
        Qk = Q(:, (k-1)*nx+1 : k*nx);
        Sk = S(:, (k-1)*nx+1 : k*nx);
        Rk = R(:, (k-1)*nu+1 : k*nu);
        H = blkdiag(H, [Qk, Sk.'; Sk, Rk]);
        h = [h; q(:, k); r(:, k)];
    end
    H = blkdiag(H, Q(:, end-nx+1 : end));
    h = [h; q(:, N+1)];
endfunction

function [G, g] = ocp_equality_constraints(N, nx, nu, A, B, b, x0)
    % equality constraint: G * w + g = 0
    G = [-eye(nx), zeros(nx, N*(nx+nu))];
    g = x0;
    for k=1:N
        Ak = A(:, (k-1)*nx+1 : k*nx);
        Bk = B(:, (k-1)*nu+1 : k*nu);
        G = [G; zeros(nx, (k-1)*(nx+nu)), Ak, Bk, -eye(nx), zeros(nx, (N-k)*(nx+nu))];
        g = [g; b(:, k)];
    end
endfunction

function [lbw, ubw] = ocp_simple_bounds(N, nx, nu, xl, ul, xu, uu)
    % simple bounds: lbw <= w <= ubw
    lbxu = [xl; ul, zeros(nu, 1)];
    lbw = lbxu(1 : end-nu).';
    ubxu = [xu; uu, zeros(nu, 1)];
    ubw = ubxu(1 : end-nu).';
endfunction

function [A_ineq, b_ineq] = ocp_general_constraints(N, nx, nu, nc, Cx, Cu, cl, cu)
    % linear inequality constraint: A_ineq * w <= b_ineq
    lin_ineq = [];
    for k=1:N
        Cxk = Cx(:, (k-1)*nx+1 : k*nx);
        Cuk = Cu(:, (k-1)*nu+1 : k*nu);
        lin_ineq = blkdiag(lin_ineq, [Cxk, Cuk]);
    end
    lin_ineq = blkdiag(lin_ineq, Cx(:, N*nx+1 : (N+1)*nx));
    A_ineq = [lin_ineq; -lin_ineq];
    b_ineq = [cu(:); -cl(:)];
endfunction

function [w_star,pi_star] = solve_structured_ocp(N, nx, nu, nc, A, B, b, x0, Q, S, R, q, r, xl, xu, ul, uu,
    Cx, Cu, cl, cu)

    global TOLERANCE;

    [H, h] = ocp_cost_function(N, nx, nu, Q, S, R, q, r);
    [G, g] = ocp_equality_constraints(N, nx, nu, A, B, b, x0);
    [lbw, ubw] = ocp_simple_bounds(N, nx, nu, xl, ul, xu, uu);
    [A_ineq, b_ineq] = ocp_general_constraints(N, nx, nu, nc, Cx, Cu, cl, cu);

    % Solve OCP
    [w_star, ~, exit_flag, ~, all_multipliers] = quadprog(H, h, A_ineq, b_ineq, G, -g, lbw, ubw);

    if(~(exit_flag == 1))
        Z = null(G);
        disp(['convex QP? : ', num2str(all(eig(Z.'*H*Z) > 1e-4))])
        error(['QP solution failed with code: ', num2str(exit_flag)]);
    end

    % Check consistency of solution with KKT system
    lambda_star = all_multipliers.eqlin;
    mu_star = -all_multipliers.lower + all_multipliers.upper;
    nu_star = all_multipliers.ineqlin;

    pi_star = lambda_star(nx+1:end);

    nc_all = (N+1)*nc;
    KKT_system = [];
    for k=1:N
        xk = w_star((k-1)*(nx+nu)+1 : (k-1)*(nx+nu)+nx);
        uk = w_star((k-1)*(nx+nu)+nx+1 : k*(nx+nu));
        lambdak = lambda_star((k-1)*nx+1 : k*nx);
        lambdak_1 = lambda_star(k*nx+1 : (k+1)*nx);
        mukx = mu_star((k-1)*(nx+nu)+1 : (k-1)*(nx+nu)+nx);
        muku = mu_star((k-1)*(nx+nu)+nx+1 : k*(nx+nu));
        nuk_ub = nu_star((k-1)*(nx+nu)+1 : k*(nx+nu));
        nuk_lb = nu_star(nc_all+(k-1)*(nx+nu)+1 : nc_all+k*(nx+nu));
        Ak = A(:, (k-1)*nx+1 : k*nx);
        Bk = B(:, (k-1)*nu+1 : k*nu);
        Qk = Q(:, (k-1)*nx+1 : k*nx);
        Sk = S(:, (k-1)*nx+1 : k*nx);
        Rk = R(:, (k-1)*nu+1 : k*nu);
        qk = q(:, k);
        rk = r(:, k);
        Cxk = Cx(:, (k-1)*nx+1 : k*nx);
        Cuk = Cu(:, (k-1)*nu+1 : k*nu);
        KKT_system = [KKT_system; Qk*xk + qk + Sk.'*uk - lambdak + Ak.'*lambdak_1 + mukx + Cxk.'*nuk_ub - Cxk.'*nuk_lb];
        KKT_system = [KKT_system; Rk*uk + rk + Sk*xk + Bk.'*lambdak_1 + muku + Cuk.'*nuk_ub - Cuk.'*nuk_lb];
    end
    lambdaN = lambda_star(N*(nx)+1 : end);
    muN = mu_star(N*(nx+nu)+1 : end);
    nuN_ub = nu_star(N*(nx+nu)+1 : (N+1)*(nx+nu));
    nuN_lb = nu_star(nc_all+N*(nx+nu)+1 : nc_all+(N+1)*(nx+nu));
    xN = w_star(N*(nx+nu)+1 : end);
    QN = Q(:, end-nx+1 : end);
    CxN = Cx(:, N*nx+1 : end);
    KKT_system = [KKT_system; QN*xN + q(:, N+1) - lambdaN + muN + CxN.'*nuN_ub - CxN.'*nuN_lb];
    if(norm(KKT_system) > TOLERANCE)
        KKT_system
        error(['Solution is not optimal! norm of KKT: ', num2str(norm(KKT_system))]);
    end
endfunction

% TODO(dimitris-robin): remove once function above is extended to cover empty constraints
function [w_star] = solve_structured_ocp_bounds(N, nx, nu, A, B, b, x0, Q, S, R, q, r, xl, xu, ul, uu)

    [H, h] = ocp_cost_function(N, nx, nu, Q, S, R, q, r);
    [G, g] = ocp_equality_constraints(N, nx, nu, A, B, b, x0);
    [lbw, ubw] = ocp_simple_bounds(N, nx, nu, xl, ul, xu, uu);

    % Solve OCP
    [w_star, ~, exit_flag] = quadprog(H, h, [], [], G, -g, lbw, ubw);
    
endfunction

function [w_star] = solve_structured_ocp_no_bounds(N, nx, nu, nc, A, B, b, x0, Q, S, R, q, r, Cx, Cu, cl, cu)

    [H, h] = ocp_cost_function(N, nx, nu, Q, S, R, q, r);
    [G, g] = ocp_equality_constraints(N, nx, nu, A, B, b, x0);
    [A_ineq, b_ineq] = ocp_general_constraints(N, nx, nu, nc, Cx, Cu, cl, cu);

    % Solve OCP
    [w_star, ~, exit_flag] = quadprog(H, h, A_ineq, b_ineq, G, -g, [], []);
endfunction

function [w_star] = solve_structured_ocp_unconstrained(N, nx, nu, A, B, b, x0, Q, S, R, q, r)

    [H, h] = ocp_cost_function(N, nx, nu, Q, S, R, q, r);
    [G, g] = ocp_equality_constraints(N, nx, nu, A, B, b, x0);

    % Solve OCP
    [w_star, ~, exit_flag] = quadprog(H, h, [], [], G, -g, [], []);

 endfunction

function [G, g, A_bar, B_bar] = calculate_transition_quantities(N, nx, nu, A, B, b, x0)
    A_bar = [-eye(nx),zeros(nx, (N-1)*nx)];
    B_bar = B(:, 1 : nu);
    for k=1:N-1
        Ak = A(:, k*nx+1 : (k+1)*nx);
        Bk = B(:, k*nu+1 : (k+1)*nu);
        A_bar = [A_bar; zeros(nx, (k-1)*(nx)), Ak, -eye(nx), zeros(nx, (N-k-1)*nx)];
        B_bar = blkdiag(B_bar, Bk);
    end
    b_bar = b(:);

    G = -A_bar \ B_bar;
    A0 = A(:, 1 : nx);
    g = -A_bar \ (b_bar + [A0*x0; zeros((N-1)*nx,1)]);
endfunction

function [H_bar, h_bar] = calculate_condensed_cost_function(N, nx, nu, Q, S, R, q, r, A, B, b, x0)
    global TOLERANCE;

    [G, g, A_bar, B_bar] = calculate_transition_quantities(N, nx, nu, A, B, b, x0);
    Q_bar = [];
    S_bar = [];
    R_bar = [];
    q_bar = [];
    r_bar = [];
    for k=1:N
        Q_bar = blkdiag(Q_bar, Q(:, (k-1)*nx+1 : k*nx));
        S_bar = blkdiag(S_bar, S(:, (k-1)*nx+1 : k*nx));
        R_bar = blkdiag(R_bar, R(:, (k-1)*nu+1 : k*nu));
        q_bar = [q_bar; q(:, k)];
        r_bar = [r_bar; r(:, k)];
    end
    Q_bar = blkdiag(Q_bar, Q(:, N*nx+1 : end));
    q_bar = [q_bar; q(:, N+1)];

    Q_bar_cut = Q_bar(nx+1 : end, nx+1 : end);
    q_bar_cut = q_bar(nx+1 : end);
    h_bar = r_bar + G.' * (q_bar_cut + Q_bar_cut * g) + S_bar * [x0; g(1:end-nx)];
    S_cut = [zeros(nu, N*nx); [S_bar(nu+1 : end, nx+1 : end), zeros((N-1)*nu, nx)]];
    H_bar = R_bar + G.'*Q_bar_cut*G + S_cut*G + G.'*S_cut.';

    % As a check, do condensing based on Frasch2014a
    c = b(:);
    A0 = A(:, 1 : nx);
    L = [A0; zeros((N-1)*nx, nx)];
    G2 = [zeros(nx, N*(nu)); -A_bar \ B_bar];
    g2 = [zeros(nx, 1); -A_bar \ c];
    Ge = [eye(nx); -A_bar \ L];
    Gex0 = Ge*x0;
    S_all = [S_bar, zeros(N*nu, nx)];
    Sex0 = (G2.'*Q_bar + S_all)*Gex0;
    R_condensed = R_bar + G2.'*Q_bar*G2 + S_all*G2 + G2.'*S_all.';
    r_condensed = r_bar + G2.'*(q_bar + Q_bar*g2) + S_all*g2 + Sex0;
    if(norm(H_bar - R_condensed) > TOLERANCE || norm(h_bar - r_condensed) > TOLERANCE)
        error(['Condensing methods do not match! Max difference = ', num2str(
            max(norm(H_bar - R_condensed), norm(h_bar - r_condensed))
        )])
    end
endfunction

function [u_lb, u_ub, g_lb, g_ub] = calculate_condensed_bounds(N, nx, ul, uu, xl, xu, g)
    u_lb = ul(:);
    u_ub = uu(:);
    g_lb = xl(nx+1 : end).' - g;
    g_ub = xu(nx+1 : end).' - g;
endfunction

function [C_bar, c_bar_lb, c_bar_ub] = calculate_condensed_general_constraints(N, nx, nu, nc,
    Cx, Cu, cl, cu, x0, G, g, g_lb, g_ub)

    Dx = [];
    Du = [];
    for k=1:N
        Dx = blkdiag(Dx, Cx(:, k*nx+1 : (k+1)*nx));
        Du = blkdiag(Du, Cu(:, (k-1)*nu+1 : k*nu));
    end
    Du = [Du; zeros(nc, N*nu)];
    Dx_bar = [zeros(nc, N*nu); Dx*G];

    D_bar = Du + Dx_bar;
    Cx0 = Cx(:, 1 : nx);
    dl = cl(:) - [Cx0*x0; Dx*g];
    du = cu(:) - [Cx0*x0; Dx*g];

    C_bar = D_bar(1 : nc, :);
    c_bar_lb = dl(1 : nc);
    c_bar_ub = du(1 : nc);
    for k=1:N
        Gk = G((k-1)*nx+1 : k*nx, :);
        D_bark = D_bar(k*nc+1 : (k+1)*nc, :);
        C_bar = [C_bar; Gk; D_bark];
        c_bar_lb = [c_bar_lb; g_lb((k-1)*nx+1 : k*nx); dl(k*nc+1 : (k+1)*nc)];
        c_bar_ub = [c_bar_ub; g_ub((k-1)*nx+1 : k*nx); du(k*nc+1 : (k+1)*nc)];
    end
endfunction
