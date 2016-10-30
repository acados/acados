
1;

function w_star = solve_structured_ocp(N, nx, nu, A, B, b, x0, Q, S, R, q, r, xl, xu, ul, uu,
    Cx, Cu, cl, cu)

    global TOLERANCE;
    % Cost function: min 0.5 * w^T H w + h^T w
    H = blkdiag(kron(eye(N), [Q, S.'; S, R]), Q);
    h = [repmat([q; r], N, 1); q];

    % Equality constraint: G * w + g = 0
    G = [-eye(nx), zeros(nx, N*(nx+nu))];
    for i=1:N
        G = [G; zeros(nx, (i-1)*(nx+nu)), A, B, -eye(nx), zeros(nx, (N-i)*(nx+nu))];
    end
    g = [x0; repmat(b, N, 1)];

    % Simple bounds: lbw <= w <= ubw
    lbw = [repmat([xl; ul], N, 1); xl];
    ubw = [repmat([xu; uu], N, 1); xu];

    % linear inequality constraint: A_ineq * w <= b_ineq
    lin_ineq = blkdiag(kron(eye(N), [Cx, Cu]), Cx);
    A_ineq = [lin_ineq; -lin_ineq];
    b_ineq = [repmat(cu, N+1, 1); -repmat(cl, N+1, 1)];

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

    nc = (N+1)*length(cl);
    KKT_system = [];
    for i=0:N-1
        xk = w_star(i*(nx+nu)+1 : i*(nx+nu)+nx);
        uk = w_star(i*(nx+nu)+nx+1 : (i+1)*(nx+nu));
        lambdak = lambda_star(i*nx+1 : (i+1)*nx);
        lambdak_1 = lambda_star((i+1)*nx+1 : (i+2)*nx);
        mukx = mu_star(i*(nx+nu)+1 : i*(nx+nu)+nx);
        muku = mu_star(i*(nx+nu)+nx+1 : (i+1)*(nx+nu));
        nuk_ub = nu_star(i*(nx+nu)+1 : (i+1)*(nx+nu));
        nuk_lb = nu_star(nc+i*(nx+nu)+1 : nc+(i+1)*(nx+nu));

        KKT_system = [KKT_system; Q*xk + q + S.'*uk - lambdak + A.'*lambdak_1 + mukx + Cx.'*nuk_ub - Cx.'*nuk_lb];
        KKT_system = [KKT_system; R*uk + r + S*xk + B.'*lambdak_1 + muku + Cu.'*nuk_ub - Cu.'*nuk_lb];
    end
    lambdaN = lambda_star(N*(nx)+1 : end);
    muN = mu_star(N*(nx+nu)+1 : end);
    nuN_ub = nu_star(N*(nx+nu)+1 : (N+1)*(nx+nu));
    nuN_lb = nu_star(nc+N*(nx+nu)+1 : nc+(N+1)*(nx+nu));
    xN = w_star(N*(nx+nu)+1 : end);

    KKT_system = [KKT_system; Q*xN + q - lambdaN + muN + Cx.'*nuN_ub - Cx.'*nuN_lb];
    if(norm(KKT_system) > TOLERANCE)
        KKT_system
        error(['Solution is not optimal! norm of KKT: ', num2str(norm(KKT_system))]);
    end
endfunction

function [G, g, A_bar, B_bar] = calculate_transition_quantities(N, nx, nu, A, B, b, x0)
    A_bar = [-eye(nx),zeros(nx,(N-1)*nx)];
    for i=1:N-1
        A_bar = [A_bar; zeros(nx, (i-1)*(nx)), A, -eye(nx), zeros(nx, (N-i-1)*nx)];
    end
    B_bar = kron(eye(N),B);
    b_bar = repmat(b, N, 1);

    G = -A_bar \ B_bar;
    g = -A_bar \ (b_bar + [A*x0; zeros((N-1)*nx,1)]);
endfunction

function [H_bar, h_bar] = calculate_condensed_cost_function(N, nx, nu, Q, S, R, q, r, A, B, b, x0)
    global TOLERANCE;

    [G, g, A_bar, B_bar] = calculate_transition_quantities(N, nx, nu, A, B, b, x0);
    h_bar = repmat(r, N, 1) + G.' * (repmat(q, N, 1) + kron(eye(N), Q) * g) + kron(eye(N), S) * [x0; g(1:end-nx)];
    S_cut = [zeros(nu, N*nx); [kron(eye(N-1), S), zeros(nu, nx)]];
    H_bar = kron(eye(N), R) + G.'*kron(eye(N), Q)*G + S_cut*G + G.'*S_cut.';

    % As a check, do condensing based on Frasch2014a
    c = repmat(b, N, 1);
    L = [A; zeros((N-1)*nx, nx)];
    G2 = [zeros(nx, N*(nu)); -A_bar \ B_bar];
    g2 = [zeros(nx, 1); -A_bar \ c];
    Ge = [eye(nx); -A_bar \ L];
    Gex0 = Ge*x0;
    Q_all = kron(eye(N+1), Q);
    S_all = [kron(eye(N), S), zeros(N*nu, nx)];
    R_all = kron(eye(N), R);
    r_all = repmat(r, N, 1);
    q_all = repmat(q, N+1, 1);
    Sex0 = (G2.'*Q_all + S_all)*Gex0;
    R_condensed = R_all + G2.'*Q_all*G2 + S_all*G2 + G2.'*S_all.';
    r_condensed = r_all + G2.'*(q_all + Q_all*g2) + S_all*g2 + Sex0;
    if(norm(H_bar - R_condensed) > TOLERANCE || norm(h_bar - r_condensed) > TOLERANCE)
        error('Condensing methods do not match!')
    end
endfunction

function [u_lb, u_ub, G_lb, G_ub] = calculate_condensed_bounds(N, ul, uu, xl, xu, g)
    u_lb = repmat(ul, N, 1);
    u_ub = repmat(uu, N, 1);
    G_lb = repmat(xl, N, 1) - g;
    G_ub = repmat(xu, N, 1) - g;
endfunction

function [C_bar, c_bar_lb, c_bar_ub] = calculate_condensed_general_constraints(N, nx, nu, nc,
    Cx, Cu, cl, cu, x0, G, g, g_lb, g_ub)

    Dx = kron(eye(N), Cx);
    Du = [kron(eye(N), Cu); zeros(nc, N*nu)];
    Dx_bar = [zeros(nc, N*nu); Dx*G];

    D_bar = Du + Dx_bar;
    dl = repmat(cl, N+1, 1) - [Cx*x0; Dx*g];
    du = repmat(cu, N+1, 1) - [Cx*x0; Dx*g];

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
