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

% Condensing routine that outputs data against which acados is tested.
clear

% Check if package 'optim' is installed, used for quadprog
list = pkg('list', 'optim');
if(isempty(list))
    try
        pkg install -forge struct -verbose
        pkg install -forge optim  -verbose
    catch
        % Try install without downloading
        pkg install /home/travis/octave_download/struct-1.0.14.tar.gz
        pkg install /home/travis/octave_download/optim-1.5.2.tar.gz
    end
end
pkg load optim

% Tolerance used to determine optimality, to compare matrices, etc.
global TOLERANCE = 1e-10;

for mode = {'LTI', 'LTV'}

    prefix = mode{1};
    eval([prefix, '_generation_functions']);
    eval('condensing_functions');

    [N, nx, nu, nb, nc] = generate_dimensions();
    [A, B, b, x0] = generate_dynamics();
    [Q, S, R, q, r] = generate_cost_function();
    [xl, xu, ul, uu] = generate_bounds(x0);
    [Cx, Cu, cl, cu] = generate_general_constraints();

    mkdir(prefix);

    saveIntVector(N,'N',prefix);
    saveIntVector(repmat(nx,1,N+1),'nx',prefix);
    saveIntVector(repmat(nu,1,N),'nu',prefix);
    saveIntVector([repmat(nx+nu,1,N) nx],'nb',prefix);
    saveIntVector(repmat(nc,1,N+1),'nc',prefix);
    saveDoubleVector(x0,'x0',prefix);

    for ii = 1:N+1
    
      if ii <= N
        Rloc = R(:,(ii-1)*nu+1:ii*nu);  
        Sloc = S(:,(ii-1)*nx+1:ii*nx);
        rloc = r(:,ii);
        Aloc = A(:,(ii-1)*nx+1:ii*nx);
        Bloc = B(:,(ii-1)*nu+1:ii*nu);
        bloc = b(:,ii);
        Culoc = Cu(:,(ii-1)*nu+1:ii*nu);
        ubloc = [xu(:,ii); uu(:,ii)];
        lbloc = [xl(:,ii); ul(:,ii)];
        idxb = (1:nx+nu)-1;
        
        saveDoubleMatrix(Aloc,['A' num2str(ii-1)],prefix);
        saveDoubleMatrix(Bloc,['B' num2str(ii-1)],prefix);
        saveDoubleVector(bloc,['bv' num2str(ii-1)],prefix);
        saveDoubleMatrix(Rloc,['R' num2str(ii-1)],prefix);
        saveDoubleMatrix(Sloc,['S' num2str(ii-1)],prefix);
        saveDoubleVector(rloc,['rv' num2str(ii-1)],prefix);
        saveDoubleMatrix(Culoc,['Cu' num2str(ii-1)],prefix);  
      else
        ubloc = xu(:,ii);
        lbloc = xl(:,ii);
        idxb = (1:nx)-1;
      end
      
      Qloc = Q(:,(ii-1)*nx+1:ii*nx);
      qloc = q(:,ii);
      Cxloc = Cx(:,(ii-1)*nx+1:ii*nx);
      ucloc = cu(:,ii);
      lcloc = cl(:,ii);
            
      saveDoubleMatrix(Qloc,['Q' num2str(ii-1)],prefix);
      saveDoubleVector(qloc,['qv' num2str(ii-1)],prefix);
      saveDoubleMatrix(Cxloc,['Cx' num2str(ii-1)],prefix);
      saveDoubleVector(ucloc,['uc' num2str(ii-1)],prefix);
      saveDoubleVector(lcloc,['lc' num2str(ii-1)],prefix);
      saveDoubleVector(ubloc,['ub' num2str(ii-1)],prefix);
      saveDoubleVector(lbloc,['lb' num2str(ii-1)],prefix);
      saveIntVector(idxb,['idxb' num2str(ii-1)],prefix);

    end
    
    cd(prefix);
    save('N.dat', 'N', '-ascii', '-double');
    save('nx.dat', 'nx', '-ascii', '-double');
    save('nu.dat', 'nu', '-ascii', '-double');
    save('nb.dat', 'nb', '-ascii', '-double')
    save('nc.dat', 'nc', '-ascii', '-double')

    save('x0.dat', 'x0', '-ascii', '-double');
    save('A.dat', 'A', '-ascii', '-double');
    save('B.dat', 'B', '-ascii', '-double');
    save('bv.dat', 'b', '-ascii', '-double');

    save('Q.dat', 'Q', '-ascii', '-double');
    save('S.dat', 'S', '-ascii', '-double');
    save('R.dat', 'R', '-ascii', '-double');
    save('qv.dat', 'q', '-ascii', '-double');
    save('rv.dat', 'r', '-ascii', '-double');

    save('upper_bound_x.dat', 'xu', '-ascii', '-double');
    save('lower_bound_x.dat', 'xl', '-ascii', '-double');
    save('upper_bound_u.dat', 'uu', '-ascii', '-double');
    save('lower_bound_u.dat', 'ul', '-ascii', '-double');

    save('general_constraint_x.dat', 'Cx', '-ascii', '-double');
    save('general_constraint_u.dat', 'Cu', '-ascii', '-double');
    save('general_constraint_ub.dat', 'cu', '-ascii', '-double');
    save('general_constraint_lb.dat', 'cl', '-ascii', '-double');

    [w_star_ocp,pi_star_ocp] = solve_structured_ocp(N, nx, nu, nc, A, B, b, x0, Q, S, R, q, r, xl, xu, ul, uu,
        Cx, Cu, cl, cu);
     
    % TODO(dimitris-robin): update solve_structured_ocp function to allow for empty bounds/constraints and eliminate functions below 
    w_star_ocp_unconstrained = solve_structured_ocp_unconstrained(N, nx, nu, A, B, b, x0, Q, S, R, q, r);
    w_star_ocp_bounds = solve_structured_ocp_bounds(N, nx, nu, A, B, b, x0, Q, S, R, q, r, xl, xu, ul, uu);
    w_star_ocp_no_bounds = solve_structured_ocp_no_bounds(N, nx, nu, nc, A, B, b, x0, Q, S, R, q, r, Cx, Cu, cl, cu);   
   
    % Do condensing
    [G, g, A_bar, B_bar] = calculate_transition_quantities(N, nx, nu, A, B, b, x0);
    [H_bar, h_bar] = calculate_condensed_cost_function(N, nx, nu, Q, S, R, q, r, A, B, b, x0);
    [u_lb, u_ub, g_lb, g_ub] = calculate_condensed_bounds(N, nx, ul, uu, xl, xu, g);
    [C_bar, c_bar_lb, c_bar_ub] = calculate_condensed_general_constraints(N, nx, nu, nc,
        Cx, Cu, cl, cu, x0, G, g, g_lb, g_ub);

    % Solve condensed QP
    [w_star_condensed_quadprog, ~, exit_flag, ~, condensed_multipliers] = quadprog(H_bar, h_bar,
        [C_bar; -C_bar], [c_bar_ub; -c_bar_lb], [], [], u_lb, u_ub);
    if(~(exit_flag == 1))
        error(['Condensed QP solution failed with code: ', num2str(exit_flag)]);
    end

    % Compare sparse and condensed solution
    XU = reshape([w_star_ocp;zeros(nu,1)], nx+nu, N+1);
    x_star_ocp = XU(1:nx, :);
    x_star_ocp = x_star_ocp(:);
    u_star_ocp = XU(nx+1:end, 1:end-1);
    u_star_ocp = u_star_ocp(:);
    if (norm(u_star_ocp - w_star_condensed_quadprog) > TOLERANCE)
        [u_star_ocp w_star_condensed_quadprog]
        error(['Difference between condensed and sparse solution: ',
            num2str(norm(u_star_ocp - w_star_condensed_quadprog))])
    end

    % Calculate unconstrained solution
    [w_star_condensed_quadprog_unconstrained, ~, exit_flag] = quadprog(H_bar, h_bar);
    if(~(exit_flag == 1))
        error(['Unconstrained condensed QP solution failed with code: ', num2str(exit_flag)]);
    end
     
    % Save data to file
    save('transition_vector.dat', 'g', '-ascii', '-double');
    save('transition_matrix.dat', 'G', '-ascii', '-double');
    save('condensed_hessian.dat', 'H_bar', '-ascii', '-double');
    save('condensed_gradient.dat', 'h_bar', '-ascii', '-double');
    % TODO(dimitris-robin): why do we save both this and the lower_bound_xu?
    save('u_lower_bound.dat', 'u_lb', '-ascii', '-double');
    save('u_upper_bound.dat', 'u_ub', '-ascii', '-double');
    save('condensed_lower_bound.dat', 'g_lb', '-ascii', '-double');
    save('condensed_upper_bound.dat', 'g_ub', '-ascii', '-double');
    save('condensed_general_constraint_matrix.dat', 'C_bar', '-ascii', '-double');
    save('condensed_general_constraint_lb.dat', 'c_bar_lb', '-ascii', '-double');
    save('condensed_general_constraint_ub.dat', 'c_bar_ub', '-ascii', '-double');
    save('w_star_ocp_constrained.dat', 'w_star_ocp', '-ascii', '-double');
    save('pi_star_ocp_constrained.dat', 'pi_star_ocp', '-ascii', '-double');
    save('w_star_ocp_bounds.dat', 'w_star_ocp_bounds', '-ascii', '-double');
    save('w_star_ocp_no_bounds.dat', 'w_star_ocp_no_bounds', '-ascii', '-double'); % only with affine constraints
    save('w_star_ocp_unconstrained.dat', 'w_star_ocp_unconstrained', '-ascii', '-double');
    cd('..');
    
    saveDoubleVector(w_star_ocp,'sol_constrained',prefix);
    saveDoubleVector(w_star_ocp_bounds,'sol_only_bounds',prefix);
    saveDoubleVector(w_star_ocp_no_bounds,'sol_only_ineq',prefix);
    saveDoubleVector(w_star_ocp_unconstrained,'sol_only_x0',prefix);

    
end
