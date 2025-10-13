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

% NOT SUPPORTED FOR NOW..
%  slack constraints
% nsbu;  // number of softened input bounds
% nsbx;  // number of softened state bounds
% nsg;  // number of softened general linear constraints
% nsh;  // number of softened nonlinear constraints

function model = detect_constraint_structure(model, constraints, stage_type)

    import casadi.*

    x = model.x;
    u = model.u;
    z = model.z;

    nx = length(x);
    nu = length(u);

    if isa(x, 'casadi.SX')
        isSX = true;
    else
        error('Constraint detection only works for casadi.SX!');
    end

    if strcmp(stage_type, 'initial')
        expr_constr = model.con_h_expr_0;
        LB = constraints.lh_0;
        UB = constraints.uh_0;
        fprintf('\nConstraint detection for initial constraints.\n');
    elseif strcmp(stage_type, 'path')
        expr_constr = model.con_h_expr;
        LB = constraints.lh;
        UB = constraints.uh;
        fprintf('\nConstraint detection for path constraints.\n');
    elseif strcmp(stage_type, 'terminal')
        expr_constr = model.con_h_expr_e;
        LB = constraints.lh_e;
        UB = constraints.uh_e;
        fprintf('\nConstraint detection for terminal constraints.\n');
    else
        error('Constraint detection: Wrong stage_type.')
    end

    if isempty(expr_constr)
        expr_constr = SX.sym('con_h_expr', 0, 0);
    end

    if ~(isa(expr_constr, 'casadi.SX') || isa(expr_constr, 'casadi.MX'))
        disp('expr_constr =')
        disp(expr_constr)
        error("Constraint type detection requires definition of constraints as CasADi SX or MX.")
    end

    % initialize
    constr_expr_h = SX.sym('con_h_expr', 0, 0);
    lh = [];
    uh = [];

    C = zeros(0, nx);
    D = zeros(0, nu);
    lg = [];
    ug = [];

    Jbx = zeros(0, nx);
    lbx = [];
    ubx = [];

    Jbu = [];
    lbu = [];
    ubu = [];

    % loop over CasADi formulated constraints
    for ii = 1:length(expr_constr)
        c = expr_constr(ii);
        if any(c.which_depends(z)) || ~c.is_linear([ x; u ]) || any(c.which_depends(model.p)) || any(c.which_depends(model.p_global))
            % external constraint
            constr_expr_h = vertcat(constr_expr_h, c);
            lh = [ lh; LB(ii)];
            uh = [ uh; UB(ii)];
            disp(['Constraint ', num2str(ii), ' is kept as a nonlinear constraint.']);
            disp('Constraint expression: ');
            disp(c)
            disp(' ')
        else % c is linear in x and u
            Jc_fun = Function('Jc_fun', {x(1)}, {jacobian(c, [x;u])});
            Jc = full(Jc_fun(0));

            if length( nonzeros(Jc) ) == 1
                % c is bound
                idb = find(Jc);
                if idb <= nx
                    % bound on x
                    Jbx = [Jbx; zeros(1, nx)];
                    Jbx(end, idb) = 1;
                    lbx = [lbx; LB(ii)/Jc(idb)];
                    ubx = [ubx; UB(ii)/Jc(idb)];
                    disp(['Constraint ', num2str(ii),...
                          ' is reformulated as a bound on x.']);
                    disp('Constraint expression: ');
                    disp(c)
                    disp(' ')
                else
                    % bound on u;
                    Jbu = [Jbu; zeros(1,nu)];
                    Jbu(end, idb-nx) = 1;
                    lbu = [lbu; LB(ii)/Jc(idb)];
                    ubu = [ubu; UB(ii)/Jc(idb)];
                    disp(['Constraint ', num2str(ii),...
                          ' is reformulated as a bound on u.']);
                    disp('Constraint expression: ');
                    disp(c)
                    disp(' ')
                end
            else
                % c is general linear constraint
                C = [C; Jc(1:nx)];
                D = [D; Jc(nx+1:end)];
                lg = [ lg; LB(ii)];
                ug = [ ug; UB(ii)];
                disp(['Constraint ', num2str(ii),...
                      ' is reformulated as a general linear constraint.']);
                disp('Constraint expression: ');
                disp(c)
                disp(' ')
            end
        end
    end


    if strcmp(stage_type, 'terminal')
        % checks
        if any(expr_constr.which_depends(u)) || ~isempty(lbu) || (~isempty(D) && any(D))
            error('Terminal constraint may not depend on control input.');
        end
        % h
        constraints.constr_type_e = 'BGH';
        if ~isempty(lh)
            model.con_h_expr_e = constr_expr_h;
            constraints.lh_e = lh;
            constraints.uh_e = uh;
        else
            model.con_h_expr_e = [];
            constraints.lh_e = [];
            constraints.uh_e = [];
        end
        % g
        if ~isempty(lg)
            constraints.C_e = C;
            constraints.lg_e = lg;
            constraints.ug_e = ug;
        end
        % bounds x
        if ~isempty(lbx)
            constraints.idxbx_e = J_to_idx(Jbx);
            constraints.lbx_e = lbx;
            constraints.ubx_e = ubx;
        end

    elseif strcmp(stage_type, 'initial')
        warning("At initial stage, only h constraints are detected.")
        constraints.constr_type_0 = 'BGH';
        % h
        if ~isempty(lh)
            model.con_h_expr_0 = constr_expr_h;
            constraints.lh_0 = lh;
            constraints.uh_0 = uh;
        else
            model.con_h_expr_0 = [];
            constraints.lh_0 = [];
            constraints.uh_0 = [];
        end
    else % path
        constraints.constr_type = 'BGH';
        % h
        if ~isempty(lh)
            model.con_h_expr = constr_expr_h;
            constraints.lh = lh;
            constraints.uh = uh;
        else
            model.con_h_expr = [];
            constraints.lh = [];
            constraints.uh = [];
        end
        % g
        if ~isempty(lg)
            constraints.C = C;
            constraints.D = D;
            constraints.lg = lg;
            constraints.ug = ug;
        end
        % bounds x
        if ~isempty(lbx)
            constraints.idxbx = J_to_idx(Jbx);
            constraints.lbx = lbx;
            constraints.ubx = ubx;
        end
        % bounds u
        if ~isempty(lbu)
            constraints.idxbu = J_to_idx(Jbu);
            constraints.lbu = lbu;
            constraints.ubu = ubu;
        end
    end
end

% TODO directly detect idxbx, etc. instead of Jbx, etc.

function idx = J_to_idx(J)
    size_J = size(J);
    nrows = size_J(1);
    idx = zeros(nrows,1);
    for i = 1:nrows
        this_idx = find(J(i,:));
        if length(this_idx) ~= 1
            error(['J_to_idx: Invalid J matrix. Exiting. Found more than one nonzero in row ' num2str(i)]);
        end
        if J(i,this_idx) ~= 1
            error(['J_to_idx: J matrices can only contain 1s, got J(' num2str(i) ', ' num2str(this_idx) ') = ' num2str(J(i,this_idx)) ]);
        end
        idx(i) = this_idx - 1; % store 0-based index
    end
end
