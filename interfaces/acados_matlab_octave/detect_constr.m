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

function model = detect_constr(model, is_e, is_init)

    if ~exist('is_init','var')
        is_init = false;
    end
   

    import casadi.*

    x = model.sym_x;
    u = model.sym_u;

    % check type
    if strcmp(class(x(1)), 'casadi.SX')
        isSX = true;
    else
        disp('constraint detection only works for SX CasADi type!!!');
        keyboard;
    end

    if isfield(model, 'sym_z')
        z = model.sym_z;
    else
        z = SX.sym('z', 0, 0);
    end

    nx = length(x);
    nu = length(u);
    nz = length(z);

    % z = model.sym_z;
    if is_e
        expr_constr = model.constr_expr_h_e;
        LB = model.constr_lh_e;
        UB = model.constr_uh_e;
        fprintf('\nConstraint detection for terminal constraints.\n');
    elseif is_init
        expr_constr = model.constr_expr_h_0;
        LB = model.constr_lh_0;
        UB = model.constr_uh_0;
        fprintf('\nConstraint detection for initial constraints.\n');
    else
        expr_constr = model.constr_expr_h;
        LB = model.constr_lh;
        UB = model.constr_uh;
        fprintf('\nConstraint detection for path constraints.\n');
    end
    constr_fun = Function('constr_fun', {x, u, z}, {expr_constr});

    % initialize
    constr_expr_h = SX.sym('constr_expr_h', 0, 0);
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
        if any(c.which_depends(z)) || ~c.is_linear([ x; u ]) || any(c.which_depends(model.sym_p))
            % external constraint
            constr_expr_h = vertcat(constr_expr_h, c);
            lh = [ lh; LB(ii)];
            uh = [ uh; UB(ii)];
            disp(['constraint ', num2str(ii), ' is kept as nonlinear constraint.']);
            disp(c);
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
                    disp(['constraint ', num2str(ii),...
                          ' is reformulated as bound on x.']);
                    disp(c);
                    disp(' ')

                else
                    % bound on u;
                    Jbu = [Jbu; zeros(1,nu)];
                    Jbu(end, idb-nx) = 1;
                    lbu = [lbu; LB(ii)/Jc(idb)];
                    ubu = [ubu; UB(ii)/Jc(idb)];
                    disp(['constraint ', num2str(ii),...
                          ' is reformulated as bound on u.']);
                    disp(c);
                    disp(' ')

                end
            else
                % c is general linear constraint
                C = [C; Jc(1:nx)];
                D = [D; Jc(nx+1:end)];
                lg = [ lg; LB(ii)];
                ug = [ ug; UB(ii)];
                disp(['constraint ', num2str(ii),...
                      ' is reformulated as general linear constraint.']);
                disp(c);
                disp(' ')
                
            end
        end
    end

    
    if is_e
        % checks
        if any(expr_constr.which_depends(u)) || ~isempty(lbu) || (~isempty(D) && any(D))
            error('terminal constraint may not depend on control input.');
        end
        % h
        model.constr_type_e = 'bgh';
        model.dim_nh_e = length(lh);
        if ~isempty(lh)
            model.constr_expr_h_e = constr_expr_h;
            model.constr_lh_e = lh;
            model.constr_uh_e = uh;
        else
            model = rmfield(model, 'constr_expr_h_e');
            model = rmfield(model, 'constr_lh_e');
            model = rmfield(model, 'constr_uh_e');
        end
        % g
        if ~isempty(lg)
            model.constr_C_e = C;
            model.constr_lg_e = lg;
            model.constr_ug_e = ug;
            model.dim_ng_e = length(lg);
        end
        % bounds x
        if ~isempty(lbx)
            model.constr_Jbx_e = Jbx;
            model.constr_lbx_e = lbx;
            model.constr_ubx_e = ubx;
            model.dim_nbx_e = length(lbx);
        end

    elseif is_init
        model.constr_type = 'bgh';
        % h
        model.dim_nh_0 = length(constr_lh_0);
        if ~isempty(lh)
            model.constr_expr_h_0 = constr_expr_h;
            model.constr_lh_0 = lh;
            model.constr_uh_0 = uh;
        else
            model = rmfield(model, 'constr_expr_h_0');
            model = rmfield(model, 'constr_lh_0');
            model = rmfield(model, 'constr_uh_0');
        end
        % g
        model.dim_ng = length(lg);
        if ~isempty(lg)
            model.constr_C = C;
            model.constr_D = D;
            model.constr_lg = lg;
            model.constr_ug = ug;
        end
        % bounds x
        if ~isempty(lbx)
            model.constr_Jbx = Jbx;
            model.constr_lbx = lbx;
            model.constr_ubx = ubx;
            model.dim_nbx = length(lbx);
        end
        % bounds u
        if ~isempty(lbu)
            model.constr_Jbu = Jbu;
            model.constr_lbu = lbu;
            model.constr_ubu = ubu;
            model.dim_nbu = length(lbu);
        end
        
    else
        model.constr_type = 'bgh';
        % h
        model.dim_nh = length(lh);
        if ~isempty(lh)
            model.constr_expr_h = constr_expr_h;
            model.constr_lh = lh;
            model.constr_uh = uh;
        else
            model = rmfield(model, 'constr_expr_h');
            model = rmfield(model, 'constr_lh');
            model = rmfield(model, 'constr_uh');
        end
        % g
        model.dim_ng = length(lg);
        if ~isempty(lg)
            model.constr_C = C;
            model.constr_D = D;
            model.constr_lg = lg;
            model.constr_ug = ug;
        end
        % bounds x
        if ~isempty(lbx)
            model.constr_Jbx = Jbx;
            model.constr_lbx = lbx;
            model.constr_ubx = ubx;
            model.dim_nbx = length(lbx);
        end
        % bounds u
        if ~isempty(lbu)
            model.constr_Jbu = Jbu;
            model.constr_lbu = lbu;
            model.constr_ubu = ubu;
            model.dim_nbu = length(lbu);
        end
        
    end
%     model

end
