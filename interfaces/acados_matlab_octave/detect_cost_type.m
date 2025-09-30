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
%   Author: Jonathan Frey: jonathanpaulfrey(at)gmail.com

function detect_cost_type(model, cost, dims, stage_type)

    import casadi.*

    x = model.x;
    u = model.u;
    z = model.z;
    p = model.p;

    nx = length(x);
    nu = length(u);
    nz = length(z);

    disp('--------------------------------------------------------------');
    if strcmp(stage_type, 'terminal')
        expr_cost = model.cost_expr_ext_cost_e;
        disp('Structure detection for terminal cost term');
    elseif strcmp(stage_type, 'path')
        expr_cost = model.cost_expr_ext_cost;
        disp('Structure detection for path cost');
    elseif strcmp(stage_type, 'initial')
        expr_cost = model.cost_expr_ext_cost_0;
        disp('Structure detection for initial cost term');
    end

    if ~(isa(expr_cost, 'casadi.SX') || isa(expr_cost, 'casadi.MX'))
        disp('expr_cost =')
        disp(expr_cost)
        error("Cost type detection require definition of cost term as CasADi SX or MX.")
    end


    if expr_cost.is_quadratic(x) && expr_cost.is_quadratic(u) && expr_cost.is_quadratic(z) ...
            && ~any(expr_cost.which_depends(p)) && ~any(expr_cost.which_depends(model.p_global)) ...
            && ~any(expr_cost.which_depends(model.t))

        if expr_cost.is_zero()
            fprintf('Cost function is zero -> Reformulating as LINEAR_LS cost.\n');
            cost_type = 'LINEAR_LS';
            ny = 0;
            Vx = []; Vu = []; Vz = []; W = []; y_ref = []; y = [];
        else
            cost_fun = Function('cost_fun', {x, u, z}, {expr_cost});
            dummy = SX.sym('dummy', 1, 1);

            fprintf('Cost function is quadratic -> Reformulating as LINEAR_LS cost.\n');

            Hxuz_fun = Function('Hxuz_fun', {dummy}, {hessian(expr_cost, [x; u; z])});
            H_xuz = full(Hxuz_fun(0));

            xuz_idx = [];
            for i = 1:(nx+nu+nz)
                if ~isempty(find(H_xuz(i,:), 1) )
                    xuz_idx = union(xuz_idx, i);
                end
            end
            x_idx = intersect(1:nx, xuz_idx);
            u_idx = intersect(1+nx:nx+nu, xuz_idx);
            z_idx = intersect(1+nx+nu : nx+nu+nz, xuz_idx);

            ny = length(xuz_idx);

            Vx = zeros(ny, nx);
            Vu = zeros(ny, nu);
            Vz = zeros(ny, nz);
            W = zeros(ny);

            i = 1;
            for id = x_idx
                Vx(i, id) = 1;
                W(i, :) = H_xuz(id, xuz_idx)/2;
                i = i+1;
            end

            for id = u_idx
                iu = id - nx;
                Vu(i, iu) = 1;
                W(i, :) = H_xuz(id, xuz_idx)/2;
                i = i+1;
            end

            for id = z_idx
                iz = id - nx - nu;
                Vz(i, iz) = 1;
                W(i, :) = H_xuz(id, xuz_idx)/2;
                i = i+1;
            end

            xuz = [x; u; z];
            y = xuz(xuz_idx);
            jac_fun = Function('jac_fun', {y}, {jacobian(expr_cost, y)'});
            y_ref = -W \ ( .5 * full(jac_fun(zeros(ny,1))) );

            y = -y_ref + Vx * x + Vu * u;
            if nz > 0
                y = y + Vz * z;
            end
            lls_cost_fun = Function('lls_cost_fun', {x, u, z}, {y' * W * y});

            rel_err_tol = 1e-13;
            for jj = 1:5
                x0 = rand(nx,1);
                u0 = rand(nu,1);
                z0 = rand(nz,1);

                val1 = full(lls_cost_fun(x0, u0, z0));
                val2 = full(cost_fun(x0, u0, z0));
                diff_eval = abs(val1-val2);
                rel_error = diff_eval / max(abs(val1), abs(val2));
                if rel_error > rel_err_tol
                    disp(['something went wrong when reformulating with linear least square cost',...
                    ' got relative error ', num2str(rel_error, '%e'), ' should be < ', num2str(rel_err_tol, '%e')]);
                    keyboard
                end
            end

            %% take into account 1/2 factor in linear least square module
            W = 2 * W;
        end

        %% extract output
        if strcmp(stage_type, 'terminal')
            if ~isempty(find(Vu,1))
                error('Cost mayer term cannot depend on control input u!');
            end
            if ~isempty(find(Vz,1))
                error('Cost mayer term cannot depend on z!');
            end
            cost.cost_type_e = 'LINEAR_LS';
            dims.ny_e = ny;
            cost.Vx_e = Vx;
            cost.W_e = W;
            cost.yref_e = y_ref;
        elseif strcmp(stage_type, 'path')
            cost.cost_type = 'LINEAR_LS';
            dims.ny = ny;
            cost.Vx = Vx;
            cost.Vu = Vu;
            cost.Vz = Vz;
            cost.W = W;
            cost.yref = y_ref;
        elseif strcmp(stage_type, 'initial')
            cost.cost_type_0 = 'LINEAR_LS';
            dims.ny_0 = ny;
            cost.Vx_0 = Vx;
            cost.Vu_0 = Vu;
            cost.Vz_0 = Vz;
            cost.W_0 = W;
            cost.yref_0 = y_ref;
        end
        fprintf('\n\nreformulated cost term in linear least squares form with:')
        fprintf('\ncost = 0.5 * || Vx * x + Vu * u + Vz * z - y_ref ||_W\n');
        fprintf('\nVx\n');
        disp(Vx);
        fprintf('\nVu\n');
        disp(Vu);
        fprintf('\nVz\n');
        disp(Vz);
        fprintf('\nW\n');
        disp(W);
        fprintf('\ny_ref\n');
        disp(y_ref);
        fprintf('\ny (symbolic)\n');
        disp(y);
        fprintf('\nNOTE: These numerical values can be updated online using the appropriate setters.\n');
% elseif
    %  TODO: can nonLINEAR_LS be detected?!
    else
        fprintf('\n\nCost function is not quadratic or depends on parameters -> Using external cost\n\n');
        if strcmp(stage_type, 'terminal')
            cost.cost_type_e = 'EXTERNAL';
        elseif strcmp(stage_type, 'path')
            cost.cost_type = 'EXTERNAL';
        elseif strcmp(stage_type, 'initial')
            cost.cost_type_0 = 'EXTERNAL';
        end
    end
    disp('--------------------------------------------------------------');

end
