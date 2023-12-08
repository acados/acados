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

function [model, opts] = detect_dims_ocp(model, opts)

    %% general
    model.dim_nx = length(model.sym_x);

    if isfield(model, 'sym_u')
        model.dim_nu = length(model.sym_u);
    else
        model.dim_nu = 0;
    end

    if isfield(model, 'sym_z')
        model.dim_nz = length(model.sym_z);
    else
        model.dim_nz = 0;
    end

    if isfield(model, 'sym_p')
        model.dim_np = length(model.sym_p);
    else
        model.dim_np = 0;
    end

    %% cost
    % initial
    if strcmp( model.cost_type_0, 'linear_ls')
        if isfield(model, 'cost_W_0') && isfield(model, 'cost_Vx_0') && isfield(model, 'cost_Vu_0')
            ny = length(model.cost_W_0);
            if ny ~= size(model.cost_Vx_0, 1) || ny ~= size(model.cost_Vu_0, 1)
                error('inconsistent dimension ny, regarding W, Vx, Vu.');
            end
        else
            error('setting linear least square cost: need W, Vx, Vu, at least one missing.')
        end
        model.dim_ny_0 = ny;
    elseif strcmp( model.cost_type_0, 'nonlinear_ls')
        if isfield(model, 'cost_W_0') && isfield(model, 'cost_expr_y_0')
            ny = length(model.cost_W_0);
            if ny ~= length(model.cost_expr_y_0)
                error('inconsistent dimension ny, regarding W, expr_y.');
            end
        else
            error('setting nonlinear least square cost: need W_0, cost_expr_y_0, at least one missing.')
        end
        model.dim_ny_0 = ny;
    end

    % path
    if strcmp( model.cost_type, 'linear_ls')
        if isfield(model, 'cost_W') && isfield(model, 'cost_Vx') && isfield(model, 'cost_Vu')
            ny = length(model.cost_W);
            if ny ~= size(model.cost_Vx, 1) || ny ~= size(model.cost_Vu, 1)
                error('inconsistent dimension ny, regarding W, Vx, Vu.');
            end
        else
            error('setting linear least square cost: need W, Vx, Vu, at least one missing.')
        end
        model.dim_ny = ny;
    elseif strcmp( model.cost_type, 'nonlinear_ls')
        if isfield(model, 'cost_W') && isfield(model, 'cost_expr_y')
            ny = length(model.cost_W);
            if ny ~= length(model.cost_expr_y)
                error('inconsistent dimension ny, regarding W, expr_y.');
            end
        else
            error('setting nonlinear least square cost: need W, cost_expr_y, at least one missing.')
        end
        model.dim_ny = ny;
    end

    % terminal
    if strcmp( model.cost_type_e, 'linear_ls')
        if isfield(model, 'cost_W_e') && isfield(model, 'cost_Vx_e')
            ny_e = length(model.cost_W_e);
            if ny_e ~= size(model.cost_Vx_e, 1)
                error('inconsistent dimension ny_e, regarding W_e, Vx_e.');
            end
        elseif ~isfield(model, 'cost_W_e') && ~isfield(model, 'cost_Vx_e')
            ny_e = 0;
            warning('Fields cost_W_e and cost_Vx_e not provided. Using empty ls terminal cost.')
        else
            error('setting linear least square cost: need W_e, Vx_e, at least one missing.')
        end
        model.dim_ny_e = ny_e;
    elseif strcmp( model.cost_type_e, 'nonlinear_ls')
        if isfield(model, 'cost_W_e') && isfield(model, 'cost_expr_y_e')
            ny_e = length(model.cost_W_e);
            if ny_e ~= length(model.cost_expr_y_e)
                error('inconsistent dimension ny_e, regarding W_e, expr_y_e.');
            end
        else
            error('setting nonlinear least square cost: need W_e, cost_expr_y_e, at least one missing.')
        end
        model.dim_ny_e = ny_e;
    end


    %% constraints
    % initial
    if isfield(model, 'constr_Jbx_0') && isfield(model, 'constr_lbx_0') && isfield(model, 'constr_ubx_0')
        nbx_0 = length(model.constr_lbx_0);
        if nbx_0 ~= length(model.constr_ubx_0) || nbx_0 ~= size(model.constr_Jbx_0, 1)
            error('inconsistent dimension nbx_0, regarding Jbx_0, lbx_0, ubx_0.');
        end
    elseif isfield(model, 'constr_Jbx_0') || isfield(model, 'constr_lbx_0') || isfield(model, 'constr_ubx_0')
        error('setting bounds on x: need Jbx_0, lbx_0, ubx_0, at least one missing.');
    else
        % no initial state constraint
        disp("detect_dims_ocp: OCP without constraints on initial state detected.")
        nbx_0 = 0;
    end
    model.dim_nbx_0 = nbx_0;

    if isfield(model, 'constr_idxbxe_0')
        model.dim_nbxe_0 = length(model.constr_idxbxe_0);
    else
        % no equalities on initial state.
        model.constr_idxbxe_0 = [];
        model.dim_nbxe_0 = 0;
    end

    % path
    if isfield(model, 'constr_Jbx') && isfield(model, 'constr_lbx') && isfield(model, 'constr_ubx')
        nbx = length(model.constr_lbx);
        if nbx ~= length(model.constr_ubx) || nbx ~= size(model.constr_Jbx, 1)
            error('inconsistent dimension nbx, regarding Jbx, lbx, ubx.');
        end
    elseif isfield(model, 'constr_Jbx') || isfield(model, 'constr_lbx') || isfield(model, 'constr_ubx')
        error('setting bounds on x: need Jbx, lbx, ubx, at least one missing.');
    else
        nbx = 0;
    end
    model.dim_nbx = nbx;

    if isfield(model, 'constr_Jbu') && isfield(model, 'constr_lbu') && isfield(model, 'constr_ubu')
        nbu = length(model.constr_lbu);
        if nbu ~= length(model.constr_ubu) || nbu ~= size(model.constr_Jbu, 1)
            error('inconsistent dimension nbu, regarding Jbu, lbu, ubu.');
        end
    elseif isfield(model, 'constr_Jbu') || isfield(model, 'constr_lbu') || isfield(model, 'constr_ubu')
        error('setting bounds on u: need Jbu, lbu, ubu, at least one missing.');
    else
        nbu = 0;
    end
    model.dim_nbu = nbu;

    if isfield(model, 'constr_C') && isfield(model, 'constr_D') && ...
       isfield(model, 'constr_lg') && isfield(model, 'constr_ug')
        ng = length(model.constr_lg);
        if ng ~= length(model.constr_ug) || ng ~= size(model.constr_C, 1) || ng ~= size(model.constr_D, 1)
            error('inconsistent dimension ng, regarding C, D, lg, ug.');
        end
    elseif isfield(model, 'constr_C') || isfield(model, 'constr_D') || ...
           isfield(model, 'constr_lg') || isfield(model, 'constr_ug')
        error('setting general linear constraints: need C, D, lg, ug, at least one missing.');
    else
        ng = 0;
    end
    model.dim_ng = ng;

    if isfield(model, 'constr_expr_h') && ...
             isfield(model, 'constr_lh') && isfield(model, 'constr_uh')
        nh = length(model.constr_lh);
        if nh ~= length(model.constr_uh) || nh ~= length(model.constr_expr_h)
            error('inconsistent dimension nh, regarding expr_h, lh, uh.');
        end
    elseif isfield(model, 'constr_expr_h') || ...
           isfield(model, 'constr_lh') || isfield(model, 'constr_uh')
        error('setting external constraint function h: need expr_h, lh, uh at least one missing.');
    else
        nh = 0;
    end
    model.dim_nh = nh;

    if isfield(model, 'constr_expr_h_0') && ...
            isfield(model, 'constr_lh_0') && isfield(model, 'constr_uh_0')
    nh_0 = length(model.constr_lh_0);
    if nh_0 ~= length(model.constr_uh_0) || nh_0 ~= length(model.constr_expr_h_0)
        error('inconsistent dimension nh_0, regarding expr_h_0, lh_0, uh_0.');
    end
    elseif isfield(model, 'constr_expr_h_0') || ...
        isfield(model, 'constr_lh_0') || isfield(model, 'constr_uh_0')
    error('setting external constraint function h: need expr_h_0, lh_0, uh_0 at least one missing.');
    else
        nh_0 = 0;
    end
    model.dim_nh_0 = nh_0;

    % terminal
    if isfield(model, 'constr_Jbx_e') && isfield(model, 'constr_lbx_e') && isfield(model, 'constr_ubx_e')
        nbx_e = length(model.constr_lbx_e);
        if nbx_e ~= length(model.constr_ubx_e) || nbx_e ~= size(model.constr_Jbx_e, 1)
            error('inconsistent dimension nbx_e, regarding Jbx_e, lbx_e, ubx_e.');
        end
    elseif isfield(model, 'constr_Jbx_e') || isfield(model, 'constr_lbx_e') || isfield(model, 'constr_ubx_e')
        error('setting bounds on x: need Jbx_e, lbx_e, ubx_e, at least one missing.');
    else
        nbx_e = 0;
    end
    model.dim_nbx_e = nbx_e;

    if isfield(model, 'constr_C_e') && ...
       isfield(model, 'constr_lg_e') && isfield(model, 'constr_ug_e')
        ng_e = length(model.constr_lg_e);
        if ng_e ~= length(model.constr_ug_e) || ng_e ~= size(model.constr_C_e, 1)
            error('inconsistent dimension ng_e, regarding C_e, lg_e, ug_e.');
        end
    elseif isfield(model, 'constr_C_e') || ...
           isfield(model, 'constr_lg_e') || isfield(model, 'constr_ug_e')
        error('setting general linear constraints: need C_e, lg_e, ug_e, at least one missing.');
    else
        ng_e = 0;
    end
    model.dim_ng_e = ng_e;

    if isfield(model, 'constr_expr_h_e') && ...
             isfield(model, 'constr_lh_e') && isfield(model, 'constr_uh_e')
        nh_e = length(model.constr_lh_e);
        if nh_e ~= length(model.constr_uh_e) || nh_e ~= length(model.constr_expr_h_e)
            error('inconsistent dimension nh_e, regarding expr_h_e, lh_e, uh_e.');
        end
    elseif isfield(model, 'constr_expr_h_e') || ...
           isfield(model, 'constr_lh_e') || isfield(model, 'constr_uh_e')
        error('setting external constraint function h: need expr_h_e, lh_e, uh_e at least one missing.');
    else
        nh_e = 0;
    end
    model.dim_nh_e = nh_e;

    %% slack dimensions
    if isfield(model, 'constr_Jsbx')
        nsbx = size(model.constr_Jsbx, 2);
    else
        nsbx = 0;
    end

    if isfield(model, 'constr_Jsbu')
        nsbu = size(model.constr_Jsbu, 2);
    else
        nsbu = 0;
    end

    if isfield(model, 'constr_Jsg')
        nsg = size(model.constr_Jsg, 2);
    else
        nsg = 0;
    end
    if isfield(model, 'constr_Jsh')
        nsh = size(model.constr_Jsh, 2);
    else
        nsh = 0;
    end
    if isfield(model, 'constr_Jsphi')
        nsphi = size(model.constr_Jsphi, 2);
    else
        nsphi = 0;
    end

    ns = nsbx + nsbu + nsg + nsh + nsphi;
    wrong_field = '';
    if isfield(model, 'cost_Zl') && ~all(size(model.cost_Zl) == [ns, ns])
        wrong_field = 'Zl';
        dim = size(model.cost_Zl);
    elseif isfield(model, 'cost_Zu') && ~all(size(model.cost_Zu) == [ns, ns])
        wrong_field = 'Zu';
        dim = size(model.cost_Zu);
    elseif isfield(model, 'cost_zl') && ~all(size(model.cost_zl) == [ns, 1])
        wrong_field = 'zl';
        dim = size(model.cost_zl);
    elseif isfield(model, 'cost_zu') && ~all(size(model.cost_zu) == [ns, 1])
        wrong_field = 'zu';
        dim = size(model.cost_zu);
    end

    if ~strcmp(wrong_field, '')
        error(['Inconsistent size for field ', wrong_field, ' with dimension ', num2str(dim),...
              '. Detected ns = ', num2str(ns), ' = nsbx + nsbu + nsg + nsh + nsphi.',...
              ' With nsbx = ', num2str(nsbx), ', nsbu = ', num2str(nsbu), ' nsg = ', num2str(nsg),...
              ' nsh = ', num2str(nsh), ', nsphi = ', num2str(nsphi), '.'])
    end

    model.dim_ns = ns;
    model.dim_nsbx = nsbx;
    model.dim_nsbu = nsbu;
    model.dim_nsg = nsg;
    model.dim_nsh = nsh;
    model.dim_nsphi = nsphi;

    %% terminal slack dimensions
    if isfield(model, 'constr_Jsbx_e')
        nsbx_e = size(model.constr_Jsbx_e, 2);
    else
        nsbx_e = 0;
    end

    if isfield(model, 'constr_Jsg_e')
        nsg_e = size(model.constr_Jsg_e, 2);
    else
        nsg_e = 0;
    end
    if isfield(model, 'constr_Jsh_e')
        nsh_e = size(model.constr_Jsh_e, 2);
    else
        nsh_e = 0;
    end
    if isfield(model, 'constr_Jsphi_e')
        nsphi_e = size(model.constr_Jsphi_e, 2);
    else
        nsphi_e = 0;
    end

    ns_e = nsbx_e + nsg_e + nsh_e + nsphi_e;
    wrong_field = '';
    if isfield(model, 'cost_Zl_e') && ~all(size(model.cost_Zl_e) == [ns_e, ns_e])
        wrong_field = 'Zl_e';
        dim = size(model.cost_Zl_e);
    elseif isfield(model, 'cost_Zu_e') && ~all(size(model.cost_Zu_e) == [ns_e, ns_e])
        wrong_field = 'Zu_e';
        dim = size(model.cost_Zu_e);
    elseif isfield(model, 'cost_zl_e') && ~all(size(model.cost_zl_e) == [ns_e, 1])
        wrong_field = 'zl_e';
        dim = size(model.cost_zl_e);
    elseif isfield(model, 'cost_zu_e') && ~all(size(model.cost_zu_e) == [ns_e, 1])
        wrong_field = 'zu_e';
        dim = size(model.cost_zu_e);
    end

    if ~strcmp(wrong_field, '')
        error(['Inconsistent size for field', wrong_field, ' with dimension ', num2str(dim),...
                '. Detected ns_e = ', num2str(ns_e), ' = nsbx_e + nsg_e + nsh_e + nsphi_e.',...
                ' With nsbx_e = ', num2str(nsbx_e), ' nsg_e = ', num2str(nsg_e),...
                ' nsh_e = ', num2str(nsh_e), ', nsphi_e = ', num2str(nsphi_e), '.'])
    end

    model.dim_ns_e = ns_e;
    model.dim_nsbx_e = nsbx_e;
    model.dim_nsg_e = nsg_e;
    model.dim_nsh_e = nsh_e;
    model.dim_nsphi_e = nsphi_e;

    % shooting nodes -> time_steps
    N = opts.param_scheme_N;
    if ~isempty(opts.shooting_nodes)
        if N + 1 ~= length(opts.shooting_nodes)
            error('inconsistent dimension N regarding shooting nodes.');
        end
        for i=1:N
            opts.time_steps(i) = opts.shooting_nodes(i+1) - opts.shooting_nodes(i);
        end
        sum_time_steps = sum(opts.time_steps);
        if abs((sum_time_steps - model.T) / model.T) > 1e-14
            warning('shooting nodes are not consistent with time horizon T, rescaling automatically');
            opts.time_steps = opts.time_steps * model.T / sum_time_steps;
        end
    elseif ~isempty(opts.time_steps)
        if N ~= length(opts.time_steps)
            error('inconsistent dimension N regarding time steps.');
        end
        sum_time_steps = sum(opts.time_steps);
        if abs((sum_time_steps - model.T) / model.T) > 1e-14
            error(['time steps are not consistent with time horizon T, ', ...
                'got T = ' num2str(model.T) '; sum(time_steps) = ' num2str(sum_time_steps) '.']);
        end
        % just to have them available, e.g. for plotting;
        opts.shooting_nodes = zeros(N+1, 1);
        for i = 1:N
            opts.shooting_nodes(i+1) = sum(opts.time_steps(1:i));
        end
    else
        opts.time_steps = model.T/N * ones(N,1);
    end
    if any(opts.time_steps < 0)
        error(['ocp discretization: time_steps between shooting nodes must all be > 0', ...
            ' got: ' num2str(opts.time_steps)])
    end
    if ~isempty(opts.sim_method_num_stages)
        if(strcmp(opts.sim_method, "erk"))
            if (any(opts.sim_method_num_stages < 1) || any(opts.sim_method_num_stages > 4))
                error(['ERK: num_stages = ', num2str(opts.sim_method_num_stages) ' not available. Only number of stages = {1,2,3,4} implemented!']);
            end
        end
    end

    % qp_dunes
    if ~isempty(strfind(opts.qp_solver,'qpdunes'))
        model.constr_idxbxe_0 = [];
        model.dim_nbxe_0 = 0;
    end

end
