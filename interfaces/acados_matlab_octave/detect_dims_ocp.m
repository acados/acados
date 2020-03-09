%
% Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
% Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
% Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
% Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
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

function model = detect_dims_ocp(model)

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
        nbx_0 = 0;
    end

    model.dim_nbx_0 = nbx_0;

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

end
