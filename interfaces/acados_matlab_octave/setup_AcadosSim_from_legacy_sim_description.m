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

function sim = setup_AcadosSim_from_legacy_sim_description(model_old, opts_old)
    % create
    sim = AcadosSim();

    % aliases of frontend obj
    model = model_old.model_struct;
    opts = opts_old.opts_struct;

    % dimensions
    dims_fields_map = struct(...
        'dim_nx', 'nx', ...
        'dim_nu', 'nu', ...
        'dim_np', 'np', ...
        'dim_nz', 'nz', ...
        'dim_gnsf_nx1', 'gnsf_nx1', ...
        'dim_gnsf_nx2', 'gnsf_nx2', ...
        'dim_gnsf_nz1', 'gnsf_nz1', ...
        'dim_gnsf_nz2', 'gnsf_nz2', ...
        'dim_gnsf_ny', 'gnsf_ny', ...
        'dim_gnsf_nuhat', 'gnsf_nuhat', ...
        'dim_gnsf_nout', 'gnsf_nout' ...
    );
    fields = fieldnames(dims_fields_map);
    num_fields = length(fields);
    for n=1:num_fields
        old_field = fields{n};
        if isfield(model, old_field)
            new_field = dims_fields_map.(old_field);
            sim.dims.(new_field) = model.(old_field);
        end
    end

    if isfield(model, 'parameter_values')
        sim.parameter_values = opts.parameter_values(:);
    end

    % model dynamics
    sim.model.name = model.name;
    if strcmp(opts.method, 'erk')
        sim.model.f_expl_expr = model.dyn_expr_f;
    elseif strcmp(opts.method, 'irk') || strcmp(opts.method, 'irk_gnsf')
        if strcmp(model.ext_fun_type, 'casadi')
            sim.model.f_impl_expr = model.dyn_expr_f;
        elseif strcmp(model.ext_fun_type, 'generic')
            sim.model.dyn_generic_source = model.dyn_generic_source;
        end
    elseif strcmp(opts.method, 'discrete')
        sim.model.dyn_ext_fun_type = model.ext_fun_type;
        if strcmp(model.ext_fun_type, 'casadi')
            sim.model.f_phi_expr = model.dyn_expr_phi;
        elseif strcmp(model.ext_fun_type, 'generic')
            sim.model.dyn_generic_source = model.dyn_generic_source;
            if isfield(model, 'dyn_disc_fun_jac_hess')
                sim.model.dyn_disc_fun_jac_hess = model.dyn_disc_fun_jac_hess;
            end
            if isfield(model, 'dyn_disc_fun_jac')
                sim.model.dyn_disc_fun_jac = model.dyn_disc_fun_jac;
            end
            sim.model.dyn_disc_fun = model.dyn_disc_fun;
        end
    elseif strcmp(opts.method, 'irk_gnsf')
        gnsf_fields_map = struct(...
        'dyn_gnsf_A', 'dyn_gnsf_A', ...
        'dyn_gnsf_B', 'dyn_gnsf_B', ...
        'dyn_gnsf_C', 'dyn_gnsf_C', ...
        'dyn_gnsf_E', 'dyn_gnsf_E', ...
        'dyn_gnsf_c', 'dyn_gnsf_c', ...
        'dyn_gnsf_A_LO', 'dyn_gnsf_A_LO', ...
        'dyn_gnsf_B_LO', 'dyn_gnsf_B_LO', ...
        'dyn_gnsf_E_LO', 'dyn_gnsf_E_LO', ...
        'dyn_gnsf_c_LO', 'dyn_gnsf_c_LO', ...
        'dyn_gnsf_L_x', 'dyn_gnsf_L_x', ...
        'dyn_gnsf_L_u', 'dyn_gnsf_L_u', ...
        'dyn_gnsf_L_xdot', 'dyn_gnsf_L_xdot', ...
        'dyn_gnsf_L_z', 'dyn_gnsf_L_z', ...
        'dyn_gnsf_expr_phi', 'dyn_gnsf_expr_phi', ...
        'dyn_gnsf_expr_f_lo', 'dyn_gnsf_expr_f_lo', ...
        'dyn_gnsf_ipiv_x', 'dyn_gnsf_ipiv_x', ...
        'dyn_gnsf_idx_perm_x', 'dyn_gnsf_idx_perm_x', ...
        'dyn_gnsf_ipiv_z', 'dyn_gnsf_ipiv_z', ...
        'dyn_gnsf_idx_perm_z', 'dyn_gnsf_idx_perm_z', ...
        'dyn_gnsf_ipiv_f', 'dyn_gnsf_ipiv_f', ...
        'dyn_gnsf_idx_perm_f', 'dyn_gnsf_idx_perm_f', ...
        'dyn_gnsf_nontrivial_f_LO', 'dyn_gnsf_nontrivial_f_LO', ...
        'dyn_gnsf_purely_linear', 'dyn_gnsf_purely_linear', ...
        'sym_gnsf_y', 'sym_gnsf_y', ...
        'sym_gnsf_uhat', 'sym_gnsf_uhat' ...
        );
        fields = fieldnames(gnsf_fields_map);
        num_fields = length(fields);
        for n=1:num_fields
            old_field = fields{n};
            if isfield(model, old_field)
                new_field = gnsf_fields_map.(old_field);
                sim.model.(new_field) = model.(old_field);
            end
        end
    else
        error(['integrator ', opts.method, ' not support for templating backend.'])
    end
    model_fields_map = struct(...
        'sym_x', 'x', ...
        'sym_xdot', 'xdot', ...
        'sym_u', 'u', ...
        'sym_p', 'p', ...
        'sym_z', 'z', ...
        'sym_t', 't' ...
    );
    fields = fieldnames(model_fields_map);
    num_fields = length(fields);
    for n=1:num_fields
        old_field = fields{n};
        if isfield(model, old_field)
            new_field = model_fields_map.(old_field);
            sim.model.(new_field) = model.(old_field);
        end
    end

    % options
    if strcmp(upper(opts.method), 'IRK_GNSF')
        sim.solver_options.integrator_type = 'GNSF';
    else
        sim.solver_options.integrator_type = upper(opts.method);
    end
    sim.solver_options.collocation_type = upper(opts.collocation_type);
    sim.solver_options.num_stages = opts.num_stages;
    sim.solver_options.num_steps = opts.num_steps;
    sim.solver_options.newton_iter = opts.newton_iter;
    sim.solver_options.newton_tol = opts.newton_tol;
    sim.solver_options.Tsim = model.T;
    sim.solver_options.sens_forw = str2bool(opts.sens_forw);
    sim.solver_options.sens_adj = str2bool(opts.sens_adj);
    sim.solver_options.sens_algebraic = str2bool(opts.sens_algebraic);
    sim.solver_options.sens_hess = str2bool(opts.sens_hess);
    sim.solver_options.output_z = str2bool(opts.output_z);
    sim.solver_options.jac_reuse = str2bool(opts.jac_reuse);
    sim.solver_options.ext_fun_compile_flags = opts.ext_fun_compile_flags;
    sim.solver_options.ext_fun_expand_dyn = opts.ext_fun_expand_dyn;

end


function b = str2bool(s)
    if strcmp(lower(s), 'false')
        b = false;
    elseif strcmp(lower(s), 'true')
        b = true;
    else
        error('either true or false.');
    end
end
