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

function sim_json = set_up_acados_sim_json(obj)
    % create
    sim_json = acados_template_mex.acados_sim_json();

    % aliases of frontend obj
    model = obj.model_struct;
    opts = obj.opts_struct;

    % file structures
    acados_folder = getenv('ACADOS_INSTALL_DIR');
    sim_json.acados_include_path = [acados_folder, '/include'];
    sim_json.acados_lib_path = [acados_folder, '/lib'];
    sim_json.shared_lib_ext = '.so';
    if ismac()
        sim_json.shared_lib_ext = '.dylib';
    end
    sim_json.json_file = 'acados_sim.json';
    sim_json.cython_include_dirs = [];  % only useful for python interface
    sim_json.code_export_directory = fullfile(pwd(), 'c_generated_code');

    % dimensions
    sim_json.dims.nx = model.dim_nx;
    sim_json.dims.nu = model.dim_nu;
    sim_json.dims.nz = model.dim_nz;
    sim_json.dims.np = model.dim_np;
    if isfield(model, 'dim_gnsf_nx1')
        sim_json.dims.gnsf_nx1 = model.dim_gnsf_nx1;
        sim_json.dims.gnsf_nz1 = model.dim_gnsf_nz1;
        sim_json.dims.gnsf_nout = model.dim_gnsf_nout;
        sim_json.dims.gnsf_ny = model.dim_gnsf_ny;
        sim_json.dims.gnsf_nuhat = model.dim_gnsf_nuhat;
    end

    % model dynamics
    sim_json.model.name = model.name;
    if strcmp(opts.method, 'erk')
        sim_json.model.f_expl_expr = model.dyn_expr_f;
    elseif strcmp(opts.method, 'irk')
        if strcmp(model.ext_fun_type, 'casadi')
            sim_json.model.f_impl_expr = model.dyn_expr_f;
        elseif strcmp(model.ext_fun_type, 'generic')
            sim_json.model.dyn_generic_source = model.dyn_generic_source;
        end
    elseif strcmp(opts.method, 'discrete')
        sim_json.model.dyn_ext_fun_type = model.ext_fun_type;
        if strcmp(model.ext_fun_type, 'casadi')
            sim_json.model.f_phi_expr = model.dyn_expr_phi;
        elseif strcmp(model.ext_fun_type, 'generic')
            sim_json.model.dyn_generic_source = model.dyn_generic_source;
            if isfield(model, 'dyn_disc_fun_jac_hess')
                sim_json.model.dyn_disc_fun_jac_hess = model.dyn_disc_fun_jac_hess;
            end
            if isfield(model, 'dyn_disc_fun_jac')
                sim_json.model.dyn_disc_fun_jac = model.dyn_disc_fun_jac;
            end
            sim_json.model.dyn_disc_fun = model.dyn_disc_fun;
        end
    elseif strcmp(opts.method, 'irk_gnsf')
        sim_json.model.gnsf.A = model.dyn_gnsf_A;
        sim_json.model.gnsf.B = model.dyn_gnsf_B;
        sim_json.model.gnsf.C = model.dyn_gnsf_C;
        sim_json.model.gnsf.E = model.dyn_gnsf_E;
        sim_json.model.gnsf.c = model.dyn_gnsf_c;
        sim_json.model.gnsf.A_LO = model.dyn_gnsf_A_LO;
        sim_json.model.gnsf.B_LO = model.dyn_gnsf_B_LO;
        sim_json.model.gnsf.E_LO = model.dyn_gnsf_E_LO;
        sim_json.model.gnsf.c_LO = model.dyn_gnsf_c_LO;

        sim_json.model.gnsf.L_x = model.dyn_gnsf_L_x;
        sim_json.model.gnsf.L_u = model.dyn_gnsf_L_u;
        sim_json.model.gnsf.L_xdot = model.dyn_gnsf_L_xdot;
        sim_json.model.gnsf.L_z = model.dyn_gnsf_L_z;

        sim_json.model.gnsf.expr_phi = model.dyn_gnsf_expr_phi;
        sim_json.model.gnsf.expr_f_lo = model.dyn_gnsf_expr_f_lo;

        sim_json.model.gnsf.ipiv_x = model.dyn_gnsf_ipiv_x;
        sim_json.model.gnsf.idx_perm_x = model.dyn_gnsf_idx_perm_x;
        sim_json.model.gnsf.ipiv_z = model.dyn_gnsf_ipiv_z;
        sim_json.model.gnsf.idx_perm_z = model.dyn_gnsf_idx_perm_z;
        sim_json.model.gnsf.ipiv_f = model.dyn_gnsf_ipiv_f;
        sim_json.model.gnsf.idx_perm_f = model.dyn_gnsf_idx_perm_f;

        sim_json.model.gnsf.nontrivial_f_LO = model.dyn_gnsf_nontrivial_f_LO;
        sim_json.model.gnsf.purely_linear = model.dyn_gnsf_purely_linear;

        sim_json.model.gnsf.y = model.sym_gnsf_y;
        sim_json.model.gnsf.uhat = model.sym_gnsf_uhat;
    else
        error(['integrator ', opts.method, ' not support for templating backend.'])
    end

    sim_json.model.x = model.sym_x;
    sim_json.model.u = model.sym_u;
    sim_json.model.z = model.sym_z;
    sim_json.model.xdot = model.sym_xdot;
    sim_json.model.p = model.sym_p;

    % options
    if strcmp(upper(opts.method), 'IRK_GNSF')
        sim_json.sim_options.integrator_type = 'GNSF';
    else
        sim_json.sim_options.integrator_type = upper(opts.method);
    end
    sim_json.sim_options.collocation_type = upper(opts.collocation_type);
    sim_json.sim_options.sim_method_num_stages = opts.num_stages;
    sim_json.sim_options.sim_method_num_steps = opts.num_steps;
    sim_json.sim_options.sim_method_newton_iter = opts.newton_iter;
    sim_json.sim_options.sim_method_newton_tol = opts.newton_tol;
    sim_json.sim_options.Tsim = model.T;
    sim_json.sim_options.sens_forw = str2bool(opts.sens_forw);
    sim_json.sim_options.sens_adj = str2bool(opts.sens_adj);
    sim_json.sim_options.sens_algebraic = str2bool(opts.sens_algebraic);
    sim_json.sim_options.sens_hess = str2bool(opts.sens_hess);
    sim_json.sim_options.output_z = str2bool(opts.output_z);
    sim_json.sim_options.sim_method_jac_reuse = str2bool(opts.jac_reuse);
    sim_json.sim_options.ext_fun_compile_flags = opts.ext_fun_compile_flags;

    % parameters
    if model.dim_np > 0
        if isempty(opts.parameter_values)
            warning(['opts_struct.parameter_values are not set.', ...
                        10 'Using zeros(np,1) by default.' 10 'You can update them later using set().']);
            sim_json.parameter_values = zeros(model.dim_np,1);
        else
            sim_json.parameter_values = opts.parameter_values(:);
        end
    end
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
