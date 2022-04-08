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

function ocp_generate_casadi_ext_fun(model_struct, opts_struct)

model_name = model_struct.name;

% get acados folder
acados_folder = getenv('ACADOS_INSTALL_DIR');

% select files to compile
c_files = {};
% dynamics
if (strcmp(model_struct.dyn_type, 'explicit'))
    % generate c for function and derivatives using casadi
    if (strcmp(opts_struct.codgen_model, 'true'))
        generate_c_code_explicit_ode(model_struct, opts_struct);
    end
    % sources list
    c_files{end+1} = [model_name, '_dyn_expl_ode_fun.c'];
    c_files{end+1} = [model_name, '_dyn_expl_vde_forw.c'];
    c_files{end+1} = [model_name, '_dyn_expl_vde_adj.c'];
    c_files{end+1} = [model_name, '_dyn_expl_ode_hess.c'];
elseif (strcmp(model_struct.dyn_type, 'implicit'))
    if (strcmp(opts_struct.sim_method, 'irk'))
        % generate c for function and derivatives using casadi
        if (strcmp(opts_struct.codgen_model, 'true'))
            generate_c_code_implicit_ode(model_struct, opts_struct);
        end
        % sources list
        c_files{end+1} = [model_name, '_dyn_impl_dae_fun.c'];
        c_files{end+1} = [model_name, '_dyn_impl_dae_fun_jac_x_xdot_z.c'];
        c_files{end+1} = [model_name, '_dyn_impl_dae_fun_jac_x_xdot_u.c'];
        c_files{end+1} = [model_name, '_dyn_impl_dae_jac_x_xdot_u_z.c'];
        if strcmp(opts_struct.nlp_solver_exact_hessian, 'true')
            c_files{end+1} = [model_name, '_dyn_impl_dae_hess.c'];
        end
    elseif (strcmp(opts_struct.sim_method, 'irk_gnsf'))
        % generate c for function and derivatives using casadi
        if (strcmp(opts_struct.codgen_model, 'true'))
            generate_c_code_gnsf(model_struct, opts_struct);
        end
        % sources list
        c_files{end+1} = [model_name, '_dyn_gnsf_get_matrices_fun.c'];
        if ~model_struct.dyn_gnsf_purely_linear
            if model_struct.dyn_gnsf_nontrivial_f_LO
                c_files{end+1} = [model_name, '_dyn_gnsf_f_lo_fun_jac_x1k1uz.c'];
            end
            c_files{end+1} = [model_name, '_dyn_gnsf_phi_fun.c'];
            c_files{end+1} = [model_name, '_dyn_gnsf_phi_fun_jac_y.c'];
            c_files{end+1} = [model_name, '_dyn_gnsf_phi_jac_y_uhat.c'];
        end
    else
        fprintf('\nocp_generate_casadi_ext_fun: sim_method not supported: %s\n', opts_struct.sim_method);
        return;
    end
elseif (strcmp(model_struct.dyn_type, 'discrete'))
    if (strcmp(model_struct.dyn_ext_fun_type, 'casadi'))
        % generate c for function and derivatives using casadi
        if (strcmp(opts_struct.codgen_model, 'true'))
            generate_c_code_disc_dyn(model_struct, opts_struct);
        end
        % sources list
        c_files{end+1} = [model_name, '_dyn_disc_phi_fun.c'];
        c_files{end+1} = [model_name, '_dyn_disc_phi_fun_jac.c'];
        c_files{end+1} = [model_name, '_dyn_disc_phi_fun_jac_hess.c'];
    end
else
    fprintf('\ncodegen_model: dyn_type not supported: %s\n', model_struct.dyn_type);
    return;
end
% nonlinear constraints
if (strcmp(model_struct.constr_type, 'bgh') && (isfield(model_struct, 'constr_expr_h') || isfield(model_struct, 'constr_expr_h_e')))
    % generate c for function and derivatives using casadi
    if (strcmp(opts_struct.codgen_model, 'true'))
        generate_c_code_nonlinear_constr(model_struct, opts_struct);
    end
    % sources list
    if isfield(model_struct, 'constr_expr_h')
        c_files{end+1} = [model_name, '_constr_h_fun.c'];
        c_files{end+1} = [model_name, '_constr_h_fun_jac_uxt_zt.c'];
        c_files{end+1} = [model_name, '_constr_h_fun_jac_uxt_zt_hess.c'];
    end
    if isfield(model_struct, 'constr_expr_h_e')
        c_files{end+1} = [model_name, '_constr_h_e_fun.c'];
        c_files{end+1} = [model_name, '_constr_h_e_fun_jac_uxt_zt.c'];
        c_files{end+1} = [model_name, '_constr_h_e_fun_jac_uxt_zt_hess.c'];
    end
end
% nonlinear least squares
if (strcmp(model_struct.cost_type, 'nonlinear_ls') || strcmp(model_struct.cost_type_e, 'nonlinear_ls'))
    % generate c for function and derivatives using casadi
    if (strcmp(opts_struct.codgen_model, 'true'))
        generate_c_code_nonlinear_least_squares(model_struct, opts_struct);
    end
    % sources list
    if isfield(model_struct, 'cost_expr_y_0')
        c_files{end+1} = [model_name, '_cost_y_0_fun.c'];
        c_files{end+1} = [model_name, '_cost_y_0_fun_jac_ut_xt.c'];
        c_files{end+1} = [model_name, '_cost_y_0_hess.c'];
    end
    if isfield(model_struct, 'cost_expr_y')
        c_files{end+1} = [model_name, '_cost_y_fun.c'];
        c_files{end+1} = [model_name, '_cost_y_fun_jac_ut_xt.c'];
        c_files{end+1} = [model_name, '_cost_y_hess.c'];
    end
    if isfield(model_struct, 'cost_expr_y_e')
        c_files{end+1} = [model_name, '_cost_y_e_fun.c'];
        c_files{end+1} = [model_name, '_cost_y_e_fun_jac_ut_xt.c'];
        c_files{end+1} = [model_name, '_cost_y_e_hess.c'];
    end
end

% external cost
if (strcmp(opts_struct.codgen_model, 'true') && ...
    ((strcmp(model_struct.cost_ext_fun_type, 'casadi') && strcmp(model_struct.cost_type, 'ext_cost')) || ...
    (strcmp(model_struct.cost_ext_fun_type_e, 'casadi') && strcmp(model_struct.cost_type_e, 'ext_cost')) || ...
    (strcmp(model_struct.cost_ext_fun_type_0, 'casadi') && strcmp(model_struct.cost_type_0, 'ext_cost'))))
    % generate c for function and derivatives using casadi    
    generate_c_code_ext_cost(model_struct, opts_struct);    
end
% external cost sources list 
if (strcmp(model_struct.cost_type, 'ext_cost') && strcmp(model_struct.cost_ext_fun_type, 'casadi') && isfield(model_struct, 'cost_expr_ext_cost'))       
    c_files{end+1} = [model_name, '_cost_ext_cost_fun.c'];
    c_files{end+1} = [model_name, '_cost_ext_cost_fun_jac.c'];
    c_files{end+1} = [model_name, '_cost_ext_cost_fun_jac_hess.c'];
end
if (strcmp(model_struct.cost_type_e, 'ext_cost') && strcmp(model_struct.cost_ext_fun_type_e, 'casadi') && isfield(model_struct, 'cost_expr_ext_cost_e'))
    c_files{end+1} = [model_name, '_cost_ext_cost_e_fun.c'];
    c_files{end+1} = [model_name, '_cost_ext_cost_e_fun_jac.c'];
    c_files{end+1} = [model_name, '_cost_ext_cost_e_fun_jac_hess.c'];    
end
if (strcmp(model_struct.cost_type_0, 'ext_cost') && strcmp(model_struct.cost_ext_fun_type_0, 'casadi') && isfield(model_struct, 'cost_expr_ext_cost_0'))
    c_files{end+1} = [model_name, '_cost_ext_cost_0_fun.c'];
    c_files{end+1} = [model_name, '_cost_ext_cost_0_fun_jac.c'];
    c_files{end+1} = [model_name, '_cost_ext_cost_0_fun_jac_hess.c'];    
end

if (strcmp(opts_struct.codgen_model, 'true'))
	for k=1:length(c_files)
		movefile(c_files{k}, opts_struct.output_dir);
	end
end

c_files_path = {};
for k=1:length(c_files)
	c_files_path{k} = fullfile(opts_struct.output_dir, c_files{k});
end

% generic external cost
if (strcmp(model_struct.cost_type, 'ext_cost') && strcmp(model_struct.cost_ext_fun_type, 'generic') &&...
    isfield(model_struct, 'cost_source_ext_cost') && isfield(model_struct, 'cost_function_ext_cost'))
    c_files_path{end+1} = model_struct.cost_source_ext_cost;        
end
if (strcmp(model_struct.cost_type_e, 'ext_cost') && strcmp(model_struct.cost_ext_fun_type_e, 'generic') &&...
    isfield(model_struct, 'cost_source_ext_cost_e') && isfield(model_struct, 'cost_function_ext_cost_e'))
    c_files_path{end+1} = model_struct.cost_source_ext_cost_e;        
end
if (strcmp(model_struct.cost_type_0, 'ext_cost') && strcmp(model_struct.cost_ext_fun_type_0, 'generic') &&...
    isfield(model_struct, 'cost_source_ext_cost_0') && isfield(model_struct, 'cost_function_ext_cost_0'))
    c_files_path{end+1} = model_struct.cost_source_ext_cost_0;        
end

% generic discrete dynamics
if (strcmp(model_struct.dyn_type, 'discrete') && strcmp(model_struct.dyn_ext_fun_type, 'generic') && ...
    isfield(model_struct, 'dyn_generic_source'))
    c_files_path{end+1} = model_struct.dyn_generic_source;
end

% check compiler
use_msvc = false;
if ~is_octave()
    mexOpts = mex.getCompilerConfigurations('C', 'Selected');
    if contains(mexOpts.ShortName, 'MSVC')
        use_msvc = true;
    end
end

ext_fun_compile_flags = opts_struct.ext_fun_compile_flags;

if use_msvc
    % get env vars for MSVC
    msvc_env = fullfile(mexOpts.Location, 'VC\Auxiliary\Build\vcvars64.bat');
    assert(isfile(msvc_env), 'Cannot find definition of MSVC env vars.');

    % assemble build command for MSVC
    out_obj_dir = [fullfile(opts_struct.output_dir), '\\'];
    out_lib = fullfile(opts_struct.output_dir, [model_name, '.dll']);
    build_cmd = sprintf('cl /O2 /EHsc /I %s /I %s /LD %s /Fo%s /Fe%s', ...
        acados_folder, fullfile(acados_folder, 'external' , 'blasfeo', 'include'), ...
        strjoin(unique(c_files_path), ' '), out_obj_dir, out_lib);

    % build
    compile_command = sprintf('"%s" & %s', msvc_env, build_cmd);
else % gcc
    % set includes
    acados_include = ['-I' acados_folder];
    blasfeo_include = ['-I' fullfile(acados_folder, 'external' , 'blasfeo', 'include')];

    if ispc
        out_lib = fullfile(opts_struct.output_dir, ['lib', model_name, '.lib']);
    else
        out_lib = fullfile(opts_struct.output_dir, ['lib', model_name, '.so']);
    end
    compile_command = ['gcc ', ext_fun_compile_flags, ' -fPIC -shared ', acados_include, ' ', blasfeo_include,...
                       ' ', strjoin(unique(c_files_path), ' '), ' -o ', out_lib];
end

compile_status = system(compile_command);
if compile_status ~= 0
    error('Compilation of model functions failed! %s %s\n%s\n\n', ...
        'Please check the compile command above and the flags therein closely.',...
        'Compile command was:', compile_command);
end

end
