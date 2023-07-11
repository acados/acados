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

function sim_generate_casadi_ext_fun(model_struct, opts_struct)

model_name = model_struct.name;

c_files = {};
if (strcmp(opts_struct.method, 'erk'))
    % generate c for function and derivatives using casadi
    if (strcmp(opts_struct.codgen_model, 'true'))
        generate_c_code_explicit_ode(model_struct, opts_struct);
    end
    % compile the code in a shared library
    c_files{end+1} = [model_name, '_dyn_expl_ode_fun.c'];
    c_files{end+1} = [model_name, '_dyn_expl_vde_forw.c'];
    c_files{end+1} = [model_name, '_dyn_expl_vde_adj.c'];
    c_files{end+1} = [model_name, '_dyn_expl_ode_hess.c'];
elseif (strcmp(opts_struct.method, 'irk'))
    % generate c for function and derivatives using casadi
    if (strcmp(opts_struct.codgen_model, 'true'))
        generate_c_code_implicit_ode(model_struct, opts_struct);
    end
    % compile the code in a shared library
    c_files{end+1} = [model_name, '_dyn_impl_dae_fun.c'];
    c_files{end+1} = [model_name, '_dyn_impl_dae_fun_jac_x_xdot_z.c'];
    c_files{end+1} = [model_name, '_dyn_impl_dae_fun_jac_x_xdot_u.c'];
    c_files{end+1} = [model_name, '_dyn_impl_dae_jac_x_xdot_u_z.c'];
    if strcmp(opts_struct.sens_hess, 'true')
        c_files{end+1} = [model_name, '_dyn_impl_dae_hess.c'];
    end
elseif (strcmp(opts_struct.method, 'irk_gnsf'))
    % generate c for function and derivatives using casadi
    if (strcmp(opts_struct.codgen_model, 'true'))
        generate_c_code_gnsf(model_struct); %, opts_struct);
    end
    % compile the code in a shared library
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
    fprintf('\nsim_generate_casadi_ext_fun: method not supported: %s\n', opts_struct.method);
    return;
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

% Store the current PATH environment variable value
origEnvPath = getenv('PATH');

% check compiler
use_msvc = false;
if ~is_octave()
    mexOpts = mex.getCompilerConfigurations('C', 'Selected');
    if contains(mexOpts.ShortName, 'MSVC')
        use_msvc = true;
    else
        % Get mex C compiler configuration and extract the location
        pathToCompilerLocation = mexOpts.Location;
        % Set environment PATH variable for this Matlab session such that
        % configured mex C compiler is prioritized
        setenv('PATH', [fullfile(pathToCompilerLocation,'bin') ';' origEnvPath]);
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
    build_cmd = sprintf('cl /O2 /EHsc /LD %s /Fo%s /Fe%s', ...
        strjoin(unique(c_files_path), ' '), out_obj_dir, out_lib);

    compile_command = sprintf('"%s" & %s', msvc_env, build_cmd);
else % gcc
    if ispc
        out_lib = fullfile(opts_struct.output_dir, ['lib', model_name, '.lib']);
    else
        out_lib = fullfile(opts_struct.output_dir, ['lib', model_name, '.so']);
    end
    compile_command = ['gcc ', ext_fun_compile_flags, ' -fPIC -shared ', strjoin(unique(c_files_path), ' '), ' -o ', out_lib];
end

% Store the PATH environment variable used during compilation for error reporting
envPath = getenv('PATH');

compile_status = system(compile_command);

% Reset the environment PATH variable to its original value before potentially
% raising an error to ensure that the path environment variable is clean
setenv('PATH',origEnvPath);

if compile_status ~= 0
    error('Compilation of model functions failed! %s %s\n%s\n\n', ...
        'Please check the compile command above and the flags therein closely.',...
        'Compile command was:', compile_command, '\n', ...
        'Environment path was: ', envPath);
end

end

