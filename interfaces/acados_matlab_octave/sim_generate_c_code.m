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

function sim_generate_c_code(sim, context)

    if nargin < 2
        % options for code generation
        code_gen_opts = struct();
        code_gen_opts.generate_hess = sim.solver_options.sens_hess;
        code_gen_opts.code_export_directory = 'c_generated_code'; % TODO: for OCP this is part of OCP class
        context = GenerateContext(sim.model.p_global, sim.model.name, code_gen_opts);
    else
        code_gen_opts = context.code_gen_opts;
    end

    model_dir = fullfile(pwd, code_gen_opts.code_export_directory, [sim.model.name '_model']);
    check_dir_and_create(model_dir);

    if strcmp(sim.model.dyn_ext_fun_type, 'generic')
        copyfile(fullfile(pwd, sim.model.dyn_generic_source), model_dir);
        context.add_external_function_file(ocp.model.dyn_generic_source, model_dir);

    elseif strcmp(sim.model.dyn_ext_fun_type, 'casadi')
        import casadi.*
        check_casadi_version();
        switch sim.solver_options.integrator_type
            case 'ERK'
                generate_c_code_explicit_ode(context, sim.model, model_dir);
            case 'IRK'
                generate_c_code_implicit_ode(context, sim.model, model_dir);
            case 'GNSF'
                generate_c_code_gnsf(context, sim.model, model_dir);
            case 'DISCRETE'
                generate_c_code_discrete_dynamics(context, sim.model, model_dir);
            otherwise
                error('Unknown integrator type.')
        end
    else
        error('Unknown dyn_ext_fun_type.')
    end

    context.finalize();
    sim.external_function_files_model = context.get_external_function_file_list(false);

    %% remove CasADi objects from model
    model.name = sim.model.name;
    model.dyn_ext_fun_type = sim.model.dyn_ext_fun_type;
    model.dyn_generic_source = sim.model.dyn_generic_source;
    model.dyn_disc_fun_jac_hess = sim.model.dyn_disc_fun_jac_hess;
    model.dyn_disc_fun_jac = sim.model.dyn_disc_fun_jac;
    model.dyn_disc_fun = sim.model.dyn_disc_fun;
    model.gnsf_nontrivial_f_LO = sim.model.gnsf_nontrivial_f_LO;
    model.gnsf_purely_linear = sim.model.gnsf_purely_linear;
    sim.model = model;
    %% post process numerical data (mostly cast scalars to 1-dimensional cells)
    dims = sim.dims;

    %% load JSON layout
    acados_folder = getenv('ACADOS_INSTALL_DIR');
    addpath(fullfile(acados_folder, 'external', 'jsonlab'))

    % parameter values
    sim.parameter_values = reshape(num2cell(sim.parameter_values), [1, dims.np]);

    %% dump JSON file
    sim_json_struct = sim.struct();
    sim_json_struct.dims = sim.dims.struct();
    sim_json_struct.solver_options = sim.solver_options.struct();

    % add compilation information to json
    libs = loadjson(fileread(fullfile(acados_folder, 'lib', 'link_libs.json')));
    sim_json_struct.acados_link_libs = libs;
    if ismac
        sim_json_struct.os = 'mac';
    elseif isunix
        sim_json_struct.os = 'unix';
    else
        sim_json_struct.os = 'pc';
    end

    json_string = savejson('', sim_json_struct, 'ForceRootName', 0);

    fid = fopen(sim.json_file, 'w');
    if fid == -1, error('Cannot create JSON file'); end
    fwrite(fid, json_string, 'char');
    fclose(fid);
    %% render templated code
    acados_template_mex.render_acados_sim_templates(sim.json_file)
    acados_template_mex.compile_sim_shared_lib(sim.code_export_directory)
end

