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

function sim_generate_c_code(sim)
    %% create folders
    check_dir_and_create(fullfile(pwd,'c_generated_code'));

    model_dir = fullfile(pwd, 'c_generated_code', [sim.model.name '_model']);
    check_dir_and_create(model_dir);

    code_gen_opts = struct('generate_hess', sim.solver_options.sens_hess);
    %% generate C code for CasADi functions / copy external functions
    % dynamics
    if (strcmp(sim.solver_options.integrator_type, 'ERK'))
        generate_c_code_explicit_ode(sim.model, code_gen_opts, model_dir);
    elseif (strcmp(sim.solver_options.integrator_type, 'IRK'))
        generate_c_code_implicit_ode(sim.model, code_gen_opts, model_dir);
    elseif (strcmp(sim.solver_options.integrator_type, 'GNSF'))
        generate_c_code_gnsf(sim.model, code_gen_opts, model_dir);
    elseif (strcmp(sim.solver_options.integrator_type, 'DISCRETE'))
        generate_c_code_discrete_dynamics(sim.model, code_gen_opts, model_dir);
    end
    if strcmp(sim.model.dyn_ext_fun_type, 'generic')
        copyfile(fullfile(pwd, sim.model.dyn_generic_source), model_dir);
    end


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

