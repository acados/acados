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

function sim_generate_c_code(obj)
    %% create folder
    if ~exist(fullfile(pwd,'c_generated_code'), 'dir')
        mkdir(fullfile(pwd, 'c_generated_code'))
    end

    %% generate C code for CasADi functions / copy external functions
    % dynamics
    if (strcmp(obj.model_struct.dyn_type, 'explicit'))
        generate_c_code_explicit_ode(obj.sim.model);
    elseif (strcmp(obj.model_struct.dyn_type, 'implicit'))
        if (strcmp(obj.opts_struct.method, 'irk'))
            opts.sens_hess = 'true';
            generate_c_code_implicit_ode(...
                obj.sim.model, opts);
        elseif (strcmp(obj.opts_struct.method, 'irk_gnsf'))
            generate_c_code_gnsf(...
                obj.sim.model);
        end
    elseif (strcmp(obj.model_struct.dyn_type, 'discrete'))
        generate_c_code_disc_dyn(obj.sim.model);
    end
    if strcmp(obj.sim.model.dyn_ext_fun_type, 'generic')
        copyfile( fullfile(pwd, obj.sim.model.dyn_generic_source),...
            fullfile(pwd, 'c_generated_code', [obj.model_struct.name '_model']));
    end


    %% remove CasADi objects from model
    model.name = obj.sim.model.name;
    model.dyn_ext_fun_type = obj.sim.model.dyn_ext_fun_type;
    model.dyn_generic_source = obj.sim.model.dyn_generic_source;
    model.dyn_disc_fun_jac_hess = obj.sim.model.dyn_disc_fun_jac_hess;
    model.dyn_disc_fun_jac = obj.sim.model.dyn_disc_fun_jac;
    model.dyn_disc_fun = obj.sim.model.dyn_disc_fun;
    model.gnsf.nontrivial_f_LO = obj.sim.model.gnsf.nontrivial_f_LO;
    model.gnsf.purely_linear = obj.sim.model.gnsf.purely_linear;
    obj.sim.model = model;
    %% post process numerical data (mostly cast scalars to 1-dimensional cells)
    dims = obj.sim.dims;

    %% load JSON layout
    acados_folder = getenv('ACADOS_INSTALL_DIR');
    json_layout_filename = fullfile(acados_folder, 'interfaces',...
                                   'acados_matlab_octave', ...
                                   'acados_template_mex', '+acados_template_mex','acados_sim_layout.json');
    % if is_octave()
    addpath(fullfile(acados_folder, 'external', 'jsonlab'))
    acados_sim_layout = loadjson(fileread(json_layout_filename));
    % else % Matlab
    %     acados_sim_layout = jsondecode(fileread(json_layout_filename));
    % end

    %% reshape opts
    opts = obj.sim.sim_options;
    opts_layout = acados_sim_layout.solver_options;
    fields = fieldnames(opts_layout);
    for i = 1:numel(fields)
        if strcmp(opts_layout.(fields{i}){1}, 'ndarray')
            property_dim_names = opts_layout.(fields{i}){2};
            if length(property_dim_names) == 1 % vector
                this_dims = [1, dims.(property_dim_names{1})];
            else % matrix
                this_dims = [dims.(property_dim_names{1}), dims.(property_dim_names{2})];
            end
            try
                opts.(fields{i}) = reshape(opts.(fields{i}), this_dims);
            catch e
                error(['error while reshaping opts.' fields{i} ...
                    ' to dimension ' num2str(this_dims), ', got ',...
                    num2str( size(opts.(fields{i}) )) , 10,...
                    e.message ]);
            end
            if this_dims(1) == 1 && length(property_dim_names) ~= 1 % matrix with 1 row
                opts.(fields{i}) = {opts.(fields{i})};
            end
        end
    end
    obj.sim.sim_options = opts;

    % parameter values
    obj.sim.parameter_values = reshape(num2cell(obj.sim.parameter_values), [ 1, dims.np]);

    %% dump JSON file
    % if is_octave()
        % savejson does not work for classes!
        % -> consider making the sim properties structs directly.
        sim_json_struct = obj.sim.struct();
        sim_json_struct.dims = obj.sim.dims.struct();
        sim_json_struct.solver_options = obj.sim.sim_options;

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
    % else % Matlab
    %     json_string = jsonencode(obj.sim);
    % end
    fid = fopen(obj.sim.json_file, 'w');
    if fid == -1, error('Cannot create JSON file'); end
    fwrite(fid, json_string, 'char');
    fclose(fid);
    %% render templated code
    acados_template_mex.render_acados_sim_templates(obj.sim.json_file)
    acados_template_mex.compile_sim_shared_lib(obj.sim.code_export_directory)
end
