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

function ocp_generate_c_code(obj)
    %% create folder
    if ~exist(fullfile(pwd,'c_generated_code'), 'dir')
        mkdir(fullfile(pwd, 'c_generated_code'))
    end
    %% generate C code for CasADi functions / copy external functions
    % dynamics
    if (strcmp(obj.model_struct.dyn_type, 'explicit'))
        generate_c_code_explicit_ode(obj.ocp.model);
    elseif (strcmp(obj.model_struct.dyn_type, 'implicit'))
        if (strcmp(obj.opts_struct.sim_method, 'irk'))
            opts.sens_hess = 'true';
            generate_c_code_implicit_ode(...
                obj.ocp.model, opts);
        elseif (strcmp(obj.opts_struct.sim_method, 'irk_gnsf'))
            generate_c_code_gnsf(...
                obj.ocp.model);
        end
    elseif (strcmp(obj.model_struct.dyn_type, 'discrete'))
        generate_c_code_disc_dyn(obj.ocp.model);
    end
    if strcmp(obj.ocp.model.dyn_ext_fun_type, 'generic')
        copyfile( fullfile(pwd, obj.ocp.model.dyn_generic_source),...
            fullfile(pwd, 'c_generated_code', [obj.model_struct.name '_model']));
    end

    % cost
    if (strcmp(obj.model_struct.cost_type, 'nonlinear_ls') || ...
        strcmp(obj.model_struct.cost_type_0, 'nonlinear_ls') || strcmp(obj.model_struct.cost_type_e, 'nonlinear_ls'))
        generate_c_code_nonlinear_least_squares( obj.model_struct, obj.opts_struct,...
              fullfile(pwd, 'c_generated_code', [obj.model_struct.name '_cost']) );
    elseif (strcmp(obj.model_struct.cost_type, 'ext_cost') || ...
            strcmp(obj.model_struct.cost_type_0, 'ext_cost') || strcmp(obj.model_struct.cost_type_e, 'ext_cost'))
            generate_c_code_ext_cost( obj.model_struct, obj.opts_struct,...
              fullfile(pwd, 'c_generated_code', [obj.model_struct.name '_cost']) );
    end
    if (strcmp(obj.ocp.cost.cost_ext_fun_type_0, 'generic'))
        copyfile(fullfile(pwd, obj.ocp.cost.cost_source_ext_cost_0), ...
            fullfile(pwd, 'c_generated_code', [obj.model_struct.name '_cost']));
    end
    if (strcmp(obj.ocp.cost.cost_ext_fun_type, 'generic'))
        copyfile(fullfile(pwd, obj.ocp.cost.cost_source_ext_cost), ...
            fullfile(pwd, 'c_generated_code', [obj.model_struct.name '_cost']));
    end
    if (strcmp(obj.ocp.cost.cost_ext_fun_type_e, 'generic'))
        copyfile(fullfile(pwd, obj.ocp.cost.cost_source_ext_cost_e), ...
            fullfile(pwd, 'c_generated_code', [obj.model_struct.name '_cost']));
    end
    % constraints
    if ((strcmp(obj.model_struct.constr_type, 'bgh') && obj.model_struct.dim_nh > 0) || ...
        (strcmp(obj.model_struct.constr_type_0, 'bgh') && obj.model_struct.dim_nh_0 > 0) || ...
        (strcmp(obj.model_struct.constr_type_e, 'bgh') && obj.model_struct.dim_nh_e > 0))
        generate_c_code_nonlinear_constr( obj.model_struct, obj.opts_struct,...
              fullfile(pwd, 'c_generated_code', [obj.model_struct.name '_constraints']) );
    end

    % set include and lib path
    acados_folder = getenv('ACADOS_INSTALL_DIR');
    obj.ocp.acados_include_path = [acados_folder, '/include'];
    obj.ocp.acados_lib_path = [acados_folder, '/lib'];

    %% remove CasADi objects from model
    model.name = obj.ocp.model.name;
    model.dyn_ext_fun_type = obj.ocp.model.dyn_ext_fun_type;
    model.dyn_generic_source = obj.ocp.model.dyn_generic_source;
    model.dyn_disc_fun_jac_hess = obj.ocp.model.dyn_disc_fun_jac_hess;
    model.dyn_disc_fun_jac = obj.ocp.model.dyn_disc_fun_jac;
    model.dyn_disc_fun = obj.ocp.model.dyn_disc_fun;
    model.gnsf.nontrivial_f_LO = obj.ocp.model.gnsf.nontrivial_f_LO;
    model.gnsf.purely_linear = obj.ocp.model.gnsf.purely_linear;
    obj.ocp.model = model;
    %% post process numerical data (mostly cast scalars to 1-dimensional cells)
    dims = obj.ocp.dims;

    constr = obj.ocp.constraints;
    props = fieldnames(constr);
    disable_last_warning();  % show warning for struct conversion only once
    for iprop = 1:length(props)
        this_prop = props{iprop};
        % add logic here if you want to work with select properties
        this_prop_value = constr.(this_prop);
        % add logic here if you want to do something based on the property's value
        if all(size(this_prop_value) == [1 1])
            constr.(this_prop) = num2cell(constr.(this_prop));
        end
    end
    obj.ocp.constraints = constr;

    cost = obj.ocp.cost;
    props = fieldnames(cost);
    for iprop = 1:length(props)
        this_prop = props{iprop};
        % add logic here if you want to work with select properties
        this_prop_value = cost.(this_prop);
        % add logic here if you want to do something based on the property's value
        if all(size(this_prop_value) == [1, 1])
            cost.(this_prop) = num2cell(cost.(this_prop));
        end
    end

    % for cost type not LINEAR_LS, fill matrices with zeros
    if ~strcmp(cost.cost_type, 'LINEAR_LS')
        cost.Vx_0 = zeros(dims.ny_0, dims.nx);
        cost.Vu_0 = zeros(dims.ny_0, dims.nu);
        cost.Vz_0 = zeros(dims.ny_0, dims.nz);
    end
    if ~strcmp(cost.cost_type, 'LINEAR_LS')
        cost.Vx = zeros(dims.ny, dims.nx);
        cost.Vu = zeros(dims.ny, dims.nu);
        cost.Vz = zeros(dims.ny, dims.nz);
    end
    if ~strcmp(cost.cost_type_e, 'LINEAR_LS')
        cost.Vx_e = zeros(dims.ny_e, dims.nx);
    end
    obj.ocp.cost = cost;

    %% load JSON layout
    acados_folder = getenv('ACADOS_INSTALL_DIR');
    json_layout_filename = fullfile(acados_folder, 'interfaces',...
                                   'acados_matlab_octave', ...
                                   'acados_template_mex', '+acados_template_mex','acados_ocp_layout.json');
    % if is_octave()
    addpath(fullfile(acados_folder, 'external', 'jsonlab'))
    acados_layout = loadjson(fileread(json_layout_filename));
    % else % Matlab
    %     acados_layout = jsondecode(fileread(json_layout_filename));
    % end

    %% reshape constraints
    constr = obj.ocp.constraints;
    constr_layout = acados_layout.constraints;
    fields = fieldnames(constr_layout);
    for i = 1:numel(fields)
        if strcmp(constr_layout.(fields{i}){1}, 'ndarray')
            property_dim_names = constr_layout.(fields{i}){2};
            if length(property_dim_names) == 1 % vector
                this_dims = [1, dims.(property_dim_names{1})];
            else % matrix
                this_dims = [dims.(property_dim_names{1}), dims.(property_dim_names{2})];
            end
            try
                constr.(fields{i}) = reshape(constr.(fields{i}), this_dims);
            catch e
                error(['error while reshaping constr.' fields{i} ...
                    ' to dimension ' num2str(this_dims), ', got ',...
                    num2str( size(constr.(fields{i}) )) , 10,...
                    e.message ]);
            end
            if this_dims(1) == 1 && length(property_dim_names) ~= 1 % matrix with 1 row
                constr.(fields{i}) = {constr.(fields{i})};
            end
        end
    end
    obj.ocp.constraints = constr;

    %% reshape cost
    cost = obj.ocp.cost;
    cost_layout = acados_layout.cost;
    fields = fieldnames(cost_layout);
    for i = 1:numel(fields)
        if strcmp(cost_layout.(fields{i}){1}, 'ndarray')
            property_dim_names = cost_layout.(fields{i}){2};
            if length(property_dim_names) == 1 % vector
                this_dims = [1, dims.(property_dim_names{1})];
            else % matrix
                this_dims = [dims.(property_dim_names{1}), dims.(property_dim_names{2})];
            end
            try
                cost.(fields{i}) = reshape(cost.(fields{i}), this_dims);
            catch e
                error(['error while reshaping cost.' fields{i} ...
                    ' to dimension ' num2str(this_dims), ', got ',...
                    num2str( size(cost.(fields{i}) )) , 10,...
                    e.message ]);
            end
            if this_dims(1) == 1 && length(property_dim_names) ~= 1 % matrix with 1 row
                cost.(fields{i}) = {cost.(fields{i})};
            end
        elseif strcmp(cost_layout.(fields{i}){1}, 'int')
            cost.(fields{i}) = cost.(fields{i}){1};
        end
    end
    obj.ocp.cost = cost;

    %% reshape opts
    opts = obj.ocp.solver_options;
    opts_layout = acados_layout.solver_options;
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
    opts.time_steps = reshape(num2cell(opts.time_steps), [1, dims.N]);
    opts.sim_method_num_stages = reshape(num2cell(opts.sim_method_num_stages), [1, dims.N]);
    opts.sim_method_num_steps = reshape(num2cell(opts.sim_method_num_steps), [1, dims.N]);
    opts.sim_method_jac_reuse = reshape(num2cell(opts.sim_method_jac_reuse), [1, dims.N]);
    obj.ocp.solver_options = opts;

    % parameter values
    obj.ocp.parameter_values = reshape(num2cell(obj.ocp.parameter_values), [ 1, dims.np]);

    %% dump JSON file
    % if is_octave()
        % savejson does not work for classes!
        % -> consider making the ocp properties structs directly.
        ocp_json_struct = orderfields(obj.ocp.struct());
        ocp_json_struct.dims = orderfields(ocp_json_struct.dims.struct());
        ocp_json_struct.cost = orderfields(ocp_json_struct.cost.struct());
        ocp_json_struct.constraints = orderfields(ocp_json_struct.constraints.struct());
        ocp_json_struct.solver_options = orderfields(ocp_json_struct.solver_options.struct());

        % add compilation information to json
        libs = loadjson(fileread(fullfile(acados_folder, 'lib', 'link_libs.json')));
        ocp_json_struct.acados_link_libs = orderfields(libs);
        if ismac
            ocp_json_struct.os = 'mac';
        elseif isunix
            ocp_json_struct.os = 'unix';
        else
            ocp_json_struct.os = 'pc';
        end

        json_string = savejson('',ocp_json_struct, 'ForceRootName', 0);
    % else % Matlab
    %     json_string = jsonencode(obj.ocp);
    % end
    fid = fopen(obj.ocp.json_file, 'w');
    if fid == -1, error('Cannot create JSON file'); end
    fwrite(fid, json_string, 'char');
    fclose(fid);
    %% render templated code
    acados_template_mex.render_acados_templates(obj.ocp.json_file)
    acados_template_mex.compile_ocp_shared_lib(obj.ocp.code_export_directory)
end
