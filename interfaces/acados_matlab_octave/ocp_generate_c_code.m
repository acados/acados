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

function ocp_generate_c_code(ocp)
    %% create folder

    check_dir_and_create(fullfile(pwd,'c_generated_code'));

    %% generate C code for CasADi functions / copy external functions
    cost = ocp.cost;
    solver_opts = ocp.solver_options;
    constraints = ocp.constraints;
    dims = ocp.dims;

    % options for code generation
    code_gen_opts = struct();
    code_gen_opts.generate_hess = strcmp(solver_opts.hessian_approx, 'EXACT');
    code_gen_opts.with_solution_sens_wrt_params = solver_opts.with_solution_sens_wrt_params;
    code_gen_opts.with_value_sens_wrt_params = solver_opts.with_value_sens_wrt_params;

    % dynamics
    % model dir is always need, other dirs are  only created if necessary
    model_dir = fullfile(pwd, 'c_generated_code', [ocp.name '_model']);
    check_dir_and_create(model_dir);

    if (strcmp(solver_opts.integrator_type, 'ERK'))
        generate_c_code_explicit_ode(ocp.model, code_gen_opts, model_dir);
    elseif (strcmp(solver_opts.integrator_type, 'IRK')) && strcmp(ocp.model.dyn_ext_fun_type, 'casadi')
        generate_c_code_implicit_ode(ocp.model, code_gen_opts, model_dir);
    elseif (strcmp(solver_opts.integrator_type, 'GNSF'))
        generate_c_code_gnsf(ocp.model, code_gen_opts, model_dir);
    elseif (strcmp(solver_opts.integrator_type, 'DISCRETE')) && strcmp(ocp.model.dyn_ext_fun_type, 'casadi')
        generate_c_code_disc_dyn(ocp.model, code_gen_opts, model_dir);
    end
    if strcmp(ocp.model.dyn_ext_fun_type, 'generic')
        copyfile( fullfile(pwd, ocp.model.dyn_generic_source), model_dir);
    end

    stage_types = {'initial', 'path', 'terminal'};

    % cost
    cost_types = {cost.cost_type_0, cost.cost_type, cost.cost_type_e};
    cost_ext_fun_types = {cost.cost_ext_fun_type_0, cost.cost_ext_fun_type, cost.cost_ext_fun_type_e};
    cost_dir = fullfile(pwd, 'c_generated_code', [ocp.name '_cost']);

    for i = 1:3
        if strcmp(cost_types{i}, 'NONLINEAR_LS')
            generate_c_code_nonlinear_least_squares( ocp.model, cost_dir, stage_types{i} );

        elseif strcmp(cost_types{i}, 'CONVEX_OVER_NONLINEAR')
            % TODO
            error("Convex-over-nonlinear cost is not implemented yet.")

        elseif strcmp(cost_types{i}, 'EXTERNAL')
            if strcmp(cost_ext_fun_types{i}, 'casadi')
                generate_c_code_ext_cost(ocp.model, cost_dir, stage_types{i});
            elseif strcmp(cost_ext_fun_types{i}, 'generic')
                setup_generic_cost(cost, cost_dir, stage_types{i})
            else
                error('Unknown value for cost_ext_fun_types %s', cost_ext_fun_types{i});
            end
        end
    end


    % constraints
    constraints_types = {constraints.constr_type_0, constraints.constr_type, constraints.constr_type_e};
    constraints_dims = {dims.nh_0, dims.nh, dims.nh_e};
    constraints_dir = fullfile(pwd, 'c_generated_code', [ocp.name '_constraints']);

    for i = 1:3
        if strcmp(constraints_types{i}, 'BGH') && constraints_dims{i} > 0
            generate_c_code_nonlinear_constr( ocp.model, constraints_dir, stage_types{i} );
        end
    end

    %% remove CasADi objects from model
    model_without_expr = struct();
    model_without_expr.name = ocp.model.name;
    model_without_expr.dyn_ext_fun_type = ocp.model.dyn_ext_fun_type;
    model_without_expr.dyn_generic_source = ocp.model.dyn_generic_source;
    model_without_expr.dyn_disc_fun_jac_hess = ocp.model.dyn_disc_fun_jac_hess;
    model_without_expr.dyn_disc_fun_jac = ocp.model.dyn_disc_fun_jac;
    model_without_expr.dyn_disc_fun = ocp.model.dyn_disc_fun;
    model_without_expr.gnsf_nontrivial_f_LO = ocp.model.gnsf_nontrivial_f_LO;
    model_without_expr.gnsf_purely_linear = ocp.model.gnsf_purely_linear;
    ocp.model = model_without_expr;

    %% post process numerical data (mostly cast scalars to 1-dimensional cells)
    props = fieldnames(ocp.constraints);
    disable_last_warning();  % show warning for struct conversion only once
    for iprop = 1:length(props)
        this_prop = props{iprop};
        % add logic here if you want to do something based on the property's value
        if strcmp(this_prop, 'x0')
            continue;
        end
        % add logic here if you want to work with select properties
        this_prop_value = ocp.constraints.(this_prop);
        if all(size(this_prop_value) == [1 1])
            ocp.constraints.(this_prop) = num2cell(ocp.constraints.(this_prop));
        end
    end

    props = fieldnames(ocp.cost);
    for iprop = 1:length(props)
        this_prop = props{iprop};
        % add logic here if you want to work with select properties
        this_prop_value = ocp.cost.(this_prop);
        % add logic here if you want to do something based on the property's value
        if all(size(this_prop_value) == [1, 1])
            ocp.cost.(this_prop) = num2cell(ocp.cost.(this_prop));
        end
    end

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
    constr_layout = acados_layout.constraints;
    fields = fieldnames(constr_layout);
    for i = 1:numel(fields)
        if strcmp(constr_layout.(fields{i}){1}, 'ndarray')
            property_dim_names = constr_layout.(fields{i}){2};
            if length(property_dim_names) == 1 % vector
                this_dims = [1, ocp.dims.(property_dim_names{1})];
            else % matrix
                this_dims = [ocp.dims.(property_dim_names{1}), ocp.dims.(property_dim_names{2})];
            end
            try
                ocp.constraints.(fields{i}) = reshape(ocp.constraints.(fields{i}), this_dims);
            catch e
                keyboard
                error(['error while reshaping constraints.' fields{i} ...
                    ' to dimension ' num2str(this_dims), ', got ',...
                    num2str( size(ocp.constraints.(fields{i}) )) , 10,...
                    e.message ]);
            end
            if this_dims(1) == 1 && length(property_dim_names) ~= 1 % matrix with 1 row
                ocp.constraints.(fields{i}) = {ocp.constraints.(fields{i})};
            end
        end
    end

    %% reshape cost
    cost_layout = acados_layout.cost;
    fields = fieldnames(cost_layout);
    for i = 1:numel(fields)
        if strcmp(cost_layout.(fields{i}){1}, 'ndarray')
            property_dim_names = cost_layout.(fields{i}){2};
            if length(property_dim_names) == 1 % vector
                this_dims = [1, ocp.dims.(property_dim_names{1})];
            else % matrix
                this_dims = [ocp.dims.(property_dim_names{1}), ocp.dims.(property_dim_names{2})];
            end
            if ~isempty(ocp.cost.(fields{i}))
                try
                    ocp.cost.(fields{i}) = reshape(ocp.cost.(fields{i}), this_dims);
                catch e
                    error(['error while reshaping cost.' fields{i} ...
                        ' to dimension ' num2str(this_dims), ', got ',...
                        num2str( size(ocp.cost.(fields{i}) )) , 10,...
                        e.message ]);
                end
                if this_dims(1) == 1 && length(property_dim_names) ~= 1 % matrix with 1 row
                    ocp.cost.(fields{i}) = {ocp.cost.(fields{i})};
                end
            end
        elseif strcmp(cost_layout.(fields{i}){1}, 'int')
            ocp.cost.(fields{i}) = ocp.cost.(fields{i}){1};
        end
    end

    %% reshape opts
    opts = ocp.solver_options;
    opts_layout = acados_layout.solver_options;
    fields = fieldnames(opts_layout);
    for i = 1:numel(fields)
        if strcmp(opts_layout.(fields{i}){1}, 'ndarray')
            property_dim_names = opts_layout.(fields{i}){2};
            if length(property_dim_names) == 1 % vector
                this_dims = [1, ocp.dims.(property_dim_names{1})];
            else % matrix
                this_dims = [ocp.dims.(property_dim_names{1}), ocp.dims.(property_dim_names{2})];
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
    opts.time_steps = reshape(num2cell(opts.time_steps), [1, ocp.dims.N]);
    opts.sim_method_num_stages = reshape(num2cell(opts.sim_method_num_stages), [1, ocp.dims.N]);
    opts.sim_method_num_steps = reshape(num2cell(opts.sim_method_num_steps), [1, ocp.dims.N]);
    opts.sim_method_jac_reuse = reshape(num2cell(opts.sim_method_jac_reuse), [1, ocp.dims.N]);

    ocp.solver_options = opts;

    % parameter values
    try
        ocp.parameter_values = reshape(num2cell(ocp.parameter_values), [1, ocp.dims.np]);
    catch e
        error(['error while reshaping parameter_values'  ...
                ' to dimension ' num2str([1, ocp.dims.np]) ', got ',...
                num2str(size(ocp.parameter_values)) , 10,...
                e.message ]);
    end
    %% dump JSON file
    % if is_octave()
        % savejson does not work for classes!
        % -> consider making the ocp properties structs directly.
        ocp_json_struct = orderfields(ocp.struct());
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
    %     json_string = jsonencode(ocp);
    % end
    fid = fopen(ocp.json_file, 'w');
    if fid == -1, error('Cannot create JSON file'); end
    fwrite(fid, json_string, 'char');
    fclose(fid);
    %% render templated code
    ocp.render_templates()
    acados_template_mex.compile_ocp_shared_lib(ocp.code_export_directory)
end


function setup_generic_cost(cost, target_dir, stage_type)

    if strcmp(stage_type, 'initial')
        cost_source_ext_cost = cost.cost_source_ext_cost_0;
    elseif strcmp(stage_type, 'path')
        cost_source_ext_cost = cost.cost_source_ext_cost;
    elseif strcmp(stage_type, 'terminal')
        cost_source_ext_cost = cost.cost_source_ext_cost_e;
    else
        error('Unknown stage_type.')
    end

    check_dir_and_create(target_dir);
    copyfile(fullfile(pwd, cost_source_ext_cost), target_dir);
end

