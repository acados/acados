function ocp_generate_c_code(obj)
    %% create folder
    if ~exist(fullfile(pwd,'c_generated_code'), 'dir')
        mkdir(fullfile(pwd, 'c_generated_code'))
    end

    %% generate C code for CasADi functions / copy external functions
    % dynamics
    if (strcmp(obj.model_struct.dyn_type, 'explicit'))
        generate_c_code_explicit_ode(obj.acados_sim_json.model);
    elseif (strcmp(obj.model_struct.dyn_type, 'implicit'))
        if (strcmp(obj.opts_struct.method, 'irk'))
            opts.sens_hess = 'true';
            generate_c_code_implicit_ode(...
                obj.acados_sim_json.model, opts);
        elseif (strcmp(obj.opts_struct.method, 'irk_gnsf'))
            generate_c_code_gnsf(...
                obj.acados_sim_json.model);
        end
    elseif (strcmp(obj.model_struct.dyn_type, 'discrete'))
        generate_c_code_disc_dyn(obj.acados_sim_json.model);
    end
    if strcmp(obj.acados_sim_json.model.dyn_ext_fun_type, 'generic')
        copyfile( fullfile(pwd, obj.acados_sim_json.model.dyn_generic_source),...
            fullfile(pwd, 'c_generated_code', [obj.model_struct.name '_model']));
    end


    %% remove CasADi objects from model
    model.name = obj.acados_sim_json.model.name;
    model.dyn_ext_fun_type = obj.acados_sim_json.model.dyn_ext_fun_type;
    model.dyn_generic_source = obj.acados_sim_json.model.dyn_generic_source;
    model.dyn_disc_fun_jac_hess = obj.acados_sim_json.model.dyn_disc_fun_jac_hess;
    model.dyn_disc_fun_jac = obj.acados_sim_json.model.dyn_disc_fun_jac;
    model.dyn_disc_fun = obj.acados_sim_json.model.dyn_disc_fun;
    model.gnsf.nontrivial_f_LO = obj.acados_sim_json.model.gnsf.nontrivial_f_LO;
    model.gnsf.purely_linear = obj.acados_sim_json.model.gnsf.purely_linear;
    obj.acados_sim_json.model = model;
    %% post process numerical data (mostly cast scalars to 1-dimensional cells)
    dims = obj.acados_sim_json.dims;

    %% load JSON layout
    acados_folder = getenv('ACADOS_INSTALL_DIR');
    json_layout_filename = fullfile(acados_folder, 'interfaces',...
                                   'acados_template','acados_template','acados_sim_layout.json');
    % if is_octave()
    addpath(fullfile(acados_folder, 'external', 'jsonlab'))
    acados_sim_layout = loadjson(fileread(json_layout_filename));
    % else % Matlab
    %     acados_sim_layout = jsondecode(fileread(json_layout_filename));
    % end

    %% reshape opts
    opts = obj.acados_sim_json.sim_options;
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
    obj.acados_sim_json.sim_options = opts;

    % parameter values
    obj.acados_sim_json.parameter_values = reshape(num2cell(obj.acados_sim_json.parameter_values), [ 1, dims.np]);

    %% dump JSON file
    % if is_octave()
        % savejson does not work for classes!
        % -> consider making the acados_sim_json properties structs directly.
        sim_json_struct = obj.acados_sim_json.struct();
        sim_json_struct.dims = sim_json_struct.dims;
        sim_json_struct.sim_options = sim_json_struct.sim_options;

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
    %     json_string = jsonencode(obj.acados_sim_json);
    % end
    fid = fopen(obj.acados_sim_json.json_file, 'w');
    if fid == -1, error('Cannot create JSON file'); end
    fwrite(fid, json_string, 'char');
    fclose(fid);
    %% render templated code
    acados_template_mex.render_acados_templates(obj.acados_sim_json.json_file)
    acados_template_mex.compile_main(obj.acados_sim_json.code_export_directory)
end
