function ocp_generate_c_code(obj)
    %% generate C code for CasADi functions
    % dynamics
    if (strcmp(obj.model_struct.dyn_type, 'explicit'))
        generate_c_code_explicit_ode(obj.acados_ocp_nlp_json.model);
    elseif (strcmp(obj.model_struct.dyn_type, 'implicit'))
        if (strcmp(obj.opts_struct.sim_method, 'irk'))
            opts.generate_hess = 1;
            generate_c_code_implicit_ode(...
                obj.acados_ocp_nlp_json.model, opts);
        end
    end
    % add checks for
    if ~strcmp( obj.model_struct.cost_type, 'linear_ls' ) || ...
        ~strcmp( obj.model_struct.cost_type_e, 'linear_ls' )
        error(['mex templating does only support linear_ls cost for now.',...
            ' Got cost_type: %s, cost_type_e: %s.\nNotice that it might still',...
            'be possible to solve the OCP from MATLAB.'], obj.model_struct.cost_type,...
            obj.model_struct.cost_type_e);
        % TODO: add
        % nonlinear least-squares
        % external cost
    elseif ~strcmp( obj.opts_struct.nlp_solver_exact_hessian, 'false')
        error(['mex templating does only support nlp_solver_exact_hessian = "false",',...
            'i.e. Gauss-Newton Hessian approximation for now.\n',...
            'Got nlp_solver_exact_hessian = "%s"\n. Notice that it might still be possible to solve the OCP from MATLAB.'],...
            obj.opts_struct.nlp_solver_exact_hessian);
        % TODO: add exact Hessian
    elseif strcmp( obj.model_struct.dyn_type, 'discrete')
        error('mex templating does only support discrete dynamics for now. Notice that it might still be possible to solve the OCP from MATLAB.');
        % TODO: implement
    elseif strcmp( obj.opts_struct.sim_method, 'irk_gnsf')
        error('mex templating does not support irk_gnsf integrator yet. Notice that it might still be possible to solve the OCP from MATLAB.');
        % TODO: implement
    end

    % TODO: implement nh > 0 , nh_e > 0
    if (strcmp(obj.model_struct.constr_type, 'bgh') && obj.model_struct.dim_nh > 0)
        error(['mex templating does not support general nonlinear constraints for now.',...
            'Got dim_nh: %d, must be 0.\nNotice that it might still',...
            'be possible to solve the OCP from MATLAB.\n'], obj.model_struct.dim_nh);
    end
    if (strcmp(obj.model_struct.constr_type_e, 'bgh') && obj.model_struct.dim_nh_e > 0)
        error(['mex templating does not support general nonlinear constraints for now.',...
            'Got dim_nh_e: %d, must be 0.\nNotice that it might still',...
            'be possible to solve the OCP from MATLAB.'], obj.model_struct.dim_nh_e);
    end
    if ~(strcmp(obj.opts_struct.param_scheme, 'multiple_shooting_unif_grid'))
        error(['mex templating does only support uniform discretizations for shooting nodes']);
    end

    % set include and lib path
    acados_folder = getenv('ACADOS_INSTALL_DIR');
    obj.acados_ocp_nlp_json.acados_include_path = [acados_folder, '/include'];
    obj.acados_ocp_nlp_json.acados_lib_path = [acados_folder, '/lib'];

    %% remove non-numerical data
    model.name = obj.acados_ocp_nlp_json.model.name;
    obj.acados_ocp_nlp_json.model = model;

    con_h.name = obj.acados_ocp_nlp_json.con_h.name;
    obj.acados_ocp_nlp_json.con_h = con_h;

    con_h_e.name = obj.acados_ocp_nlp_json.con_h_e.name;
    obj.acados_ocp_nlp_json.con_h_e = con_h_e;

    con_p.name = obj.acados_ocp_nlp_json.con_p.name;
    obj.acados_ocp_nlp_json.con_p = con_p;

    con_p_e.name = obj.acados_ocp_nlp_json.con_p_e.name;
    obj.acados_ocp_nlp_json.con_p_e = con_p_e;

    %% post process numerical data (mostly cast scalars to 1-dimensional cells)
    constr = obj.acados_ocp_nlp_json.constraints;
    % props = properties(constr);
    props = fieldnames(constr);
    for iprop = 1:length(props)
        thisprop = props{iprop};
        % add logic here if you want to work with select properties
        thisprop_value = constr.(thisprop);
        % add logic here if you want to do something based on the property's value
        if size(thisprop_value) == [1 1]
            constr.(thisprop) = num2cell(constr.(thisprop));
        end
    end
    obj.acados_ocp_nlp_json.constraints = constr;

    cost = obj.acados_ocp_nlp_json.cost;
    % props = properties(cost);
    props = fieldnames(cost);
    for iprop = 1:length(props)
        thisprop = props{iprop};
        %%%Add logic here if you want to work with select properties
        thisprop_value = cost.(thisprop);
        %%%Add logic here if you want to do something based on the property's value
        if norm(size(thisprop_value) - [1, 1]) == 0
            cost.(thisprop) = num2cell(cost.(thisprop));
        end
    end
    obj.acados_ocp_nlp_json.cost = cost;

    %% load JSON layout
    acados_folder = getenv('ACADOS_INSTALL_DIR');
    json_layout_filename = fullfile(acados_folder, 'interfaces',...
                                   'acados_template','acados_template','acados_layout.json');
    % if is_octave()
    addpath(fullfile(acados_folder, 'external', 'jsonlab'))
    acados_layout = loadjson(fileread(json_layout_filename));
    % else % Matlab
    %     acados_layout = jsondecode(fileread(json_layout_filename));
    % end

    %% reshape constraints
    dims = obj.acados_ocp_nlp_json.dims;
    constr = obj.acados_ocp_nlp_json.constraints;
    constr_l = acados_layout.constraints;
    fields = fieldnames(constr_l);
    for i = 1:numel(fields)
        if strcmp(constr_l.(fields{i}){1}, 'ndarray')
            if length(constr_l.(fields{i}){2}) == 1
                this_dims = [1, dims.(constr_l.(fields{i}){2}{1})];
            else
                this_dims = [dims.(constr_l.(fields{i}){2}{1}), dims.(constr_l.(fields{i}){2}{1})];
            end
            try
                constr.(fields{i}) = reshape(constr.(fields{i}), this_dims);
            catch e
                error(['error while reshaping constr.' fields{i} ...
                    ' to dimension ' num2str(this_dims), ', got ',...
                    num2str( size(constr.(fields{i}) )) , ' .\n ',...
                    e.message ]);
            end
        end
    end
    obj.acados_ocp_nlp_json.constraints = constr;

    %% reshape cost
    cost = obj.acados_ocp_nlp_json.cost;
    cost_l = acados_layout.cost;
    fields = fieldnames(cost_l);
    for i = 1:numel(fields)
        if strcmp(cost_l.(fields{i}){1}, 'ndarray')
            if length(cost_l.(fields{i}){2}) == 1
                this_dims = [1, dims.(cost_l.(fields{i}){2}{1})];
            else
                this_dims = [dims.(cost_l.(fields{i}){2}{1}), dims.(cost_l.(fields{i}){2}{2})];
            end
            try
                cost.(fields{i}) = reshape(cost.(fields{i}), this_dims);
            catch e
                    error(['error while reshaping cost.' fields{i} ...
                        ' to dimension ' num2str(this_dims), ', got ',...
                        num2str( size(cost.(fields{i}) )) , ' .\n ',...
                        e.message ]);
            end
            % convert 1-dimensional arrays to cells
            if length(cost_l.(fields{i}){2}) == 2 && (this_dims(1) == 1 || this_dims(2) == 1)
                field_as_cell = {};
                for j = 1:max(this_dims(1), this_dims(2))
                    field_as_cell{end+1} = num2cell(cost.(fields{i})(j));
                end
                cost.(fields{i}) = field_as_cell;
            end
        end
    end
    obj.acados_ocp_nlp_json.cost = cost;

    %% dump JSON file
    % if is_octave()
        % savejson does not work for classes!
        % -> consider making the acados_ocp_nlp_json properties structs directly.
        ocp_json_struct = struct(obj.acados_ocp_nlp_json);
        ocp_json_struct.dims = struct(ocp_json_struct.dims);
        ocp_json_struct.cost = struct(ocp_json_struct.cost);
        ocp_json_struct.constraints = struct(ocp_json_struct.constraints);
        ocp_json_struct.solver_options = struct(ocp_json_struct.solver_options);

        json_string = savejson('',ocp_json_struct, 'ForceRootName', 0);
    % else % Matlab
    %     json_string = jsonencode(obj.acados_ocp_nlp_json);
    % end
    fid = fopen('acados_ocp_nlp.json', 'w');
    if fid == -1, error('Cannot create JSON file'); end
    fwrite(fid, json_string, 'char');
    fclose(fid);
    %% render templated code
    acados_template_mex.render_acados_templates('acados_ocp_nlp.json')
    %% compile main
    acados_template_mex.compile_main()
end
