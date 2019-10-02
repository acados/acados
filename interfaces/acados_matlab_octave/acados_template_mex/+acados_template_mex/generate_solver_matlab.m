function generate_solver_matlab(acados_ocp_nlp_json_file)
    
    acados_template_folder = [getenv('ACADOS_INSTALL_DIR'), '/interfaces/acados_template/acados_template'];

    % get model name from json file
    acados_ocp = jsondecode(fileread(acados_ocp_nlp_json_file));

    model_name = acados_ocp.model.name;

    % setting up loader and environment
    template_glob = [acados_template_folder, '/c_templates_tera/*'];
    chdir('c_generated_code');
    % render source template
    template_file = 'main.in.c';
    out_file = ['main_', model_name, '.c'];
    % output file
    os_cmd = [getenv('ACADOS_INSTALL_DIR'), '/bin/t_renderer ', '"', template_glob, '"', ' ', '"', ...
    template_file, '"', ' ', '"', '../', acados_ocp_nlp_json_file, ...
    '"', ' ', '"', out_file, '"'];

    system(os_cmd);

    template_file = 'acados_solver.in.c';
    out_file = ['acados_solver_', model_name, '.c'];

    % output file
    os_cmd = [getenv('ACADOS_INSTALL_DIR'), '/bin/t_renderer ', '"', template_glob, '"', ' ', '"', ...
    template_file, '"', ' ', '"', '../', acados_ocp_nlp_json_file, ...
    '"', ' ', '"', out_file, '"'];

    system(os_cmd);

    % render solver template
    template_file = 'acados_solver.in.h';
    out_file = ['acados_solver_', model_name, '.h'];
    % output file
    os_cmd = [getenv('ACADOS_INSTALL_DIR'), '/bin/t_renderer ', '"', template_glob, '"', ' ', '"', ...
            template_file, '"', ' ', '"', '../', acados_ocp_nlp_json_file, ...
            '"', ' ', '"', out_file, '"'];

    system(os_cmd);

    % render header templates
    chdir([model_name, '_model/']);
    % render source template
    template_file = 'model.in.h';
    out_file = [model_name, '_model.h'];
    % output file
    os_cmd = [getenv('ACADOS_INSTALL_DIR'), '/bin/t_renderer ', '"', template_glob, '"', ' ', '"', ...
            template_file, '"', ' ', '"', '../../', acados_ocp_nlp_json_file, ...
            '"', ' ', '"', out_file, '"'];

    system(os_cmd);
    chdir('..');

    if(acados_ocp.dims.npd > 0)
        % render header templates
        dir_name = [acados_ocp.con_p.name, '_p_constraint/'];
        if~(exist(dir_name, 'dir'))
            mkdir(dir_name);
        end
        chdir(dir_name);
        % render source template
        template_file = 'p_constraint.in.h';
        out_file = [acados_ocp.con_p.name, '_p_constraint.h'];
        % output file
        os_cmd = [getenv('ACADOS_INSTALL_DIR'), '/bin/t_renderer ', '"', template_glob, '"', ' ', '"', ...
                template_file, '"', ' ', '"', '../../', acados_ocp_nlp_json_file, ...
                '"', ' ', '"', out_file, '"'];

        system(os_cmd);
        chdir('..');
    end

    if(acados_ocp.dims.nh > 0)
    dir_name = [acados_ocp.con_h.name, '_h_constraint/'];
        if~(exist(dir_name, 'dir'))
            mkdir(dir_name);
        end
        chdir(dir_name)
        % render source template
        template_file = 'h_constraint.in.h';
        out_file = [acados_ocp.con_h.name, '_h_constraint.h'];
        % output file
        os_cmd = [getenv('ACADOS_INSTALL_DIR'), '/bin/t_renderer ', '"', template_glob, '"', ' ', '"', ...
                template_file, '"', ' ', '"', '../../', acados_ocp_nlp_json_file, ...
                '"', ' ', '"', out_file, '"'];

        system(os_cmd);
        chdir('..');
    end

    % render Makefile template
    template_file = 'Makefile.in';
    out_file = 'Makefile';
    % output file
    os_cmd = [getenv('ACADOS_INSTALL_DIR'), '/bin/t_renderer ', '"', template_glob, '"', ' ', '"', ...
            template_file, '"', ' ', '"', '../', acados_ocp_nlp_json_file, ...
            '"', ' ', '"', out_file, '"'];

    system(os_cmd);

    % render source template
    template_file = 'acados_solver_sfun.in.c';
    out_file = ['acados_solver_sfunction_' , model_name, '.c'];
    % output file
    os_cmd = [getenv('ACADOS_INSTALL_DIR'), '/bin/t_renderer ', '"', template_glob, '"', ' ', '"', ...
            template_file, '"', ' ', '"', '../', acados_ocp_nlp_json_file, ...
            '"', ' ', '"', out_file, '"'];

    system(os_cmd);

    % render MATLAB make script
    template_file = 'make_sfun.in.m';
    out_file = 'acados_solver_sfun.in.c';
    % output file
    os_cmd = [getenv('ACADOS_INSTALL_DIR'), '/bin/t_renderer ', '"', template_glob, '"', ' ', '"', ...
            template_file, '"', ' ', '"', '../', acados_ocp_nlp_json_file, ...
            '"', ' ', '"', out_file, '"'];
    system(os_cmd);
   
    fprintf('Successfully generated acados solver!\n');

    % build generated code
    if isunix || ismac 
        % compile if on Mac and Unix platforms
        system('make');
        system('make shared_lib');
    else
        disp(['Commandline compilation of generated C code not yet supported under Windows.', ...
            'Please consider building the code in the c_generated_code folder from Windows Subsystem for Linux.'])
    end
    
    chdir('..');
    fprintf('Successfully built generated code!\n');

end
