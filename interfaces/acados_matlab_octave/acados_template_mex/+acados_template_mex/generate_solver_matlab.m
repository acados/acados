function generate_solver_matlab(acados_ocp_nlp_json_file)

    acados_template_folder = [getenv('ACADOS_INSTALL_DIR'), '/interfaces/acados_template/acados_template'];

    if ~isfile([getenv('ACADOS_INSTALL_DIR'), '/bin/t_renderer'])
        error('acados/bin/t_renderer not found! Please download the t_renderer binaries from https://github.com/acados/tera_renderer/releases and place them in acados/bin (strip version and extension from the binary name)')
    end
    % get model name from json file
    acados_ocp = jsondecode(fileread(acados_ocp_nlp_json_file));

    model_name = acados_ocp.model.name;

    % setting up loader and environment
    template_dir = [acados_template_folder, '/c_templates_tera/*'];
    chdir('c_generated_code');

    %% render source template
    template_file = 'main.in.c';
    out_file = ['main_', model_name, '.c'];
    render_file( acados_ocp_nlp_json_file, template_dir, template_file, out_file )

    % render solver template
    template_file = 'acados_solver.in.c';
    out_file = ['acados_solver_', model_name, '.c'];
    render_file( acados_ocp_nlp_json_file, template_dir, template_file, out_file )

    template_file = 'acados_solver.in.h';
    out_file = ['acados_solver_', model_name, '.h'];
    render_file( acados_ocp_nlp_json_file, template_dir, template_file, out_file )

    % render sim_solver template
    template_file = 'acados_sim_solver.in.c';
    out_file = ['acados_sim_solver_', model_name, '.c'];
    render_file( acados_ocp_nlp_json_file, template_dir, template_file, out_file )

    template_file = 'acados_sim_solver.in.h';
    out_file = ['acados_sim_solver_', model_name, '.h'];
    render_file( acados_ocp_nlp_json_file, template_dir, template_file, out_file )

    % render header templates
    chdir([model_name, '_model/']);
    template_file = 'model.in.h';
    out_file = [model_name, '_model.h'];
    render_file( acados_ocp_nlp_json_file, template_dir, template_file, out_file, 1)

    chdir('..');

    if (acados_ocp.dims.npd > 0)
        % render header templates
        dir_name = [acados_ocp.con_p.name, '_p_constraint/'];
        if~(exist(dir_name, 'dir'))
            mkdir(dir_name);
        end
        chdir(dir_name);
        % render source template
        template_file = 'p_constraint.in.h';
        out_file = [acados_ocp.con_p.name, '_p_constraint.h'];
        render_file( acados_ocp_nlp_json_file, template_dir, template_file, out_file, 1)

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
        render_file( acados_ocp_nlp_json_file, template_dir, template_file, out_file, 1)

        chdir('..');
    end

    % render Makefile template
    template_file = 'Makefile.in';
    out_file = 'Makefile';
    render_file( acados_ocp_nlp_json_file, template_dir, template_file, out_file )

    % render source template
    template_file = 'acados_solver_sfun.in.c';
    out_file = ['acados_solver_sfunction_' , model_name, '.c'];
    render_file( acados_ocp_nlp_json_file, template_dir, template_file, out_file )

    % render MATLAB make script
    template_file = 'make_sfun.in.m';
    out_file = 'acados_solver_sfun.in.c';
    render_file( acados_ocp_nlp_json_file, template_dir, template_file, out_file )
   
    fprintf('Successfully generated acados solver!\n');

    % build generated code
    if isunix || ismac 
        % compile if on Mac or Unix platform
        [ status, result ] = system('make');
        if status
            cd ..
            error('building templated code failed.\nGot status %d, result: %s',...
                  status, result);
        end
        [ status, result ] = system('make shared_lib');
        if status
            cd ..
            error('building templated code as shared library failed.\nGot status %d, result: %s',...
                  status, result);
        end
    else
        disp(['Commandline compilation of generated C code not yet supported under Windows.', ...
            'Please consider building the code in the c_generated_code folder from Windows Subsystem for Linux.'])
    end
    
    chdir('..');
    fprintf('Successfully built generated code!\n');

end


%% auxilary function

function render_file( acados_ocp_nlp_json_file, template_dir, template_file, out_file, subfolder_depth )
    if nargin < 5
        subfolder_depth = 0;
    end
    json_location = repmat('../', 1, subfolder_depth+1);

    os_cmd = [getenv('ACADOS_INSTALL_DIR'), '/bin/t_renderer ', '"',...
        template_dir, '"', ' ', '"', template_file, '"', ' ', '"',...
        json_location, acados_ocp_nlp_json_file, '"', ' ', '"', out_file, '"'];
    
    [ status, result ] = system(os_cmd);
    if status
        cd ..
        error('rendering %s failed.\n command: %s\n returned status %d, got result:\n%s\n\n',...
            template_file, os_cmd, status, result);
%     else
%         disp(['Redering ' template_file ': success!']);
    end
    % this should return status != 0, maybe fix in tera renderer?
    if contains( result, 'Error' )
        cd ..
        error('rendering %s failed.\n command: %s\n returned status %d, got result: %s',...
            template_file, os_cmd, status, result);
    end
end
