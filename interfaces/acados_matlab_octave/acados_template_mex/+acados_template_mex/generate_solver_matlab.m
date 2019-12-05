function generate_solver_matlab(acados_ocp_nlp_json_file)

    acados_template_folder = fullfile(getenv('ACADOS_INSTALL_DIR'),...
                          'interfaces', 'acados_template', 'acados_template');

    %% check if t_renderer is available -> download if not
    if ispc()
        t_renderer_location = fullfile(getenv('ACADOS_INSTALL_DIR'),'bin','t_renderer.exe');
    else
        t_renderer_location = fullfile(getenv('ACADOS_INSTALL_DIR'),'bin','t_renderer');
    end
    
    if ~isfile( t_renderer_location )

        message = ['\nDear acados user, we could not find t_renderer binaries,',...
            '\n which are needed to export templated C code from ',...
            'Matlab.\n Press any key to proceed setting up the t_renderer automatically.',...
            '\n Press "n" or "N" to exit, if you wish to set up t_renderer yourself.\n',...
            '\n https://github.com/acados/tera_renderer/releases'];

        In = input(message,'s');

        if strcmpi( In, 'n')
            error('Please set up t_renderer yourself and try again');
        else
            t_renderer_version = 'v0.0.20';
            if isunix()
                suffix = '-linux';
            elseif ismac()
                suffix = '-osx';
            elseif ispc()
                suffix = '-windows.exe';
            end

            tera_url = ['https://github.com/acados/tera_renderer/releases/download/', ...
                    t_renderer_version '/t_renderer-', t_renderer_version, suffix];
            destination = fullfile(getenv('ACADOS_INSTALL_DIR'), 'bin');
            tmp_file = websave(destination, tera_url);
            
            acados_bin_dir = fullfile(getenv('ACADOS_INSTALL_DIR'),'bin');
            if ~exist(acados_bin_dir, 'dir')
                [~,~] = mkdir(acados_bin_dir);
            end
            
            movefile(tmp_file, t_renderer_location);

            if isunix() || ismac()
                % make executable
                system(['chmod a+x ', t_renderer_location]);
            end
            fprintf('\nSuccessfully set up t_renderer\n')
        end
    end

    % get model name from json file
    acados_ocp = jsondecode(fileread(acados_ocp_nlp_json_file));

    model_name = acados_ocp.model.name;

    % setting up loader and environment
    template_dir = fullfile(acados_template_folder, 'c_templates_tera','*');
    chdir('c_generated_code');

    %% render source template
    json_location = fullfile('..');
    
    % main
    template_file = 'main.in.c';
    out_file = ['main_', model_name, '.c'];
    render_file( acados_ocp_nlp_json_file, template_dir, template_file, out_file, t_renderer_location, json_location )

    % make_main_mex
    template_file = 'make_main_mex.in.m';
    out_file = ['make_main_mex_', model_name, '.m'];
    render_file( acados_ocp_nlp_json_file, template_dir, template_file, out_file, t_renderer_location, json_location )

    % main for matlab/octave
    template_file = 'main_mex.in.c';
    out_file = ['main_mex_', model_name, '.c'];
    render_file( acados_ocp_nlp_json_file, template_dir, template_file, out_file, t_renderer_location, json_location )

    % make_mex
    template_file = 'make_mex.in.m';
    out_file = ['make_mex_', model_name, '.m'];
    render_file( acados_ocp_nlp_json_file, template_dir, template_file, out_file, t_renderer_location, json_location )

    % MEX constructor
    template_file = 'acados_mex_create.in.c';
    out_file = ['acados_mex_create_', model_name, '.c'];
    render_file( acados_ocp_nlp_json_file, template_dir, template_file, out_file, t_renderer_location, json_location )

    % MEX destructor
    template_file = 'acados_mex_free.in.c';
    out_file = ['acados_mex_free_', model_name, '.c'];
    render_file( acados_ocp_nlp_json_file, template_dir, template_file, out_file, t_renderer_location, json_location )

    % MEX solve
    template_file = 'acados_mex_solve.in.c';
    out_file = ['acados_mex_solve_', model_name, '.c'];
    render_file( acados_ocp_nlp_json_file, template_dir, template_file, out_file, t_renderer_location, json_location )

    % render solver template
    template_file = 'acados_solver.in.c';
    out_file = ['acados_solver_', model_name, '.c'];
    render_file( acados_ocp_nlp_json_file, template_dir, template_file, out_file, t_renderer_location, json_location  )

    template_file = 'acados_solver.in.h';
    out_file = ['acados_solver_', model_name, '.h'];
    render_file( acados_ocp_nlp_json_file, template_dir, template_file, out_file, t_renderer_location, json_location  )

    % render sim_solver template
    template_file = 'acados_sim_solver.in.c';
    out_file = ['acados_sim_solver_', model_name, '.c'];
    render_file( acados_ocp_nlp_json_file, template_dir, template_file, out_file, t_renderer_location, json_location  )

    template_file = 'acados_sim_solver.in.h';
    out_file = ['acados_sim_solver_', model_name, '.h'];
    render_file( acados_ocp_nlp_json_file, template_dir, template_file, out_file, t_renderer_location, json_location  )

    % render header templates
    chdir([model_name, '_model']);
    json_location = fullfile('..', '..');
    
    template_file = 'model.in.h';
    out_file = [model_name, '_model.h'];
    render_file( acados_ocp_nlp_json_file, template_dir, template_file, out_file, t_renderer_location, json_location )

    chdir('..');

    if (acados_ocp.dims.npd > 0)
        % render header templates
        dir_name = [acados_ocp.con_p.name, '_p_constraint'];
        if~(exist(dir_name, 'dir'))
            mkdir(dir_name);
        end
        chdir(dir_name);
        % render source template
        template_file = 'p_constraint.in.h';
        out_file = [acados_ocp.con_p.name, '_p_constraint.h'];
        json_location = fullfile('..', '..');
        render_file( acados_ocp_nlp_json_file, template_dir, template_file, out_file, t_renderer_location, json_location)

        chdir('..');
    end

    if(acados_ocp.dims.nh > 0)
        dir_name = [acados_ocp.con_h.name, '_h_constraint'];
        if~(exist(dir_name, 'dir'))
            mkdir(dir_name);
        end
        chdir(dir_name)
        % render source template
        template_file = 'h_constraint.in.h';
        out_file = [acados_ocp.con_h.name, '_h_constraint.h'];
        json_location = fullfile('..', '..');
        render_file( acados_ocp_nlp_json_file, template_dir, template_file, out_file, t_renderer_location, json_location)

        chdir('..');
    end
    
    json_location = fullfile('..');

    % render Makefile template
    template_file = 'Makefile.in';
    out_file = 'Makefile';
    render_file( acados_ocp_nlp_json_file, template_dir, template_file, out_file, t_renderer_location, json_location )

    % render source template
    template_file = 'acados_solver_sfun.in.c';
    out_file = ['acados_solver_sfunction_' , model_name, '.c'];
    render_file( acados_ocp_nlp_json_file, template_dir, template_file, out_file, t_renderer_location, json_location )

    % render MATLAB make script
    template_file = 'make_sfun.in.m';
    out_file = 'make_sfun.m';
    render_file( acados_ocp_nlp_json_file, template_dir, template_file, out_file, t_renderer_location, json_location )
   
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

function render_file( acados_ocp_nlp_json_file, template_dir, template_file, out_file, ...
                      t_renderer_location, json_location )

    os_cmd = [t_renderer_location, ' "',...
        template_dir, '"', ' ', '"', template_file, '"', ' ', '"',...
        fullfile(json_location, acados_ocp_nlp_json_file), '"', ' ', '"', out_file, '"'];
    
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
