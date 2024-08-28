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

function render_acados_sim_templates(acados_sim_json_file)

    acados_root_dir = getenv('ACADOS_INSTALL_DIR');
    acados_template_folder = fullfile(acados_root_dir,...
                          'interfaces', 'acados_template', 'acados_template');

    t_renderer_location = get_tera();

    %% load json data
    acados_sim = loadjson(fileread(acados_sim_json_file));
    model_name = acados_sim.model.name;

    %% render templates
    template_dir = fullfile(acados_template_folder, 'c_templates_tera','*');
    matlab_template_dir = fullfile(acados_template_folder, 'c_templates_tera','matlab_templates', '*');
    json_fullfile = fullfile(pwd, acados_sim_json_file);
    main_dir = pwd;
    chdir('c_generated_code');

    % cell array with entries (template_dir, template_file, output file)
    template_dir_file_and_output = { ...
        {template_dir, 'main_sim.in.c', ['main_sim_', model_name, '.c']}, ...
        {matlab_template_dir, 'mex_sim_solver.in.m', [model_name, '_mex_sim_solver.m']}, ...
        {matlab_template_dir, 'make_mex_sim.in.m', ['make_mex_sim_', model_name, '.m']}, ...
        {matlab_template_dir, 'acados_sim_create.in.c', ['acados_sim_create_', model_name, '.c']}, ...
        {matlab_template_dir, 'acados_sim_free.in.c', ['acados_sim_free_', model_name, '.c']}, ...
        {matlab_template_dir, 'acados_sim_set.in.c', ['acados_sim_set_', model_name, '.c']}, ...
        {template_dir, 'acados_sim_solver.in.c', ['acados_sim_solver_', model_name, '.c']}, ...
        {template_dir, 'acados_sim_solver.in.h', ['acados_sim_solver_', model_name, '.h']}, ...
        {matlab_template_dir, 'acados_sim_solver_sfun.in.c', ['acados_sim_solver_sfunction_', model_name, '.c']}, ...
        {matlab_template_dir, 'make_sfun_sim.in.m', ['make_sfun_sim_', model_name, '.m']}, ...
        {template_dir, 'Makefile.in','Makefile'}, ...
        {template_dir, 'CMakeLists.in.txt','CMakeLists.txt'}};

    num_entries = length(template_dir_file_and_output);
    for n=1:num_entries
        entry = template_dir_file_and_output{n};
        render_file(json_fullfile, entry{1}, entry{2}, entry{3}, t_renderer_location);
    end

    c_dir = pwd;
    chdir([model_name, '_model']);
    render_file(json_fullfile, template_dir, 'model.in.h', [model_name, '_model.h'], t_renderer_location);
    cd(c_dir);

    fprintf('Successfully rendered acados templates!\n');
    cd(main_dir)

end


%% auxilary functions
function render_file( json_fullfile, template_dir, template_file, out_file, ...
                      t_renderer_location )

    os_cmd = [t_renderer_location, ' "',...
        template_dir, '"', ' ', '"', template_file, '"', ' ', '"',...
        json_fullfile, '"', ' ', '"', out_file, '"'];

    [ status, result ] = system(os_cmd);
    if status
        cd ..
        error('rendering %s failed.\n command: %s\n returned status %d, got result:\n%s\n\n',...
            template_file, os_cmd, status, result);
    end
    % NOTE: this should return status != 0, maybe fix in tera renderer?
    if ~isempty(strfind( result, 'Error' )) % contains not implemented in Octave
        cd ..
        error('rendering %s failed.\n command: %s\n returned status %d, got result: %s',...
            template_file, os_cmd, status, result);
    end
end

