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
    matlab_template_path = 'matlab_templates';
    json_fullfile = fullfile(pwd, acados_sim_json_file);
    main_dir = pwd;
    chdir('c_generated_code');

    % cell array with entries (template_file, output file)
    template_list = { ...
        {'main_sim.in.c', ['main_sim_', model_name, '.c']}, ...
        {fullfile(matlab_template_path, 'mex_sim_solver.in.m'), [model_name, '_mex_sim_solver.m']}, ...
        {fullfile(matlab_template_path, 'make_mex_sim.in.m'), ['make_mex_sim_', model_name, '.m']}, ...
        {fullfile(matlab_template_path, 'acados_sim_create.in.c'), ['acados_sim_create_', model_name, '.c']}, ...
        {fullfile(matlab_template_path, 'acados_sim_free.in.c'), ['acados_sim_free_', model_name, '.c']}, ...
        {fullfile(matlab_template_path, 'acados_sim_set.in.c'), ['acados_sim_set_', model_name, '.c']}, ...
        {'acados_sim_solver.in.c', ['acados_sim_solver_', model_name, '.c']}, ...
        {'acados_sim_solver.in.h', ['acados_sim_solver_', model_name, '.h']}, ...
        {fullfile(matlab_template_path, 'acados_sim_solver_sfun.in.c'), ['acados_sim_solver_sfunction_', model_name, '.c']}, ...
        {fullfile(matlab_template_path, 'make_sfun_sim.in.m'), ['make_sfun_sim_', model_name, '.m']}, ...
        {'Makefile.in', 'Makefile'}, ...
        {'CMakeLists.in.txt', 'CMakeLists.txt'}};

    num_entries = length(template_list);
    for n=1:num_entries
        entry = template_list{n};
        render_file( entry{1}, entry{2}, json_fullfile);
    end

    c_dir = pwd;
    chdir([model_name, '_model']);
    render_file( 'model.in.h', [model_name, '_model.h'], json_fullfile);
    cd(c_dir);

    fprintf('Successfully rendered acados templates!\n');
    cd(main_dir)
end
