%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.


function render_file( in_file, out_file, json_fullfile, template_glob )

    t_renderer_location = get_tera();

    acados_root_dir = getenv('ACADOS_INSTALL_DIR');
    if nargin < 4
        acados_template_folder = fullfile(acados_root_dir,...
                            'interfaces', 'acados_template', 'acados_template', 'c_templates_tera');
        [path, name, ext] = fileparts(in_file);
        template_glob = fullfile(acados_template_folder, path, '**', '*');
        in_file = [name, ext];
    end

    os_cmd = [t_renderer_location, ' "',...
        template_glob, '"', ' ', '"', in_file, '"', ' ', '"',...
        json_fullfile, '"', ' ', '"', out_file, '"'];

    [ status, result ] = system(os_cmd);
    if status
        cd ..
        error('rendering %s failed.\n command: %s\n returned status %d, got result:\n%s\n\n',...
            in_file, os_cmd, status, result);
    end
    % NOTE: this should return status != 0, maybe fix in tera renderer?
    if ~isempty(strfind( result, 'Error' )) % contains not implemented in Octave
        cd ..
        error('rendering %s failed.\n command: %s\n returned status %d, got result: %s',...
            in_file, os_cmd, status, result);
    end
end

