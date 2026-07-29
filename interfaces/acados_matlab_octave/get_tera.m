function t_renderer_location = get_tera()
    acados_root_dir = getenv('ACADOS_INSTALL_DIR');

    % check if t_renderer is available -> download if not
    if ispc()
        t_renderer_location = fullfile(acados_root_dir, 'bin', 't_renderer.exe');
    else
        t_renderer_location = fullfile(acados_root_dir, 'bin', 't_renderer');
    end

    if ~exist( t_renderer_location, 'file' )
        set_up_t_renderer( t_renderer_location )
    elseif ~is_tera_version_sufficient(t_renderer_location)
        fprintf('\nt_renderer found but its version does not meet the minimum requirement. Updating automatically...\n');
        set_up_t_renderer( t_renderer_location )
    end
end


function sufficient = is_tera_version_sufficient(t_renderer_location)
% is_tera_version_sufficient - Check whether an installed t_renderer binary
%   meets the minimum version requirement.
%
%   sufficient = is_tera_version_sufficient(t_renderer_location)
%
%   Parameters:
%     t_renderer_location - full path to the t_renderer executable
%
%   Returns:
%     sufficient - true if the binary is present and its version is >=
%                  TERA_MINIMUM_VERSION, false otherwise.
%
%   Versions older than v0.2.0 do not support --version; a non-zero exit
%   status from that call is therefore treated as "version too old".

    TERA_MINIMUM_VERSION = [0, 2, 0];

    % Versions older than v0.2.0 do not support --version at all.
    % A non-zero exit status therefore means the version is too old.
    [status, version_output] = system(['"' t_renderer_location '" --version']);
    if status ~= 0
        sufficient = false;
        return
    end

    % Expected output format: "t_renderer X.Y.Z" or just "X.Y.Z"
    tokens = regexp(strtrim(version_output), '(\d+)\.(\d+)\.(\d+)', 'tokens');
    if isempty(tokens)
        % Cannot parse version - assume it needs updating
        sufficient = false;
        return
    end
    version_parts = cellfun(@str2double, tokens{1});
    sufficient = compare_versions(version_parts, TERA_MINIMUM_VERSION) >= 0;
end


function result = compare_versions(v1, v2)
% compare_versions - Lexicographically compare two version vectors.
%
%   result = compare_versions(v1, v2)
%
%   Parameters:
%     v1, v2 - numeric row vectors of version components, e.g. [0, 2, 0]
%
%   Returns:
%     result -  1 if v1 > v2
%               0 if v1 == v2
%              -1 if v1 < v2
    for i = 1:min(length(v1), length(v2))
        if v1(i) > v2(i)
            result = 1; return;
        elseif v1(i) < v2(i)
            result = -1; return;
        end
    end
    result = sign(length(v1) - length(v2));
end
