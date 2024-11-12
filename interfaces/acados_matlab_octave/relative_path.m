function rel_path = relative_path(target_path, base_path)
    % RELATIVE_PATH - Computes the relative path from base_path to target_path.
    % Works across multiple operating systems (Windows, macOS, Linux) and is
    % compatible with both MATLAB and Octave.

    % Convert both paths to absolute paths
    target_path = absolute_path(target_path);
    base_path = absolute_path(base_path);

    % return '.' if paths match
    if strcmp(target_path, base_path)
        rel_path = '.';
        return;
    end

    % Split the paths into their individual components
    target_parts = strsplit(target_path, filesep);
    base_parts = strsplit(base_path, filesep);

    % Find the first non-matching part of the paths
    i = 1;
    while i <= min(length(target_parts), length(base_parts)) && strcmp(target_parts{i}, base_parts{i})
        i = i + 1;
    end

    % The number of '..' to go up from the base_path to the common root
    up_steps = length(base_parts) - i + 1;

    % Calculate the relative path
    rel_parts = [repmat({'..'}, 1, up_steps), target_parts(i:end)];

    % Join the parts into the final relative path
    if isempty(rel_parts)
        rel_path = '';
    else
        rel_path = fullfile(rel_parts{:});
    end
end


function abs_path = absolute_path(path)
    % ABSOLUTE_PATH - Returns the absolute path for a given input path.
    % Ensures compatibility across platforms and works with both Octave and MATLAB.

    if ispc
        % Windows platform - Normalize slashes
        path = strrep(path, '/', '\');
    else
        % Unix-like systems - Normalize backslashes
        path = strrep(path, '\', '/');
    end

    % Ensure absolute path
    if ~is_absolute(path)
        path = fullfile(pwd, path);
    end

    % Remove redundant parts (., ..)
    abs_path = remove_dot_segments(path);
end

function is_abs = is_absolute(path)
    % IS_ABSOLUTE - Checks if a given path is absolute.
    if ispc
        % On Windows, absolute paths start with a drive letter (e.g., C:\)
        is_abs = length(path) >= 2 && path(2) == ':';
    else
        % On Unix-like systems, absolute paths start with '/'
        is_abs = length(path) >= 1 && path(1) == '/';
    end
end

function clean_path = remove_dot_segments(path)
    % REMOVE_DOT_SEGMENTS - Removes redundant . and .. from a path.

    parts = strsplit(path, filesep);
    new_parts = {};

    for i = 1:length(parts)
        part = parts{i};
        if strcmp(part, '..')
            if ~isempty(new_parts)
                new_parts(end) = [];  % Go up one directory
            end
        elseif ~strcmp(part, '.') && ~isempty(part)
            new_parts{end+1} = part; % Add non-empty, non-current parts
        end
    end

    % Join the parts back together
    clean_path = fullfile(new_parts{:});
end
