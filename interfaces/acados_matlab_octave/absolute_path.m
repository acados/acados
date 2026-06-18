

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

    % Detect if the original path was an absolute Unix path (starts with '/')
    leading_sep = false;
    if ~ispc && ~isempty(path) && path(1) == filesep
        leading_sep = true;
    end

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

    % Join the parts back together and restore leading separator for
    % absolute Unix paths. If no parts remain and it was absolute, return '/'.
    if isempty(new_parts)
        if leading_sep
            clean_path = filesep;
        else
            clean_path = '';
        end
    else
        clean_path = fullfile(new_parts{:});
        if leading_sep
            % prepend leading slash (use char concatenation - fullfile would
            % otherwise strip a leading empty component)
            clean_path = [filesep clean_path];
        end
    end
end
