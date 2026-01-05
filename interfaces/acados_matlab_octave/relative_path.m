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
