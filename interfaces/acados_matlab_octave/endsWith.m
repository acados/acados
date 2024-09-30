function result = endsWith(in, suffix)
    % Check if the input is a string
    if ischar(in)
        % Check if the string ends with the suffix
        if length(in) >= length(suffix) && strcmp(in(end-length(suffix)+1:end), suffix)
            result = true;
        else
            result = false;
        end
    elseif iscell(in)
        % Initialize the result array
        result = false(size(in));

        % Iterate over each string in the array
        for i = 1:numel(in)
            str = in{i};

            % Check if the string ends with the suffix
            if length(str) >= length(suffix) && strcmp(str(end-length(suffix)+1:end), suffix)
                result(i) = true;
            end
        end
    else
        error('Input must be a string or a cell array of strings.');
    end
end
