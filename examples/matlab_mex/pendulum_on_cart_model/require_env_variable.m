function require_env_variable(var_name)
    try
        tmp = getenv(var_name);
        if isempty(tmp)
            error(strcat('ERROR: variable ', var_name, ' not defined'));
        end
    catch error
        exit_with_error(error);
    end
end