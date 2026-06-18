function h = hash_struct(struct)

    % remove keys that should not affect the hash
    ignored_fields = {'external_function_files_model', 'external_function_files_ocp', 'json_loaded'};
    for i = 1:length(ignored_fields)
        if isfield(struct, ignored_fields{i})
            struct = rmfield(struct, ignored_fields{i});
        end
    end
    str = savejson(struct);
    if is_octave()
        h = hash('MD2', str);
    else
        % try
        %     % Since R2022b
        %     h = keyHash(str);
        % catch
            try
                % Should work in older versions
                import java.security.*;
                import java.math.*;
                md = MessageDigest.getInstance('MD5');
                hash = md.digest(double(str));
                bi = BigInteger(1, hash);
                h = char(bi.toString(16));
            catch
                warning("Could not hash struct, returning empty string");
                h = '';
            end
        % end
    end
end
