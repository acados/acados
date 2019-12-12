function compile_main()
    cd c_generated_code
    %% build main file
    if isunix || ismac 
        % compile if on Mac or Unix platform
        [ status, result ] = system('make');
        if status
            cd ..
            error('building templated code failed.\nGot status %d, result: %s',...
                  status, result);
        end
        [ status, result ] = system('make shared_lib');
        if status
            cd ..
            error('building templated code as shared library failed.\nGot status %d, result: %s',...
                  status, result);
        end
    else
        disp(['Commandline compilation of generated C code not yet supported under Windows.', ...
            'Please consider building the code in the c_generated_code folder from Windows Subsystem for Linux.'])
    end
    cd ..
    fprintf('Successfully built main file!\n');
end