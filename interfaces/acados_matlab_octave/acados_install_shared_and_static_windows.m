function acados_install_shared_and_static_windows(varargin)
% Script for installing acados as shared and static libs on windows.
% All shared libs will be moved to bin, static libs are in `lib`.
% CmakeConfigString - config string for the CMAKE command [Default='-DCMAKE_POLICY_VERSION_MINIMUM=3.5']
    switch(nargin)
        case 0
            cmakeConfigString='-DCMAKE_POLICY_VERSION_MINIMUM=3.5';
    end
    if contains(cmakeConfigString, '-DBUILD_SHARED_LIBS')
        error('options BUILD_SHARED_LIBS should not be contained in cmake config.')
    end

    disp('installing acados with shared libs')
    acados_install_windows('-DBUILD_SHARED_LIBS=ON -DCMAKE_POLICY_VERSION_MINIMUM=3.5')
    acados_root_dir = getenv('ACADOS_INSTALL_DIR');
    acados_bin_path = fullfile(acados_root_dir, "bin");
    acados_lib_path = fullfile(acados_root_dir, "lib");
    [status, message ] = movefile(fullfile(acados_lib_path, '*.lib'), acados_bin_path);

    disp('installing acados with static libs')
    acados_install_windows('-DBUILD_SHARED_LIBS=OFF -DCMAKE_POLICY_VERSION_MINIMUM=3.5')
end
