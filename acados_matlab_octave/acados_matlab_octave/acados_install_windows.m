function acados_install_windows(varargin)
% acados_install_windows([CmakeConfigString])
% Install script for acados on windows
% CmakeConfigString - config string for the CMAKE command [Default='-DBUILD_SHARED_LIBS=OFF -DACADOS_WITH_OSQP=OFF']


    switch(nargin)
        case 0
            cmakeConfigString='-DBUILD_SHARED_LIBS=OFF -DACADOS_WITH_OSQP=OFF';
        case 1
            cmakeConfigString=varargin{1};
        otherwise 
            error('function called with %d parameters, was expecting max 1',nargin);
    end    

    % Derive the path for the acados root
    fullPath = mfilename('fullpath');
    % Extract path for this file
    [folderPath,~,~]=fileparts(fullPath);
    % Remove the two top directories to get the acados path
    [folderPath,~,~]=fileparts(folderPath);
    [acadosPath,~,~]=fileparts(folderPath);

    acadosBuildPath=fullfile(acadosPath,'build');

    %% Acados installation instructions:
    % Add the subfolders bin and x86_64-w64-mingw32\bin of the above mentioned mingw installation to your environment variable PATH.
    fprintf('Setting up environment PATH variable\n');

    % Store the current PATH environment variable value
    origEnvPath=getenv('PATH');
    % Read the mex C compiler configuration and extract the location
    cCompilerConfig = mex.getCompilerConfigurations('C');
    pathToCompilerLocation = cCompilerConfig.Location;
    % Modify the environment PATH variable for this Matlab session such that
    % the mex C compiler takes priority ensuring calls to gcc uses the
    % configured mex compiler
    setenv('PATH', [fullfile(pathToCompilerLocation,'bin') ';' origEnvPath]);

    %% Acados installation instructions:
    % $ACADOS_INSTALL_DIR=$(pwd)
    % mkdir -p build
    % cd build
    fprintf('Creating build dir in %s\n', acadosBuildPath);
    setenv('ACADOS_INSTALL_DIR', acadosPath);
    mkdir(acadosBuildPath);
    cd(acadosBuildPath);

    %% Acados installation instructions:
    % cmake.exe -G "MinGW Makefiles" -DACADOS_INSTALL_DIR="$ACADOS_INSTALL_DIR" -DBUILD_SHARED_LIBS=OFF -DACADOS_WITH_OSQP=ON ..
    % # useful options to add above:
    % # -DACADOS_WITH_QPOASES=ON/OFF -DACADOS_WITH_OSQP=ON/OFF -DACADOS_WITH_QPDUNES=ON/OFF ..
    % # -DBLASFEO_TARGET=GENERIC -DHPIPM_TARGET=GENERIC
    % # NOTE: check the output of cmake: -- Installation directory: should be <acados_root_folder>,
    % #     if this is not the case, set -DACADOS_INSTALL_DIR=<acados_root_folder> explicitly above.
    fprintf('Executing cmake configuration\n');
    % Command slightly modified to work with CMD instead of PowerShell
    cmake_cmd = sprintf('cmake.exe -G "MinGW Makefiles" -DACADOS_INSTALL_DIR=%%ACADOS_INSTALL_DIR%% %s ..',cmakeConfigString);

    status=system(cmake_cmd);
    if (status~=0)
        error('cmake command failed. Command was:\n%s\n', cmake_cmd);
    end

    %% Acados installation instructions:
    % mingw32-make.exe -j4
    % mingw32-make.exe install
    fprintf('Compiling and installing using minGW\n');

    compile_command='mingw32-make.exe -j4';
    status = system(compile_command);
    if status ~= 0
        %Store the PATH environment variable used during compile for error reporting
        envPath=getenv('PATH');
        setenv('PATH',origEnvPath);
        error('Compilation of acados failed %s\n%s\n%s\n%s\n\n', ...
            'Compile command was:', compile_command, ...
            'Environment path was: ', envPath);
    end

    install_command='mingw32-make.exe install';
    status = system(install_command);
    if status ~= 0
        %Store the PATH environment variable used during compile for error reporting
        envPath=getenv('PATH');
        setenv('PATH',origEnvPath);
        error('Installation of acados failed %s\n%s\n%s\n%s\n\n', ...
            'Install command was:', install_command, ...
            'Environment path was: ', envPath);
    end

    % Restore path to original
    setenv('PATH',origEnvPath);

    %% Download external dependencies
    % casadi
    fprintf('Downloading casadi\n');
    addpath(fullfile(acadosPath,'interfaces', 'acados_matlab_octave'));
    run('acados_env_variables_windows');
    check_acados_requirements(true);

    %% renderer
    fprintf('Install the template renderer\n');
    acados_root_dir = getenv('ACADOS_INSTALL_DIR');
    %% check if t_renderer is available -> download if not
    t_renderer_location = fullfile(acados_root_dir, 'bin', 't_renderer.exe');

    if ~exist( t_renderer_location, 'file' )
        acados_template_mex.set_up_t_renderer( t_renderer_location )
    end

end
