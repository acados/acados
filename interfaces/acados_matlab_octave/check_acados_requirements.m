function check_acados_requirements(varargin)
% check_acados_requirements([force])
% If force is not provided it is default set to true

    %The first variable parameter determines whether to force install
    switch(nargin)
        case 0
            force=true;
        case 1
            force=varargin{1};
        otherwise
            error('function called with %d parameters, was expecting max 1',nargin);       
    end


    % check environment variables
    env_run = getenv('ENV_RUN');
    if (~strcmp(env_run, 'true'))
        error('env.sh has not been sourced! Before executing this example, run: source env.sh');
    end
    acados_dir = getenv('ACADOS_INSTALL_DIR');
    
    if ~is_casadi_available
        % offer to install CasADi
        message = ['\nDear acados user, we could not import CasADi',...
            ',\n which is needed to run the acados Matlab/Octave examples.',...
            '\n Press any key to proceed setting up the CasADi automatically.',...
            '\n Press "n" or "N" to exit, if you wish to set up CasADi yourself.\n'];
        if ~force
            In = input(message,'s');
        else
            In = 'Y';
        end
        
        if strcmpi( In, 'n') 
            error('Please set up CasADi yourself and try again.');
        else
            % download CasADi
            CasADi_version = '3.5.5';
            url = strcat('https://github.com/casadi/casadi/releases/download/',...
                   CasADi_version, '/');
            external_folder = fullfile(acados_dir, 'external');
            if ~is_octave
                destination = fullfile(external_folder, 'casadi-matlab');
            else
                filename = strcat('casadi-linux-octave-5.2.0-v', CasADi_version, '.tar.gz');
                destination = fullfile(external_folder, 'casadi-octave');
                fprintf(['\nWe cannot set up CasADi for Octave.\n',...
                        'Please download\n', [url,filename],...
                        '\nand unpack it into <acados_root_dir>/external/casadi-octave.\n\n']);
                error('exiting');
            end


            if ismac
                if ~verLessThan('matlab', '8.5')
                    filename = strcat('casadi-osx-matlabR2015a-v', CasADi_version, '.tar.gz');
                elseif ~verLessThan('matlab', '8.4')
                    filename = strcat('casadi-osx-matlabR2014b-v', CasADi_version, '.tar.gz');
                elseif ~verLessThan('matlab', '8.3')
                    filename = strcat('casadi-osx-matlabR2014a-v', CasADi_version, '.tar.gz');
                end
            elseif isunix
                if verLessThan('matlab', '8.4')
                    filename = strcat('casadi-linux-matlabR2014a-v', CasADi_version, '.tar.gz');
                else % R2014b or later
                    filename = strcat('casadi-linux-matlabR2014b-v', CasADi_version, '.tar.gz');
                end

            elseif ispc
                if ~verLessThan('matlab', '9.0')
                    filename = strcat('casadi-windows-matlabR2016a-v', CasADi_version,'.zip');
                elseif ~verLessThan('matlab', '8.4')
                    filename = strcat('casadi-windows-matlabR2014b-v', CasADi_version,'.zip');
                elseif ~verLessThan('matlab', '8.3')
                    filename = strcat('casadi-windows-matlabR2014a-v', CasADi_version,'.zip');
                end

            end

            try
                disp(['trying to download and unpack: ', url,filename])
                file = websave(destination, [url,filename]);
                [~,~,ending] = fileparts(file);

                % unpack
                if strcmp(ending, '.zip')
                    unzip(file, destination)
                else
                    untar(file, destination)
                end

                addpath(destination)
            catch
                error(['Sorry, we could not set up CasADi for your system, please try manually\n.',...
                    'Instructions can be found on https://web.casadi.org/get/\n',...
                    'We recommend using CasADi version', CasADi_version]);
            end
        end
        if ~is_casadi_available
                error(['Sorry, we could not set up CasADi for your system, please try manually\n.',...
                    'Instructions can be found on https://web.casadi.org/get/\n',...
                    'We recommend using CasADi version', CasADi_version]);
        end
    end
end

function got_casadi = is_casadi_available()
    got_casadi = 1;
    try
    % check CasADi availibility
        import casadi.*
        test = SX.sym('test');
    catch e
        got_casadi = 0;
    end
end
