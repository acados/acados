function set_up_t_renderer(t_renderer_location,force)
arguments
    t_renderer_location
    force=true
end
    message = ['\nDear acados user, we could not find t_renderer binaries,',...
        '\n which are needed to export templated C code from ',...
        'Matlab.\n Press any key to proceed setting up the t_renderer automatically.',...
        '\n Press "n" or "N" to exit, if you wish to set up t_renderer yourself.\n',...
        '\n https://github.com/acados/tera_renderer/releases'];
    
    if(~force)
        In = input(message,'s');
        if strcmpi( In, 'n')
            error('Please set up t_renderer yourself and try again');
        end
    end
	t_renderer_version = 'v0.0.34';
	if ismac()
		suffix = '-osx';
	elseif isunix()
		suffix = '-linux';
	elseif ispc()
		suffix = '-windows';
	end
	acados_root_dir = getenv('ACADOS_INSTALL_DIR');

	tera_url = ['https://github.com/acados/tera_renderer/releases/download/', ...
			t_renderer_version '/t_renderer-', t_renderer_version, suffix];
	destination = fullfile(acados_root_dir, 'bin');
	tmp_file = websave(destination, tera_url);

	if ~exist(destination, 'dir')
		[~,~] = mkdir(destination);
	end

	movefile(tmp_file, t_renderer_location);

	if isunix()
		% make executable
		system(['chmod a+x ', t_renderer_location]);
	end
	fprintf('\nSuccessfully set up t_renderer\n')
end