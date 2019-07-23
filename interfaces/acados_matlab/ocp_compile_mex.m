function ocp_compile_mex(opts)

% get acados folder
acados_folder = getenv('ACADOS_INSTALL_DIR');
mex_flags = getenv('ACADOS_MEX_FLAGS');

% set paths
acados_mex_folder = fullfile(acados_folder, 'interfaces', 'acados_matlab');
acados_include = ['-I', acados_folder];
acados_interfaces_include = ['-I', fullfile(acados_folder, 'interfaces')];
external_include = ['-I', fullfile(acados_folder, 'external')];
blasfeo_include = ['-I', fullfile(acados_folder, 'external', 'blasfeo', 'include')];
acados_lib_path = ['-L', fullfile(acados_folder, 'lib')];

% compile mex
mex_files = { ...
	fullfile(acados_mex_folder, 'ocp_create.c'), ...
	fullfile(acados_mex_folder, 'ocp_destroy.c'), ...
	fullfile(acados_mex_folder, 'ocp_create_ext_fun.c'), ...
	fullfile(acados_mex_folder, 'ocp_destroy_ext_fun.c'), ...
	fullfile(acados_mex_folder, 'ocp_solve.c'), ...
	fullfile(acados_mex_folder, 'ocp_precompute.c'), ...
	fullfile(acados_mex_folder, 'ocp_set.c'), ...
	fullfile(acados_mex_folder, 'ocp_get.c'), ...
	} ;

if is_octave()
	if exist('build/cflags_octave.txt')==0
		diary 'build/cflags_octave.txt'
		diary on
		mkoctfile -p CFLAGS
		diary off
		input_file = fopen('build/cflags_octave.txt', 'r');
		cflags_tmp = fscanf(input_file, '%[^\n]s');
		fclose(input_file);
		cflags_tmp = [cflags_tmp, ' -std=c99 -fopenmp'];
		if (strcmp(opts.qp_solver, 'full_condensing_qpoases'))
			cflags_tmp = [cflags_tmp, ' -DACADOS_WITH_QPOASES'];
		end
		input_file = fopen('build/cflags_octave.txt', 'w');
		fprintf(input_file, '%s', cflags_tmp);
		fclose(input_file);
	end
	input_file = fopen('build/cflags_octave.txt', 'r');
	cflags_tmp = fscanf(input_file, '%[^\n]s');
	fclose(input_file);
	setenv('CFLAGS', cflags_tmp);
end

for ii=1:length(mex_files)
	disp(['compiling ', mex_files{ii}])
	if is_octave()
%		mkoctfile -p CFLAGS
		if (strcmp(opts.qp_solver, 'full_condensing_qpoases'))
			mex(acados_include, acados_interfaces_include, external_include, blasfeo_include, acados_lib_path, '-lacados', '-lhpipm', '-lblasfeo', '-lqpOASES_e', mex_files{ii})
		else
			mex(acados_include, acados_interfaces_include, external_include, blasfeo_include, acados_lib_path, '-lacados', '-lhpipm', '-lblasfeo', mex_files{ii})
		end
	else
		if (strcmp(opts.qp_solver, 'full_condensing_qpoases'))
			mex(mex_flags, 'CFLAGS=$CFLAGS -std=c99 -fopenmp -DACADOS_WITH_QPOASES', acados_include, acados_interfaces_include, external_include, blasfeo_include, acados_lib_path, '-lacados', '-lhpipm', '-lblasfeo', '-lqpOASES_e', mex_files{ii})
		else
			mex(mex_flags, 'CFLAGS=$CFLAGS -std=c99 -fopenmp', acados_include, acados_interfaces_include, external_include, blasfeo_include, acados_lib_path, '-lacados', '-lhpipm', '-lblasfeo', mex_files{ii})
		end
	end
end

if is_octave()
	system(['mv -f *.o build/']);
	system(['mv -f *.mex build/']);
else
	system(['mv -f *.mexa64 build/']);
end
