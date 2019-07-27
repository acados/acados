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

mex_names = { ...
  'ocp_create', ...
  'ocp_destroy', ...
  'ocp_create_ext_fun', ...
  'ocp_destroy_ext_fun', ...
  'ocp_solve', ...
  'ocp_precompute', ...
  'ocp_set', ...
  'ocp_get' ...
};
mex_files = cell(length(mex_names), 1);
for k=1:length(mex_names)
  mex_files{k} = fullfile(acados_mex_folder, [mex_names{k}, '.c']);
end

% compile mex
if is_octave()
	if exist(fullfile(opts.output_dir, 'cflags_octave.txt'), 'file')==0
		diary(fullfile(opts.output_dir, 'cflags_octave.txt'))
		diary on
		mkoctfile -p CFLAGS
		diary off
		input_file = fopen(fullfile(opts.output_dir, 'cflags_octave.txt'), 'r');
		cflags_tmp = fscanf(input_file, '%[^\n]s');
		fclose(input_file);
		cflags_tmp = [cflags_tmp, ' -std=c99 -fopenmp'];
		if (strcmp(opts.qp_solver, 'full_condensing_qpoases'))
			cflags_tmp = [cflags_tmp, ' -DACADOS_WITH_QPOASES'];
		end
		input_file = fopen(fullfile(opts.output_dir, 'cflags_octave.txt'), 'w');
		fprintf(input_file, '%s', cflags_tmp);
		fclose(input_file);
	end
	input_file = fopen(fullfile(opts.output_dir, 'cflags_octave.txt'), 'r');
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
  movefile('*.o', opts.output_dir)
end

for k=1:length(mex_names)
  clear(mex_names{k})
  movefile([mex_names{k}, '.', mexext], opts.output_dir);
end




