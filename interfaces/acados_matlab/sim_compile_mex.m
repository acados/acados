function sim_compile_mex(opts)

% get acados folder
acados_folder = getenv('ACADOS_INSTALL_DIR');
mex_flags = getenv('ACADOS_MEX_FLAGS');

% set paths
acados_mex_folder = fullfile(acados_folder, 'interfaces', 'acados_matlab');
acados_include = ['-I' acados_folder];
acados_interfaces_include = ['-I' fullfile(acados_folder, 'interfaces')];
acados_lib_path = ['-L' fullfile(acados_folder, 'lib')];

mex_names = { ...
	'sim_create', ...
	'sim_destroy', ...
	'sim_create_ext_fun', ...
	'sim_destroy_ext_fun', ...
	'sim_solve', ...
	'sim_set', ...
	'sim_get', ...
	'sim_precompute' ...
};

mex_files = cell(length(mex_names), 1);
for k=1:length(mex_names)
	mex_files{k} = fullfile(acados_mex_folder, [mex_names{k}, '.c']);
end

% compile mex
if is_octave()
	if exist(fullfile(opts.output_dir, 'cflags_octave.txt'), 'file')==0
		diary(fullfile(opts.output_dir, 'cflags_octave.txt'));
		diary on
		mkoctfile -p CFLAGS
		diary off
		input_file = fopen(fullfile(opts.output_dir, 'cflags_octave.txt'), 'r');
		cflags_tmp = fscanf(input_file, '%[^\n]s');
		fclose(input_file);
		cflags_tmp = [cflags_tmp, ' -std=c99 -fopenmp'];
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
		mex(acados_include, acados_interfaces_include, acados_lib_path, '-lacados', '-lhpipm', '-lblasfeo', mex_files{ii})
	else
		mex(mex_flags, 'CFLAGS=$CFLAGS -std=c99 -fopenmp', acados_include, acados_interfaces_include, acados_lib_path, '-lacados', '-lhpipm', '-lblasfeo', mex_files{ii})
	end
end


if is_octave()
	movefile('*.o', opts.output_dir);
end

%system(['mv -f *.mexa64 ', opts.output_dir])
for k=1:length(mex_names)
	clear(mex_names{k})
%	[status, message] = movefile([mex_names{k}, '.', mexext], opts.output_dir)
	[status, message] = copyfile([mex_names{k}, '.', mexext], opts.output_dir);
	delete([mex_names{k}, '.', mexext]);
end
