function t_renderer_location = get_tera()
    acados_root_dir = getenv('ACADOS_INSTALL_DIR');

    % check if t_renderer is available -> download if not
    if ispc()
        t_renderer_location = fullfile(acados_root_dir, 'bin', 't_renderer.exe');
    else
        t_renderer_location = fullfile(acados_root_dir, 'bin', 't_renderer');
    end

    if ~exist( t_renderer_location, 'file' )
        set_up_t_renderer( t_renderer_location )
    end
end
