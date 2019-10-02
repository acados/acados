function generate_solver(acados_ocp_nlp_json_file, python_interpreter_path)

    acados_template_mex_path = [getenv('ACADOS_INSTALL_DIR') '/interfaces/acados_matlab/acados_template_mex/+acados_template_mex' ] ;
    commandStr = [python_interpreter_path, ' ', acados_template_mex_path, '/generate_solver_matlab.py --json_file_name ', acados_ocp_nlp_json_file];
    [status, commandOut] = system(commandStr);
    disp(commandOut);
    if status
        fprintf('error calling generate_solver script!');
    end
end