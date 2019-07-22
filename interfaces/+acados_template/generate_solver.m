function generate_solver(acados_ocp_nlp_json_file, python_interpreter_path)
    import acados_template.*
    s = what('acados_template');
    acados_template_path = s.path;
    commandStr = [python_interpreter_path, ' ', acados_template_path, '/generate_solver_matlab.py --json_file_name ', acados_ocp_nlp_json_file];
    [status, commandOut] = system(commandStr);
    disp(commandOut);
    if status
        fprintf('error calling generate_solver script!');
    end
end