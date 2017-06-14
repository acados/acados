add_test(NAME ocp_qp_matlab
    COMMAND ${MATLAB_EXECUTABLE} -noFigureWindows -nodesktop -nosplash -r
        "addpath('${EXAMPLES_SRC_DIR}/matlab');try, ocp_qp_example, \
        catch ME, disp(ME.getReport()), exit(1), end, exit(0)")

add_test(NAME ocp_nlp_matlab
    COMMAND ${MATLAB_EXECUTABLE} -noFigureWindows -nodesktop -nosplash -r
        "addpath('${EXAMPLES_SRC_DIR}/matlab');try, ocp_nlp_example, \
        catch ME, disp(ME.getReport()), exit(1), end, exit(0)")