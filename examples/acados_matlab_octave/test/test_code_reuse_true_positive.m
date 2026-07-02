% clear generated code
folderName = 'c_generated_code';

[status, msg, msgID] = rmdir(folderName, 's');

if status
    fprintf('Folder %s and contents deleted.', folderName);
else
    fprintf('Failed: %s\n', msg);
end


ocp = create_pendulum_ocp();
% create solver

solver_creation_opts = struct();
% precompiled mode.
solver_creation_opts.generate = false;
solver_creation_opts.build = false;
solver_creation_opts.compile_mex_wrapper = false;

% test 1) check creation is done when solver is not there
ocp_solver = AcadosOcpSolver(ocp, solver_creation_opts);
ocp_solver.solve();

if (~ocp_solver.solver_creation_opts.generate || ~ocp_solver.solver_creation_opts.build || ~ocp_solver.solver_creation_opts.compile_mex_wrapper)
    error('Should have detected that code reuse is not possible');
end
clear ocp_solver;

% test 2) check reuse works
ocp_solver = AcadosOcpSolver(ocp, solver_creation_opts);
if (ocp_solver.solver_creation_opts.generate || ocp_solver.solver_creation_opts.build || ocp_solver.solver_creation_opts.compile_mex_wrapper)
    error('Code reuse should have worked.');
end


% test 3) check creation is done when OCPs are different
ocp.cost.yref = 1.2 * ocp.cost.yref;
text_output = evalc('ocp_solver = AcadosOcpSolver(ocp, solver_creation_opts);')
disp(text_output);
if (~ocp_solver.solver_creation_opts.generate || ~ocp_solver.solver_creation_opts.build || ~ocp_solver.solver_creation_opts.compile_mex_wrapper)
    error('Should have detected that code reuse is not possible');
end


% isempty(strfind(..) <=> ~contains() with Octave compatibility.
if isempty(strfind(text_output, 'yref'))
    error('should complain about yref');
end

if isempty(strfind(text_output, 'reuse not possible'))
    error('should say reuse not possible');
end
disp('Code reuse not possible was detected successfully.')
