%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% The 2-Clause BSD License
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.;

%

clearvars; clc; close all;


% list the examples you would like to test
targets = {
    '../p_global_example/test_p_global_without_dependencies.m';
    '../p_global_example/main.m';
    '../p_global_example/simulink_test_p_global.m';
    '../mocp_transition_example/main_parametric_mocp.m';
    '../pendulum_on_cart_model/nonlinear_constraint_test.m';
    '../dense_nlp/test_qpscaling.m';
};


pass = zeros(1, length(targets));  % keep track of test results
messages = cell(1, length(targets));  % and error messages
setenv("TEST_DIR", pwd)
for idx = 1:length(targets)
    setenv("TEST_MESSAGE", "")
    [dir, file, extension] = fileparts(targets{idx});

    testpath = getenv("TEST_DIR");
    % clear variables that might contain CasADi objects to avoid SWIG
    % warnings
    clear ocp_solver ocp ocp_model model sim sim_model sim_solver params
    save(strcat(testpath, "/test_workspace.mat"))
    setenv("LD_RUN_PATH", strcat(testpath, "/", dir, "/c_generated_code"))

    try
        run(targets{idx});
        test_val = true;
    catch exception
        setenv("TEST_MESSAGE", exception.message)
        warning(exception.message);
        clear exception
        test_val = false;
    end

    % use absolute path, since current directory depends on point of failure
    testpath = getenv("TEST_DIR");
    load(strcat(testpath, "/test_workspace.mat"));
    disp(['test', targets{idx},' success'])
    messages{idx} = getenv("TEST_MESSAGE");
    if contains(targets{idx},'simulink'); bdclose('all'); end
    delete(strcat(testpath, "/test_workspace.mat"));
    % delete generated code to avoid failure in examples using similar names
    code_gen_dir = strcat(testpath, "/", dir, "/c_generated_code");
    % if exist(code_gen_dir, 'dir')
    %     rmdir(code_gen_dir, 's')
    % end
    close all;
    % clc;
end

% clc;
fail = false;
disp('Succesful tests: ')
for idx = 1:length(targets)
    if strcmp(messages{idx},"")
        disp(targets{idx})
    end
end
disp(' ')
disp('Failed tests: ')
for idx = 1:length(targets)
    if ~strcmp(messages{idx},"")
        disp(targets{idx})
        disp(['message: ',messages{idx}])
        fail = true;
    end
end
if fail==true
    error('Test failure');
end
clearvars