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
    '../generic_dyn_disc/disc_dyn_example_ocp.m';
    '../generic_external_cost/external_cost_example_ocp.m';
    '../getting_started/extensive_example_ocp.m';
%     '../getting_started/minimal_example_closed_loop.m';
    '../getting_started/minimal_example_sim.m';
%     '../getting_started/simulink_example.m';
%     '../getting_started/simulink_example_advanced.m';
%     '../linear_mass_spring_model/example_closed_loop.m';
    '../linear_mass_spring_model/example_ocp.m';
%     '../linear_mass_spring_model/example_sim.m';
    '../linear_mpc/example_closed_loop.m';
    '../lorentz/example_mhe.m';
%     '../masses_chain_model/example_closed_loop.m';
    '../masses_chain_model/example_ocp.m';
    '../pendulum_dae/example_closed_loop.m';
    '../pendulum_dae/example_sim.m';
%     '../pendulum_on_cart_model/example_closed_loop.m';
    '../pendulum_on_cart_model/example_ocp.m';
%     '../pendulum_on_cart_model/example_ocp_custom_hess.m';
%     '../pendulum_on_cart_model/example_ocp_param_sens.m';
%     '../pendulum_on_cart_model/example_ocp_reg.m';
%     '../pendulum_on_cart_model/example_sim.m';
%     '../pendulum_on_cart_model/example_solution_sens_closed_loop.m';
%     '../pendulum_on_cart_model/experiment_dae_formulation.m';
    '../race_cars/main.m';
    '../simple_dae_model/example_ocp.m';
%     '../swarming/example_closed_loop.m';
    '../swarming/example_ocp.m';
%     '../swarming/example_sim.m';
%     '../wind_turbine_nx6/example_closed_loop.m';
    '../wind_turbine_nx6/example_ocp.m';
%     '../wind_turbine_nx6/example_sim.m';
%
%     './test_checks.m';
%     './test_mhe_lorentz.m';
%     './test_ocp_OSQP.m';
%     './test_ocp_linear_mass_spring.m';
%     './test_ocp_pendulum_dae.m';
%     './test_ocp_pendulum_on_cart.m';
%     './test_ocp_qpdunes.m';
%     './test_ocp_simple_dae.m';
%     './test_ocp_wtnx6.m';
%     './test_sens_adj.m';
%     './test_sens_forw.m';
%     './test_sens_hess.m';
%     './test_sim_dae.m';
%     './test_target_selector.m'
};

pass = zeros(1, length(targets));  % keep track of test results
messages = cell(1, length(targets));  % and error messages
setenv("TEST_DIR", pwd)
for idx = 1:length(targets)
    setenv("TEST_MESSAGE", "")
    [dir, file, extension] = fileparts(targets{idx});

    testpath = getenv("TEST_DIR");
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
    rmdir(strcat(testpath, "/", dir, "/c_generated_code"), 's')
    close all; clc;
end

clc;
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