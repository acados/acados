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

clear; clc; close all;
p = path;  % save current path for loading later
start_dir = pwd;

% list the examples you would like to test
targets = {
%     '../generic_dyn_disc/disc_dyn_example_ocp.m';
%     '../generic_external_cost/external_cost_example_ocp.m';
%     '../getting_started/extensive_example_ocp.m';
%     '../getting_started/minimal_example_closed_loop.m';
%     '../getting_started/minimal_example_sim.m';
%     '../getting_started/simulink_example.m';
%     '../getting_started/simulink_example_advanced.m';
%     '../linear_mass_spring_model/example_closed_loop.m';
%     '../linear_mass_spring_model/example_ocp.m';
%     '../linear_mass_spring_model/example_sim.m';
%     '../linear_mpc/example_closed_loop.m';
%     '../lorentz/example_mhe.m';
%     '../masses_chain_model/example_closed_loop.m';
%     '../masses_chain_model/example_ocp.m';
%     '../pendulum_dae/example_closed_loop.m';  % error while reshaping cost.Vz_0
%     '../pendulum_dae/example_sim.m';
%     '../pendulum_on_cart_model/example_closed_loop.m';
%     '../pendulum_on_cart_model/example_ocp.m';
%     '../pendulum_on_cart_model/example_ocp_custom_hess.m';
%     '../pendulum_on_cart_model/example_ocp_param_sens.m';  % Unable to resolve the name ocp.t_ocp.eval_param_sens.
%     '../pendulum_on_cart_model/example_ocp_reg.m';
%     '../pendulum_on_cart_model/example_sim.m';
%     '../pendulum_on_cart_model/example_solution_sens_closed_loop.m';  % Unable to resolve the name ocp.t_ocp.eval_param_sens.
%     '../pendulum_on_cart_model/experiment_dae_formulation.m';
%     '../race_cars/main.m';  % acados returns status 4
%     '../simple_dae_model/example_ocp.m';
%     '../swarming/example_closed_loop.m';
%     '../swarming/example_ocp.m';
%     '../swarming/example_sim.m';
%     '../wind_turbine_nx6/example_closed_loop.m';
%     '../wind_turbine_nx6/example_ocp.m';
%     '../wind_turbine_nx6/example_sim.m';

%     './test_checks.m';  % SHOULD FAIL, DOES
%     './test_mhe_lorentz.m';  % OK
%     './test_ocp_OSQP.m';  % CRASH - solver not installed
%     './test_ocp_linear_mass_spring.m';  % OK 
%     './test_ocp_pendulum_dae.m';  % OK
%     './test_ocp_pendulum_on_cart.m';  % CRASH - solver not installed
%     './test_ocp_qpdunes.m';  % CRASH - solver not installed
%     './test_ocp_simple_dae.m';  % OK
%     './test_ocp_wtnx6.m';  % FAIL - when not in installation folder?
%     './test_sens_adj.m';  % OK
%     './test_sens_forw.m';  % OK
%     './test_sens_hess.m';  % OK
%     './test_sim_dae.m';  % OK
%     './test_target_selector.m'  % OK
};

test_ok = zeros(1,length(targets));  % keep track of test results
test_messages = cell(1,length(targets));  % and error messages

for test_number = 1:length(targets)
    [target_dir,target_file,target_extension] = fileparts(targets{test_number});
    cd(target_dir)
    save test_workspace.mat  % due to "clear all" in the examples
    try
        run([target_file,target_extension])
        clear; close all; clc;  % clean up after a succesful test
        load test_workspace.mat
        if contains(target_file,'simulink'); bdclose('all'); end
        test_ok(test_number) = 1;  % save the result
        disp(['testing ',targets{test_number},' succesful'])
        pause(3)  % wait a bit before moving on
    catch exception
        clearvars -except exception;  % keep the caught error
        close all; clc;
        load test_workspace.mat
        if contains(target_file,'simulink'); bdclose('all'); end
        warning(['testing ',targets{test_number},' failed'])
        test_messages{test_number} = exception.message;  % keep the message
        clear exception;  % forget error before the next test
        pause(3)  % wait a bit before moving on
    end
    delete test_workspace.mat  % clean up the example folder
    cd(start_dir)  % go back to the starting folder
    path(p);  % restore the saved path
end

% print out test results
clc;
disp('Succesful tests: ')
for test_number = 1:length(targets)
    if test_ok(test_number)
        disp(targets{test_number})
    end
end
disp(' ')
disp('Failed tests: ')
for test_number = 1:length(targets)
    if ~test_ok(test_number)
        disp(targets{test_number})
        disp(['    message: ',test_messages{test_number}])
    end
end