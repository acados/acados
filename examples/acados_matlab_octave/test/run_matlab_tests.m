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


assert(1+1==2)

disp('assertation works')

disp('checking environment variables')

disp('MATLABPATH')
disp(getenv('MATLABPATH'))

disp('MODEL_FOLDER')
disp(getenv('MODEL_FOLDER'))


disp('ENV_RUN')
disp(getenv('ENV_RUN'))

disp('LD_LIBRARY_PATH')
disp(getenv('LD_LIBRARY_PATH'))

disp('pwd')
disp(pwd)

disp('running tests')

%% run all tests
test_names = [
    "test_code_reuse",
    "test_sim_code_reuse",
    "run_test_dim_check",
"run_test_ocp_mass_spring",
% "run_test_ocp_pendulum",
"run_test_ocp_wtnx6",
% "run_test_sim_adj",
"run_test_sim_dae",
% "run_test_sim_forw",
"run_test_sim_hess",
"param_test",
];

for k = 1:length(test_names)
    disp(strcat("running test ", test_names(k)));
    run(test_names(k))
end
