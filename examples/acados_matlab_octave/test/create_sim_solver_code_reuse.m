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
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 'AS IS'
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


function sim_solver = create_sim_solver_code_reuse(creation_mode)

    addpath('../pendulum_on_cart_model')

    check_acados_requirements()

    json_file = 'pendulum_ocp.json';
    solver_creation_opts = struct();
    solver_creation_opts.json_file = json_file;
    if strcmp(creation_mode, 'standard')
        disp('Standard creation mode');
    elseif strcmp(creation_mode, 'precompiled') || strcmp(creation_mode, 'no_sim')
        solver_creation_opts.generate = false;
        solver_creation_opts.build = false;
        solver_creation_opts.compile_mex_wrapper = false;
    else
        error('Invalid creation mode')
    end

    if strcmp(creation_mode, 'no_sim')
        sim = [];
    else
        model = get_pendulum_on_cart_model();
        sim = AcadosSim();
        sim.model = model;
        sim.solver_options.Tsim = 0.1; % simulation time
        sim.solver_options.integrator_type = 'ERK';
    end

    %% create integrator
    sim_solver = AcadosSimSolver(sim, solver_creation_opts);
end
