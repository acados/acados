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



function [p_global, m, l, coefficients, coefficient_vals, knots, p_global_values] = create_p_global(lut)

    import casadi.*
    m = MX.sym('m');
    l = MX.sym('l');
    p_global = {m, l};
    p_global_values = [0.1; 0.8];

    large_scale = false;
    if lut
        % generate random values for spline coefficients
        % knots = {[0,0,0,0,0.2,0.5,0.8,1,1,1,1],[0,0,0,0.1,0.5,0.9,1,1,1]};

        if large_scale
            % large scale lookup table
            knots = {0:200,0:200};
            coefficient_vals = 0.1*ones(38809, 1);
        else
            % small scale lookup table
            knots = {0:19,0:19};
            coefficient_vals = 0.1*ones(256, 1);
        end

        coefficients = MX.sym('coefficient', numel(coefficient_vals), 1);
        p_global{end+1} = coefficients;
        p_global_values = [p_global_values; coefficient_vals(:)];
    else
        coefficient_vals = [];
        knots = [];
        coefficients = MX.sym('coefficient', 0, 1);
    end

    p_global = vertcat(p_global{:});
end
