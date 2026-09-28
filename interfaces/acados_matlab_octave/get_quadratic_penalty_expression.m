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

function penalty = get_quadratic_penalty_expression(h_expr, lh, uh, Z_l, Z_u)
    % Returns a CasADi expression corresponding to a quadratic penalty on the
    % constraint violation with quadratic weight diag(Z_l) for lower bound
    % violations and diag(Z_u) for upper bound violations.
    %
    % Inputs:
    %   h_expr : CasADi SX/MX expression of size n x 1 or 1 x n
    %   lh     : lower bound vector of length n
    %   uh     : upper bound vector of length n
    %   Z_l    : lower-bound quadratic penalty weights, length n
    %   Z_u    : upper-bound quadratic penalty weights, length n

    if ~(isa(h_expr, 'casadi.SX') || isa(h_expr, 'casadi.MX'))
        error('h_expr must be a CasADi SX or MX expression.');
    end

    h_expr = h_expr(:);
    lh = lh(:);
    uh = uh(:);
    Z_l = Z_l(:);
    Z_u = Z_u(:);

    if ~(numel(h_expr) == numel(lh) && numel(h_expr) == numel(uh) ...
            && numel(h_expr) == numel(Z_l) && numel(h_expr) == numel(Z_u))
        error('lh, uh, Z_l, and Z_u must all have the same length as h_expr.');
    end

    lower_violation = fmax(lh - h_expr, 0);
    upper_violation = fmax(h_expr - uh, 0);

    penalty = 0.5 * sum(Z_l .* lower_violation.^2 + Z_u .* upper_violation.^2);
end
