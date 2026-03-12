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

function [] = verify_weighting_matrix(A, name, tol)
    % verify_weighting_matrix - Check if a matrix is square, symmetric, and
    % either positive definite or (diagonal and positive semi-definite)
    %and raise an error if not.
    %
    % Parameters:
    %   A   - square matrix to check
    %   name - name of the matrix for error message
    %   tol - tolerance for eigenvalue comparison (default: 1e-10)
    %
    if nargin < 3
        tol = 1e-10;
    end

    if ~ismatrix(A) || size(A, 1) ~= size(A, 2)
        error('Matrix %s is not square.', name);
    end
    if norm(A - A.', inf) > tol
        error('Matrix %s is not symmetric.', name);
    end

    if isequal(A, diag(diag(A)))
        if any(diag(A) < 0)
            error('Diagonal weighting matrix %s is not positive semi-definite.', name);
        end
    else
        E = eig(A);
        result = all(E > tol);
        if ~result
            error('Matrix %s is not positive definite. Eigenvalues: %s', name, mat2str(E));
        end
    end
end
