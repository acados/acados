%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

function [] = verify_weighting_matrix(A, name, tol)
    % verify_weighting_matrix - Check if a matrix is square, symmetric, and
    % either positive definite or (diagonal and positive semidefinite)
    % and raise an error if not.
    %
    % Parameters:
    %   A   - square matrix to check
    %   name - name of the matrix for error message
    %   tol - tolerance for eigenvalue comparison (default: 1e-10)
    %
    if nargin < 3
        tol = 1e-10;
    end

    % allow empty matrix, corresponding to ny = 0 (mirrors Python interface)
    if isempty(A)
        if size(A, 1) ~= size(A, 2)
            error('Weighting matrix %s is not square.', name);
        end
        return
    end

    if ~ismatrix(A) || size(A, 1) ~= size(A, 2)
        error('Matrix %s is not square.', name);
    end
    if norm(A - A.', inf) > tol
        warning('Matrix %s is not symmetric.', name);
    else
        if isequal(A, diag(diag(A)))
            if any(diag(A) < 0)
                warning('Diagonal weighting matrix %s is not positive semidefinite.', name);
            end
        else
            try
                E = eig(A);
                result = all(E > tol);
            catch
                warning('Eigenvalue decomposition of %s failed, matrix might not be positive definite.', name)
                result = true;
            end
            if ~result
                warning('Matrix %s is not positive definite. Eigenvalues: %s', name, mat2str(E));
            end
        end
    end
end
