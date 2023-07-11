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
%   Author: Jonathan Frey: jonathanpaulfrey(at)gmail.com

function [ gnsf ] = reformulate_with_invertible_E_mat( gnsf, model, print_info)
%% Description
% this function checks that the necessary condition to apply the gnsf
% structure exploiting integrator to a model, namely that the matrices E11,
% E22 are invertible holds.
% If this is not the case, it will make these matrices invertible and add
% corresponding terms, to the term C * phi, such that the obtained model is
% still equivalent

%% import CasADi and load models

import casadi.*

% check invertibility of E11, E22; and reformulate if needed
ind_11 = 1:gnsf.nx1;
ind_22 = gnsf.nx1+1 : gnsf.nx1+gnsf.nz1;

if print_info
    disp(' ');
    disp('----------------------------------------------------');
    disp('checking rank of E11 and E22');
    disp('----------------------------------------------------');
end

%% check if E11, E22 are invertible

if or( rank(gnsf.E( ind_11, ind_11)) ~= gnsf.nx1, ...
        rank(gnsf.E( ind_22, ind_22)) ~= gnsf.nz1 ) 
    
    % print warning (always)
    disp(['the rank of E11 or E22 is not full after the reformulation']);
    disp(' ');
    disp(['the script will try to reformulate the model with an invertible matrix instead']);
    disp(['NOTE: this feature is based on a heuristic, it should be used with care!!!']);
    
    %% load models
    xdot = gnsf.xdot;
    z = gnsf.z;
    
    % % GNSF
    % get dimensions
    nx1 = gnsf.nx1;    
    x1dot = xdot(1:nx1);

    k = [x1dot; z];
    for i = [1,2]
        if i == 1
            ind = 1:gnsf.nx1;
        else
            ind = gnsf.nx1+1 : gnsf.nx1 + gnsf.nz1;
        end
        mat = gnsf.E(ind, ind);
        while rank(mat) < length(ind)
            if print_info
                disp(' ');
                disp(['the rank of E',num2str(i),num2str(i),' is not full']);
                disp(['the algorithm will try to reformulate the model with an invertible matrix instead']);
                disp(['NOTE: this feature is not super stable and might need more testing!!!!!!']);
            end

            for sub_max = ind
                sub_ind = min(ind):sub_max; 
                % regard the submatrix mat(sub_ind, sub_ind);
                sub_mat = gnsf.E(sub_ind, sub_ind);
                if rank(sub_mat) < length(sub_ind)
                    % reformulate the model by adding a 1 to last diagonal
                    % element and changing rhs respectively.
                    gnsf.E(sub_max, sub_max) = gnsf.E(sub_max, sub_max) + 1;
                    % this means adding the term 1 * k(sub_max) to the sub_max
                    % row of the l.h.s
                    if isempty(find(gnsf.C(sub_max,:), 1))
                        % add new nonlinearity entry
                        gnsf.C(sub_max, gnsf.n_out + 1) = 1;
                        gnsf.phi_expr = [gnsf.phi_expr; k(sub_max)];
                    else
                        ind_f = find(gnsf.C(sub_max,:));
                        % add term to corresponding nonlinearity entry
                        % note: herbey we assume that C is a selection matrix,
                        % i.e. gnsf.phi_expr(ind_f) is only entering one equation;
                        if length(find(gnsf.C(:,ind_f))) ~= 1
                            error('matrix C is not a selection matrix, reformulation with invertible E11, E22 not supported!!!');
                        else
                            gnsf.phi_expr(ind_f) = gnsf.phi_expr(ind_f) + k(sub_max) / ...
                                gnsf.C(sub_max, ind_f);
                        end
                    end
                end
            end
            gnsf = determine_input_nonlinearity_function( gnsf );
        end
    end
    check_reformulation(model, gnsf, print_info);
    disp('successfully reformulated the model with invertible matrices E11, E22');
else
    if print_info
        disp(' ');
        disp('the rank of both E11 and E22 is naturally full after the reformulation ');
        disp('==>  model reformulation finished');
        disp(' ');
    end
end

if det(gnsf.E_LO) == 0
    disp('_______________________________________________________________________________________________________');
    disp(' ');
    disp('TAKE CARE ');
    disp('E_LO matrix is NOT regular after automatic transcription!');
    disp('->> this means the model CANNOT be used with the gnsf integrator');    
    disp('->> it probably means that one entry (of xdot or z) that was moved to the linear output type system');
    disp('    does not appear in the model at all (zero column in E_LO)');
    disp(' OR: the columns of E_LO are linearly dependent ');
    disp(' ');
    disp(' SOLUTIONs: a) go through your model & check equations the method wanted to move to LOS');
    disp('            b) deactivate the detect_LOS option');
    disp('_______________________________________________________________________________________________________');
end

end
