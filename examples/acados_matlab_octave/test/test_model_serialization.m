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
% Test serialization of CasADi expressions in AcadosModel
% This test verifies that CasADi expressions can be serialized and deserialized
% correctly, which is necessary for JSON dumping and loading.
%

function test_model_serialization()
    import casadi.*

    disp('Testing CasADi expression serialization in AcadosModel...');

    % Create a simple pendulum model
    model = AcadosModel();
    model.name = 'test_pendulum';

    % Define symbolic variables
    x = SX.sym('x', 4);
    xdot = SX.sym('xdot', 4);
    u = SX.sym('u', 1);

    model.x = x;
    model.xdot = xdot;
    model.u = u;

    % Add some dynamics
    model.f_expl_expr = vertcat(x(3), x(4), u(1), -9.81);

    % Add a cost expression
    model.cost_y_expr = vertcat(x, u);

    % Add a constraint expression
    model.con_h_expr = x(1)^2 + x(2)^2;

    % Test serialize method
    disp('  Testing serialize method...');
    [serialized_str, expression_names] = model.serialize();
    
    % Check that we got some expressions
    assert(~isempty(serialized_str), 'Serialized string should not be empty');
    assert(~isempty(expression_names), 'Expression names should not be empty');
    
    % Check that expected expressions are in the list
    expected_exprs = {'x', 'xdot', 'u', 'f_expl_expr', 'cost_y_expr', 'con_h_expr'};
    for i = 1:length(expected_exprs)
        assert(any(strcmp(expression_names, expected_exprs{i})), ...
            ['Expected expression ' expected_exprs{i} ' not found in expression_names']);
    end
    
    disp('  Serialize method: PASSED');

    % Test deserialize method
    disp('  Testing deserialize method...');
    model2 = AcadosModel();
    model2.deserialize(serialized_str, expression_names);
    
    % Check that expressions were restored
    assert(~isempty(model2.x), 'Deserialized x should not be empty');
    assert(~isempty(model2.u), 'Deserialized u should not be empty');
    assert(~isempty(model2.f_expl_expr), 'Deserialized f_expl_expr should not be empty');
    
    % Check dimensions match
    assert(length(model2.x) == length(model.x), 'x dimensions should match');
    assert(length(model2.u) == length(model.u), 'u dimensions should match');
    
    disp('  Deserialize method: PASSED');

    % Test convert_to_struct_for_json_dump
    disp('  Testing convert_to_struct_for_json_dump...');
    model_struct = model.convert_to_struct_for_json_dump();
    
    % Check that struct contains serialized expressions
    assert(isfield(model_struct, 'serialized_expressions'), ...
        'Struct should contain serialized_expressions field');
    assert(isfield(model_struct, 'expression_names'), ...
        'Struct should contain expression_names field');
    assert(~isempty(model_struct.serialized_expressions), ...
        'serialized_expressions should not be empty');
    assert(~isempty(model_struct.expression_names), ...
        'expression_names should not be empty');
    
    disp('  convert_to_struct_for_json_dump: PASSED');

    % Test from_struct static method
    disp('  Testing from_struct static method...');
    model3 = AcadosModel.from_struct(model_struct);
    
    % Check that expressions were restored
    assert(~isempty(model3.x), 'Restored x should not be empty');
    assert(~isempty(model3.u), 'Restored u should not be empty');
    assert(~isempty(model3.f_expl_expr), 'Restored f_expl_expr should not be empty');
    
    % Check dimensions match
    assert(length(model3.x) == length(model.x), 'Restored x dimensions should match');
    assert(length(model3.u) == length(model.u), 'Restored u dimensions should match');
    
    % Check name was preserved
    assert(strcmp(model3.name, model.name), 'Model name should be preserved');
    
    disp('  from_struct static method: PASSED');

    disp('All serialization tests PASSED!');
end
